/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/// \file

#ifndef yap__AttributeUtilities_h
#define yap__AttributeUtilities_h

#include "fwd/AttributeUtilities.h"

#include "container_utils.h"
#include "Exceptions.h"

#include <memory>
#include <type_traits>
#include <vector>

namespace yap {

/// base class implementing a return_type typedef
template <typename T>
struct with_return_type
{
    /// \typedef return_type of functors
    using return_type = T;
};

/// functor class to return an attribute of the template parameter type
/// \tparam T return type
/// \tparam U argument type
/// \defgroup Attributes attribute classes
template <typename T, typename U>
struct attribute_of<T, U> : with_return_type<T>
{
    /// U& functor
    virtual T operator()(const U&) const = 0;

    /// U* functor
    virtual T operator()(const U* const u) const
    {
        if (!u)
            throw exceptions::Exception("null ptr", "attribute_of::operator()");
        return operator()(*u);
    }

    /// shared_ptr<U> functor
    virtual T operator()(const std::shared_ptr<const U>& u) const
    { return operator()(u.get()); }
};

/// functor class to return an attribute of any of the template parameter types
/// \tparam T return type
/// \tparam U argument type
/// \tparam V argument types
/// \ingroup Attributes
template <typename T, typename U, typename ... V>
struct attribute_of : public attribute_of<T, U>, public attribute_of<T, V...>
{
    using attribute_of<T, U>::operator();
    using attribute_of<T, V...>::operator();
};

/// functor class to simply return argument
struct identity
{
    /// functor returns argument
    template <typename T>
    constexpr auto operator()(T&& t) const noexcept -> decltype(std::forward<T>(t))
    { return std::forward<T>(t); }
};

/// functor class to check an attribute of an argument
template <typename A>
class check_attribute : with_return_type<const bool>
{
public:
    /// constructor
    /// \param a Attribute functor to use
    /// \param val value to check return value of Attribute functor against
    check_attribute(const A& a, const typename A::return_type& val)
        : Attr_(a), Value_(val) {}
    
    /// constructor
    /// \param val value to check return value of attribute functor against
    /// Attribute functor is default constructed
    check_attribute(const typename A::return_type& val)
        : Value_(val) {}

    /// functor
    template <typename U>
    const bool operator()(const U& u) const
    { return Attr_(u) == Value_; }

private:
    /// attribute object
    A Attr_;

    /// value to check against
    typename A::return_type Value_;
};

/// Functor template class for checking for a member of a set by pointer
/// \ingroup Attributes
template <typename T, typename U, typename ... V>
class has_pointed_to_object : public attribute_of<const bool, U, V...>
{
public:
    /// Constructor
    explicit has_pointed_to_object(const std::vector<const T*>& tv)
        : attribute_of<const bool, U, V...>(), Objects_(tv) {}

    /// Constructor
    template <typename ... Others>
    explicit has_pointed_to_object(const T* t, Others ... others)
        : has_pointed_to_object(std::vector<const T*>({t, others...})) {}

    /// Constructor
    explicit has_pointed_to_object(const std::vector<std::shared_ptr<const T> >& tv)
        : attribute_of<const bool, U, V...>(), Objects_(raw_pointers(tv)) {}

    /// Constructor
    template <typename ... Others>
    explicit has_pointed_to_object(const std::shared_ptr<const T>& t, Others ... others)
        : has_pointed_to_object(std::vector<std::shared_ptr<const T> >({t, others...})) {}

    /// Constructor
    template <typename ... Others>
    explicit has_pointed_to_object(const T& t, Others ... others)
        : has_pointed_to_object(std::vector<const T*>({&t, &others...})) {}

    /// \return Objects_
    const std::vector<const T*>& objects() const
    { return Objects_; }

private:
    /// vector of pointers to check for
    std::vector<const T*> Objects_;
};

/// functor class to compare arguments by an attribute
template <typename A, typename C = std::less<typename A::return_type> >
class compare_by : with_return_type<const bool>
{
public:
    compare_by() = default;

    /// constructor
    /// \param a Attribute functor to use
    /// \param comp Comparitor to use
    compare_by(const A& a, const C& comp) : Attr_(a), Comp_(comp) {}

    /// constructor
    /// \param a Attribute functor to use
    compare_by(const A& a) : Attr_(a) {}

    /// constructor
    /// \param comp Comparitor to use
    compare_by(const C& comp) : Comp_(comp) {}

    /// functor
    template <typename T, typename U>
    const bool operator()(const T& t, const U& u) const
    { return Comp_(Attr_(t), Attr_(u)); }

private:
    /// attribute object
    A Attr_;

    /// value to check against
    C Comp_;
};

}

#endif
