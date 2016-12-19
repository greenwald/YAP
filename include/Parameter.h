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

#ifndef yap_Parameter_h
#define yap_Parameter_h

#include "fwd/Parameter.h"

#include "Exceptions.h"
#include "VariableStatus.h"

#include <algorithm>
#include <complex>
#include <memory>
#include <numeric>

namespace yap {

/// \class ParameterBase
/// \brief Class holding basic properties of a parameter, but not a value!
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup Parameters
class ParameterBase
{
public:

    /// Constructor
    ParameterBase() = default;

    /// virtual destructor
    virtual ~ParameterBase() = default;

    /// copy constructor
    ParameterBase(const ParameterBase&) = default;
    
    /// move constructor
    ParameterBase(ParameterBase&&) = default;

    /// copy assignment
    ParameterBase& operator=(const ParameterBase&) = default;

    /// move assignment
    ParameterBase& operator=(ParameterBase&&) = default;

    /// \return VariableStatus
    VariableStatus& variableStatus()
    { return VariableStatus_; }

    /// \return VariableStatus (const)
    const VariableStatus variableStatus() const
    { return VariableStatus_; }

    /// \return number of real elements in parameter
    virtual const size_t size() const = 0;

    /// Set value from vector
    virtual const VariableStatus setValue(const std::vector<double>& V) = 0;

private:

    /// Status of variable
    VariableStatus VariableStatus_{VariableStatus::changed};

};

/// \return VariableStatus::changed if any element in range is changed; VariableStatus::unchanged otherwise
template <typename IterType>
constexpr VariableStatus variable_status(IterType first, IterType last)
{
    return std::any_of(first, last,
                       [](const std::shared_ptr<ParameterBase>& p)
                       {return p->variableStatus() == VariableStatus::changed;})
        ? VariableStatus::changed : VariableStatus::unchanged;
}

/// \class Parameter
/// \brief Template class holding also a value for a parameter
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Parameters
/// \tparam T type to be stored in Parameter
template <typename T>
class Parameter : public ParameterBase
{
public:

    /// Value-assigning constructor
    explicit Parameter(T t) : ParameterBase(), ParameterValue_(t) {}

    /// \return value of parameter
    virtual typename std::conditional<std::is_fundamental<T>::value, const T, const T&>::type
    value() const
    { return ParameterValue_; }

    using ParameterBase::setValue;

    /// set value
    virtual const VariableStatus setValue(typename std::conditional<std::is_fundamental<T>::value, const T, const T&>::type val)
    {
        if (variableStatus() == VariableStatus::fixed)
            throw exceptions::ParameterIsFixed("", "Parameter::setValue");
        if (ParameterValue_ == val)
            return variableStatus();
        ParameterValue_ = val;
        return variableStatus() = VariableStatus::changed;
    }

    /// set value by operator
    Parameter& operator=(typename std::conditional<std::is_fundamental<T>::value, const T, const T&>::type t)
    { setValue(t); return *this; }


private:

    /// Value stored
    T ParameterValue_;

};

/// \return string of Parameter
template <typename T>
inline std::string to_string(const Parameter<T>& P)
{ using std::to_string; return to_string(P.value()) + " (" + to_string(P.variableStatus()) + ")"; }

/// \class RealParameter
/// \ingroup Parameters
class RealParameter : public Parameter<double>
{
public:

    /// constructor
    explicit RealParameter(double t) : Parameter(t) {}

    /// \return size = 1
    virtual const size_t size() const override
    { return 1; }

    using Parameter::setValue;

    /// Set value from vector
    virtual const VariableStatus setValue(const std::vector<double>& V) override
    { return setValue(V[0]); }

    using Parameter::operator=;
};

/// \class NonnegativeRealParameter
/// \ingroup Parameters
class NonnegativeRealParameter : public RealParameter
{
public:
    /// Constructor
    explicit NonnegativeRealParameter(double t) : RealParameter(t)
    {
        if (value() < 0)
            throw exceptions::Exception("value is negative", "NonnegativeRealParameter::NonnegativeRealParameter");
    }

    using RealParameter::setValue;

    /// Checks if value is non-negative
    virtual const VariableStatus setValue(double t) override
    {
        if (t < 0)
            throw exceptions::Exception("value is negative", "NonnegativeRealParameter::setValue");
        return RealParameter::setValue(t);
    }

    using RealParameter::operator=;
};

/// \class PositiveRealParameter
/// \ingroup Parameters
class PositiveRealParameter : public NonnegativeRealParameter
{
public:
    /// Constructor
    explicit PositiveRealParameter(double t) : NonnegativeRealParameter(t)
    {
        if (t <= 0)
            throw exceptions::Exception("value is not positive", "PositiveRealParameter::PositiveRealParameter");
    }
    
    using NonnegativeRealParameter::setValue;

    /// Checks if value is positive
    virtual const VariableStatus setValue(double t) override
    {
        if (t <= 0)
            throw exceptions::Exception("value is not positive", "NonnegativeRealParameter::setValue");
        return RealParameter::setValue(t);
    }

    using NonnegativeRealParameter::operator=;
};

/// \class ComplexParameter
/// \ingroup Parameters
class ComplexParameter : public Parameter<std::complex<double> >
{
public:

    /// constructor
    explicit ComplexParameter(const std::complex<double>& t) : Parameter(t) {}

    /// \return size = 2
    virtual const size_t size() const override
    { return 2; }

    using Parameter::setValue;

    /// Set value from vector
    virtual const VariableStatus setValue(const std::vector<double>& V) override
    { return setValue(std::complex<double>(V[0], V[1])); }

    using Parameter::operator=;
};


/// \class ComplexComponentParameter
/// \brief Abstract base allowing access to the components of a ComplexParameter as a RealParameter
/// \author Daniel Greenwald
class ComplexComponentParameter : public RealParameter
{
public:
    /// Constructor
    /// \param par shared_ptr to ComplexParameter to access
    explicit ComplexComponentParameter(std::shared_ptr<ComplexParameter> par) : RealParameter(0.), Parent_(par)
    { if (!Parent_) throw exceptions::Exception("Parent unset", "ComplexComponentParameter::ComplexComponentParameter"); }
        
    /// \return value of parameter by accessing parent
    virtual const double value() const override
    { return component(Parent_->value()); }

    using RealParameter::setValue;

    /// set value by accessing parent
    virtual const VariableStatus setValue(double val) override
    { return Parent_->setValue(setComponent(Parent_->value(), val)); }
        
    /// \return shared_ptr to parent
    std::shared_ptr<ComplexParameter> parent()
    { return Parent_; }

    using RealParameter::operator=;

protected:

    /// \return component value from whole value
    virtual double component(const std::complex<double>& c) const = 0;

    /// \return complex value with component changed
    virtual std::complex<double> setComponent(const std::complex<double>& c, double v) const = 0;

private:

    /// ComplexParameter this object points to a component of
    std::shared_ptr<ComplexParameter> Parent_;

};

/// \class RealComponentParameter
/// \brief RealParameter accessing real component of ComplexParameter
/// \author Daniel Greenwald
class RealComponentParameter : public ComplexComponentParameter
{
public:

    /// Constructor
    /// \param par shared_ptr to ComplexParameter to access
    explicit RealComponentParameter(std::shared_ptr<ComplexParameter> par)
        : ComplexComponentParameter(par) {}

    using ComplexComponentParameter::operator=;

protected:

    /// \return component value from whole value
    virtual double component(const std::complex<double>& c) const override
    { return real(c); }

    /// \return complex value with component changed
    virtual std::complex<double> setComponent(const std::complex<double>& c, double v) const override
    { return std::complex<double>(v, imag(c)); }

};

/// \class ImaginaryComponentParameter
/// \brief ImaginaryParameter accessing imaginary component of ComplexParameter
/// \author Daniel Greenwald
class ImaginaryComponentParameter : public ComplexComponentParameter
{
public:

    /// Constructor
    /// \param par shared_ptr to ComplexParameter to access
    explicit ImaginaryComponentParameter(std::shared_ptr<ComplexParameter> par)
        : ComplexComponentParameter(par) {}

    using ComplexComponentParameter::operator=;

protected:

    /// \return component value from whole value
    virtual double component(const std::complex<double>& c) const override
    { return imag(c); }

    /// \return complex value with component changed
    virtual std::complex<double> setComponent(const std::complex<double>& c, double v) const override
    { return std::complex<double>(real(c), v); }

};

/// \return number of real elements in ParameterVector
inline const size_t size(const ParameterVector& V)
{ return std::accumulate(V.begin(), V.end(), size_t(0), [](size_t& s, const ParameterVector::value_type & p) {return s += p->size();}); }

/// set values in ParameterVector from iterators
template <class InputIt>
void set_values(ParameterVector::iterator first_par, ParameterVector::iterator last_par, InputIt first_val, InputIt last_val)
{
    while (first_par != last_par and first_val != last_val and (first_val + (*first_par)->size() - 1) != last_val) {
        (*first_par)->setValue(std::vector<double>(first_val, first_val + (*first_par)->size()));
        first_val += (*first_par++)->size();
    }
    if (first_par != last_par)
        throw exceptions::Exception("insufficient number of values provided", "set_values");
}

/// set values in ParameterVector from values in vector<double>
inline void set_values(ParameterVector& pars, const std::vector<double>& vals)
{ set_values(pars.begin(), pars.end(), vals.begin(), vals.end()); }

}

#endif
