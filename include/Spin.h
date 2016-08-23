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

#ifndef yap_spin_h
#define yap_spin_h

#include "fwd/Spin.h"

#include "MathUtilities.h"

#include <algorithm>
#include <cstdlib>
#include <numeric>
#include <string>

namespace yap {

namespace spin_literals {
    constexpr Spin operator "" _spin(unsigned long long j);
    constexpr Spin operator "" _halfspin(unsigned long long two_j);
}

/// \class SpinProjection
/// \brief Stores a spinProjection
class SpinProjection
{
public:
    constexpr SpinProjection(const Spin& S);

protected:
    
    /// double constructor
    constexpr explicit SpinProjection(double j) : Value_(j) {}

public:

    /// cast to double
    constexpr explicit operator double() const
    { return Value_; }

    friend constexpr SpinProjection operator-(const SpinProjection& S)
    { return SpinProjection(-static_cast<double>(S)); }

    friend constexpr SpinProjection operator*(int lhs, const SpinProjection& rhs)
    { return SpinProjection(lhs * static_cast<double>(rhs)); }

    SpinProjection& operator+=(const SpinProjection& rhs)
    { Value_ += rhs.Value_; return *this; }

    friend SpinProjection operator+(SpinProjection lhs, const SpinProjection& rhs)
    { return lhs += rhs; }

    /// pre-increment operator
    virtual SpinProjection& operator++()
    { using namespace spin_literals; return *static_cast<SpinProjection*>(&(*this += 1_spin)); }

    /// post-increment operator
    SpinProjection operator++(int)
    { SpinProjection temp(*this); operator++(); return temp; }

    SpinProjection& operator-=(const SpinProjection& rhs)
    { Value_ -= rhs.Value_; return *this; }

    friend SpinProjection operator-(SpinProjection lhs, const SpinProjection& rhs)
    { return lhs -= rhs; }

    /// pre-increment operator
    virtual SpinProjection& operator--()
    { using namespace spin_literals; return *static_cast<SpinProjection*>(&(*this -= 1_spin)); }

    /// post-increment operator
    SpinProjection operator--(int)
    { SpinProjection temp(*this); operator--(); return temp; }

    /// equality
    friend constexpr bool operator==(const SpinProjection& A, const SpinProjection& B)
    { return A.Value_ == B.Value_; }

    /// inequality
    friend constexpr bool operator!=(const SpinProjection& A, const SpinProjection& B)
    { return A.Value_ != B.Value_; }
    
    /// less than
    friend constexpr bool operator<(const SpinProjection& A, const SpinProjection& B)
    { return A.Value_ < B.Value_; }

    /// less than or equal
    friend constexpr bool operator<=(const SpinProjection& A, const SpinProjection& B)
    { return A.Value_ <= B.Value_; }

    /// greater than
    friend constexpr bool operator>(const SpinProjection& A, const SpinProjection& B)
    { return A.Value_ > B.Value_; }

    /// greater than or equal
    friend constexpr bool operator>=(const SpinProjection& A, const SpinProjection& B)
    { return A.Value_ >= B.Value_; }

private:
    
    /// SpinProjection value
    double Value_;

};

constexpr bool integer(const SpinProjection& S)
{ return static_cast<double>(S) == floor(static_cast<double>(S)); }

constexpr bool half_integer(const SpinProjection& S)
{ return !integer(S); }

inline std::string to_string(const SpinProjection& S)
{
    return integer(S) ?
        std::to_string(static_cast<int>(static_cast<double>(S)))
        :
        std::to_string(static_cast<int>(2 * static_cast<double>(S))) + "/2";
}

class Spin : public SpinProjection
{
private:

    /// double constructor
    constexpr explicit Spin(unsigned int int j) : SpinProjection(0.5 * j) {}

public:
    using Spin::operator++;

    /// pre-increment operator
    virtual Spin& operator++() override
    { using namespace spin_literals; return *static_cast<Spin*>(&(*this += 1_halfspin)); }

    friend constexpr SpinProjection spin_literals::operator "" _spin(unsigned long long j)
    { return SpinProjection(2 * j); }

    friend constexpr SpinProjection spin_literals::operator "" _halfspin(unsigned long long two_j)
    { return SpinProjection(two_j); }

};


/// convert 2*J to string (e.g. 1/2, 1, 3/2, etc.)
inline std::string spin_to_string(int twoJ)
{ return is_even(twoJ) ? std::to_string(twoJ / 2) : std::to_string(twoJ) + "/2"; }

/// convert SpinVector to string
inline std::string to_string(const SpinVector& two_j)
{
    return std::accumulate(two_j.begin(), two_j.end(), std::string(""),
                           [](const std::string & s, const SpinVector::value_type & j)
    {return s + " + " + spin_to_string(j);}).erase(0, 3);
}

/// convert SpinVector to string
inline std::string to_string(const SpinProjectionVector& two_m)
{
    return std::accumulate(two_m.begin(), two_m.end(), std::string(""),
                           [](const std::string & s, const SpinProjectionVector::value_type & m)
    {return s + " + " + spin_to_string(m);}).erase(0, 3);
}

/// \return whether three spins fulfill the triangle relationship
/// \param two_a 2 * spin a
/// \param two_b 2 * spin b
/// \param two_c 2 * spin c
/// \return \f$ \Delta(abc) \f$
constexpr bool triangle(unsigned two_a, unsigned two_b, unsigned two_c)
{ return is_even(two_a + two_b + two_c) and std::abs((int)two_a - (int)two_b) <= (int)two_c and two_c <= (two_a + two_b); }

/// \return Whether angular momentum is conserved in J -> j1 + j2 with orbital angular momentum l
/// \param two_J 2 * spin of initial state
/// \param two_j1 2 * spin of first daughter
/// \param two_j2 2 * spin of second daughter
/// \param l orbital angular momentum
constexpr bool conserves(unsigned two_J, unsigned two_j1, unsigned two_j2, int l)
{
    // check that the spins are consistent, and that the triangle requirements for (Jlj) and (j1j2j) can be met simultaneously
    return is_even(two_J + two_j1 + two_j2)
           and (std::min<int>(two_j1 + two_j2, two_J + 2 * l) >= std::max(std::abs((int)two_j1 - (int)two_j2), std::abs((int)two_J - 2 * l)));
}

/// \return vector of all spin projections, from -two_j to two_j
/// \param two_j spin to make projections of
inline const SpinProjectionVector projections(unsigned two_j)
{
    SpinProjectionVector spv;
    spv.reserve(two_j + 1);
    for (int two_m = -two_j; two_m <= (int)two_j; two_m += 2)
        spv.push_back(two_m);
    return spv;
}

/// \return vector of all possible spin projection states of spins in two_J
/// \param two_J SpinVector of spins to make projections of
inline const std::vector<SpinProjectionVector> projections(const SpinVector& two_J)
{
    // initialize vector of spin projections to -two_j
    std::vector<int> two_M;
    two_M.reserve(two_J.size());
    std::transform(two_J.begin(), two_J.end(), std::back_inserter(two_M),
    [](const SpinVector::value_type & two_j) {return -two_j;});

    std::vector<SpinProjectionVector> SPV;
    // fill SPV with "odometer"-style looping
    while (two_M.back() <= (int)two_J.back()) {
        SPV.push_back(two_M);
        two_M[0] += 2;
        for (size_t i = 0; (i < two_M.size() - 1) and (two_M[i] > (int)two_J[i]); ++i) {
            two_M[i] = -two_J[i];
            two_M[i + 1] += 2;
        }
    }
    return SPV;
}

}

#endif
