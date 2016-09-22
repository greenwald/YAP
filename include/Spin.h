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

/// convert 2*J to string (e.g. 1/2, 1, 3/2, etc.)
inline std::string spin_to_string(int twoJ)
{ return is_even(twoJ) ? std::to_string(twoJ / 2) : std::to_string(twoJ) + "/2"; }

/// convert SpinVector to string
inline std::string to_string(const SpinVector& two_j)
{
    return std::accumulate(two_j.begin(), two_j.end(), std::string(""),
                           [](const std::string& s, const SpinVector::value_type& j)
                           {return s + " + " + spin_to_string(j);}).erase(0, 3);
}

/// convert SpinVector to string
inline std::string to_string(const SpinProjectionVector& two_m)
{
    return std::accumulate(two_m.begin(), two_m.end(), std::string(""),
                           [](const std::string& s, const SpinProjectionVector::value_type& m)
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
const SpinProjectionVector projections(unsigned two_j);

/// \return vector of all possible spin projection states of spins in two_J
/// \param two_J SpinVector of spins to make projections of
const std::vector<SpinProjectionVector> projections(const SpinVector& two_J);

/// \return vector of possible spins formed from the coupling of provided spins [|A-B| -> (A+B)]
/// \param two_A 2 * first spin
/// \param two_B 2 * second spin
const SpinVector triangle(unsigned two_A, unsigned two_B);

}

#endif
