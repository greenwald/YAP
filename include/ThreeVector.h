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

#ifndef yap_ThreeVector_h
#define yap_ThreeVector_h

#include "Matrix.h"
#include "Vector.h"

namespace yap {

/// \typedef ThreeVector
/// \ingroup VectorAlgebra
template <typename T>
using ThreeVector = Vector<T, 3>;

/// \return cross product
template <typename T>
constexpr ThreeVector<T> cross(const ThreeVector<T>& A, const ThreeVector<T>& B) noexcept
{ return ThreeVector<T>({A[1]* B[2] - A[2]* B[1], A[2]* B[0] - A[0]* B[2], A[0]* B[1] - A[1]* B[0]}); }

/// \return angle between two 3D vectors
template <typename T>
constexpr T angle(const ThreeVector<T>& A, const ThreeVector<T>& B)
{ return acos(A * B / abs(A) / abs(B)); }

/// skew symmetric matrix formed from a 3 vector
template<typename T>
constexpr SquareMatrix<T, 3> skewSymmetric(const ThreeVector<T> V) noexcept
{
    return {
        {    0., -V[2],  V[1] },
        {  V[2],    0., -V[0] },
        { -V[1],  V[0],    0. }
    };
}

}
#endif
