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

/// \files

#ifndef yap_FourVector_h
#define yap_FourVector_h

#include "CoordinateSystem.h"
#include "Matrix.h"
#include "ThreeVector.h"
#include "Vector.h"

#include <algorithm>
#include <type_traits>

namespace yap {

/// \class FourVector
/// \brief Four-Vector handling
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup VectorAlgebra
template <typename T>
class FourVector : public Vector<T, 4>
{
public:
    /// intializer_list constructor
    constexpr FourVector(const std::array<T, 4>& l) noexcept : Vector<T, 4>(l) {}

    /// energy + ThreeVector constructor
    /// \param E 0th component
    /// \param P #ThreeVector component
    constexpr FourVector(const T& E, const ThreeVector<T>& P) noexcept : Vector<T, 4>( {E, P[0], P[1], P[2]}) {}

    /// Default constructor
    FourVector() = default;

    /// Use std::array's assignment operators
    using Vector<T, 4>::operator=;

    /// \return inner (dot) product for 4-vectors
    constexpr T operator*(const Vector<T, 4>& B) const override
    { return this->front() * B.front() - std::inner_product(this->begin() + 1, this->end(), B.begin() + 1, (T)0); }

    /// unary minus for 4-vector,
    /// does not change sign of zero'th component
    constexpr Vector<T, 4> operator-() const override
    { return FourVector<T>((*this)[0], -vect(*this)); }
};

/// \return Spatial #ThreeVector inside #FourVector
template <typename T>
constexpr ThreeVector<T> vect(const FourVector<T>& V) noexcept
{ return ThreeVector<T>({V[1], V[2], V[3]}); }

/// \return boost vector of this #FourVector
template <typename T>
constexpr ThreeVector<T> boost(const FourVector<T>& V)
{ return (V[0] != 0) ? (T(1) / V[0]) * vect(V) : ThreeVector<T>({0, 0, 0}); }

/// multiply a 4x4 matrix times a FourVector
template <typename T>
constexpr FourVector<T> operator*(const FourMatrix<T>& R, const FourVector<T>& V)
{ return FourVector<T>(R * static_cast<Vector<T, 4> >(V)); }

/// apply a three-rotation to a FourVector (rotating only the spatial components)
template <typename T>
constexpr FourVector<T> operator*(const ThreeMatrix<T>& R, const FourVector<T>& V)
{ return FourVector<T>(V[0], R * vect(V)); }

/// Calculate helicity frame of V transformed from C,
/// with z = unit(V), y = C.z X z, x = y X z
/// \param V Fourector defining new Z direction
/// \param C CoordinateSystem aiding in defining new Y direction
template <typename T>
typename std::enable_if<std::is_arithmetic<T>::value, CoordinateSystem<T, 3> >::type
helicityFrame(const FourVector<T>& V, const CoordinateSystem<T, 3>& C)
{ return helicityFrame(vect<T>(V), C); }

}
#endif
