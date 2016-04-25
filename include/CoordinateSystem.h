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

#ifndef yap_CoordinateSystem_h
#define yap_CoordinateSystem_h

#include "Vector.h"
#include "Matrix.h"
#include "ThreeVector.h"

#include <algorithm>
#include <string>
#include <type_traits>

#define VECTOREPSILON 1e-10

namespace yap {

/// \typedef CoordinateSystem
/// \ingroup VectorAlgebra
template <typename T, size_t N>
using CoordinateSystem = std::array<Vector<T, N>, N>;

/// \return string
template <typename T, size_t N>
std::string to_string(const CoordinateSystem<T, N>& C)
{
    std::string s = "(";
    std::for_each(C.begin(), C.end(), [&](const Vector<T, N>& v) {s += to_string(v) + ", ";});
    s.erase(s.size() - 2, 2);
    s += ")";
    return s;
}

/// \return CoordinateSystem with vectors of unit norm
/// \param C CoordinateSystem to base on
template <typename T, size_t N>
const CoordinateSystem<T, N> unit(const CoordinateSystem<T, N>& C)
{
    CoordinateSystem<T, N> uC;
    std::transform(C.begin(), C.end(), uC.begin(), [](const Vector<T, N>& c) {return yap::unit<T, N>(c);});
    return uC;
}

/// multiply a SquareMatrix times each coorindate vector.
/// Useful for rotating coordinate frames
template <typename T, size_t N>
const CoordinateSystem<T, N> operator*(const SquareMatrix<T, N>& M, const CoordinateSystem<T, N>& C)
{
    CoordinateSystem<T, N> MC;
    std::transform(C.begin(), C.end(), MC.begin(), [&](const Vector<T, N>& c) {return M * c; });
    return MC;
}

/// \name Specifically for 3D systems
/// @{

/// \return Whether 3D CoordinateSystem is right handed
/// \param C 3D CoordinateSystem to check
template <typename T>
constexpr bool isRightHanded(const CoordinateSystem<T, 3>& C)
{
    // +[2] is proportional to [0] cross [1]
    return (norm(unit(C[2]) - unit(cross(C[0], C[1]))) < VECTOREPSILON);
}

/// \return Whether 3D CoordinateSystem is left handed
/// \param C 3D CoordinateSystem to check
template <typename T>
constexpr bool isLeftHanded(const CoordinateSystem<T, 3>& C)
{
    // -[2] is proportional to [0] cross [1]
    return (norm(unit(C[2]) + unit(cross(C[0], C[1]))) < VECTOREPSILON);
}

/// \return azimuthal angle of V in coordinate system C, in [-pi, pi]
/// \param V ThreeVector to calculate azimuthal angle of
/// \param C Coordinate frame to measure in
template <typename T>
constexpr T phi(const ThreeVector<T>& V, const CoordinateSystem<T, 3>& C)
{ return angle(V - (V * C[2]) * C[2], C[0]) * ( ((V * C[1]) >= 0) ? T(1) : T(-1) ); }

/// \return polar angle of V in coordinate system C, in [0, pi]
/// \param V ThreeVector to calculate azimuthal angle of
/// \param C Coordinate frame to measure in
template <typename T>
constexpr T theta(const ThreeVector<T>& V, const CoordinateSystem<T, 3>& C)
{ return angle(V, C[2]); }

/// This is the fastest to use if calculating both angles
/// \return azimuthal (0; phi) and polar (1; theta) angles of V in coordinate system C
/// \param V ThreeVector to calculate azimuthal angle of
/// \param C Coordinate frame to measure in
template <typename T>
constexpr std::array<T, 2> angles(const ThreeVector<T>& V, const CoordinateSystem<T, 3>& C)
{ return std::array<T, 2> { phi(V, C), theta(V, C) }; }

/// @}

}
#endif
