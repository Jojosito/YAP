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

#ifndef yap_Constants_h
#define yap_Constants_h

#include "CoordinateSystem.h"
#include "FourVector.h"
#include "ThreeVector.h"

#include <complex>

namespace yap {

/// \name Complex constants
/// @{
/// complex i
constexpr auto Complex_i = std::complex<double>(0, 1);
/// @}

/// \name real constants
/// @{

/// pi
template <typename T = double>
constexpr T pi()
{ return acos((T) - 1); }

/// convert degrees to radians
template <typename T = double>
constexpr T rad_per_deg()
{ return pi<T>() / T(180); }

/// convert degrees to radians
template <typename T>
constexpr T rad(const T& d)
{ return d * rad_per_deg<T>(); }

/// convert radians to degrees
template <typename T = double>
constexpr T deg_per_rad()
{ return T(180) / pi<T>(); }

/// convert radians to degrees
template <typename T>
constexpr T deg(const T& r)
{ return r * deg_per_rad<double>(); }

/// @}

/// \name #ThreeVector constants
/// @{

/// X axis (ThreeVector)
constexpr auto ThreeAxis_X = ThreeVector<double>({1, 0, 0});

/// Y axis (ThreeVector)
constexpr auto ThreeAxis_Y = ThreeVector<double>({0, 1, 0});

/// Z axis (ThreeVector)
constexpr auto ThreeAxis_Z = ThreeVector<double>({0, 0, 1});

/// Standard 3D coordinate system
constexpr auto ThreeAxes = CoordinateSystem<double, 3> {ThreeAxis_X, ThreeAxis_Y, ThreeAxis_Z};

/// 0 as ThreeVector;
constexpr auto ThreeVector_0 = ThreeVector<double>({0, 0, 0});

/// @}

/// \name #FourVector constants
/// @{

/// T axis (FourVector)
constexpr auto FourAxis_T = FourVector<double>({1, 0, 0, 0});

/// X axis (FourVector)
constexpr auto FourAxis_X = FourVector<double>({0, 1, 0, 0});

/// Y axis (FourVector)
constexpr auto FourAxis_Y = FourVector<double>({0, 0, 1, 0});

/// Z axis (FourVector)
constexpr auto FourAxis_Z = FourVector<double>({0, 0, 0, 1});

/// Standard 4D coordinate system
constexpr auto FourAxes = CoordinateSystem<double, 4>({FourAxis_T, FourAxis_X, FourAxis_Y, FourAxis_Z});

/// 0 as FourVector;
constexpr auto FourVector_0 = FourVector<double>({0, 0, 0, 0});

/// @}

}

#endif
