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

#ifndef yap_LorentzTransformation_h
#define yap_LorentzTransformation_h

#include "Constants.h"
#include "Matrix.h"
#include "ThreeVector.h"

namespace yap {

/// \return a 4D Lorentz-transformation matrix for a pure rotation
/// \param R #ThreeMatrix defining the rotation
template <typename T>
FourMatrix<T> lorentzTransformation(const ThreeMatrix<T>& R)
{
    FourMatrix<T> L = unitMatrix<T, 4>();
    for (unsigned i = 0; i < R.size(); ++i)
        std::copy(R[i].begin(), R[i].end(), L[i].begin());
    return L;
}

/// \return a 4D Lorentz-transformation matrix for a pure boost
/// \param V #FourVector defining boost
template <typename T>
FourMatrix<T> lorentzTransformation(FourVector<T> V)
{
    // jojosito
    /*FourVector<T> B = (T(1) / V[0]) * V;
    T gamma = T(1) / sqrt(T(1) - norm(vect(B)));
    B *= sqrt(gamma - T(1));
    B[0] = -sqrt(gamma + T(1));

    FourMatrix<T> L = unitMatrix<T, 4>() + outer(B, B);
    //L[0][0] = gamma;
    return L;*/

    // paolo
    V *= T(1) / V[0];
    T scalarBeta = (T) sqrt( V[1]*V[1] + V[2]*V[2] + V[3]*V[3] );
    T gamma = (T) 1 / sqrt( (T) 1 - scalarBeta*scalarBeta );
    FourVector<T> B(
            {
        - gamma*scalarBeta/(gamma-1),
                V[1]/scalarBeta,
                V[2]/scalarBeta,
                V[3]/scalarBeta
            } );
    B *= sqrt(gamma-T(1));

    FourMatrix<T> L = unitMatrix<T, 4>() + outer(B, B);
    L[0][0] = gamma;
    return L;


    /*T gamma = T(1) / abs(V);
    FourVector<T> B = (gamma / V[0]) * V;
    B[0] = 1;

    FourMatrix<T> L = unitMatrix<T, 4>() + outer(B, B);
    L[0][0] = gamma;
    return L;*/
}

/// \return a 4D Lorentz-transformation matrix for a pure boost
/// \param V #ThreeVector defining boost
template <typename T>
FourMatrix<T> lorentzTransformation(const ThreeVector<T>& V)
{
    /*T gamma = T(1) / sqrt(T(1) - norm(V));
    FourVector<T> B(T(1), gamma * V);

    FourMatrix<T> L = unitMatrix<T, 4>() + outer(B, B);
    L[0][0] = gamma;
    return L;*/

    /*T gamma = T(1) / sqrt(T(1) - norm(V));
    FourVector<T> B(T(1), V);

    FourMatrix<T> L = unitMatrix<T, 4>() + outer(B, B);
    L[0][0] = gamma;
    return L;*/

    return lorentzTransformation<T>(FourVector<T>(T(1), V));
}

/// \return a 4D Lorentz-transformation matrix for a pure boost
/// \param fourVecs the sum of these define the boost
template <typename T>
constexpr FourMatrix<T> lorentzTransformation(const std::vector<FourVector<T> >& fourVecs)
{ return lorentzTransformation(std::accumulate(fourVecs.begin(), fourVecs.end(), FourVector<T>({0, 0, 0, 0}))); }

/// \return a 4D Lorentz-transformation matrix for a rotation followed by a boost
/// \param R #ThreeMatrix defining rotation
/// \param V #FourVector defining boost
template <typename T>
constexpr FourMatrix<T> lorentzTransformation(const ThreeMatrix<T>& R, const FourVector<T>& V)
{ return lorentzTransformation<T>(V) * lorentzTransformation<T>(R); }

/// \return a 4D Lorentz-transformation matrix for a rotation followed by a boost
/// \param R #ThreeMatrix defining rotation
/// \param V #ThreeVector defining boost
template <typename T>
constexpr FourMatrix<T> lorentzTransformation(const ThreeMatrix<T>& R, const ThreeVector<T>& V)
{ return lorentzTransformation<T>(V) * lorentzTransformation<T>(R); }

/// \return a 4D Lorentz-transformation matrix for a boost followed by a rotation
/// \param R #ThreeMatrix defining rotation
/// \param V #FourVector defining boost
template <typename T>
constexpr FourMatrix<T> lorentzTransformation(const FourVector<T>& V, const ThreeMatrix<T> R)
{ return lorentzTransformation<T>(R) * lorentzTransformation<T>(V); }

/// \return a 4D Lorentz-transformation matrix for a boost followed by a rotation
/// \param R #ThreeMatrix defining rotation
/// \param V #ThreeVector defining boost
template <typename T>
constexpr FourMatrix<T> lorentzTransformation(const ThreeVector<T>& V, const ThreeMatrix<T>& R)
{ return lorentzTransformation<T>(R) * lorentzTransformation<T>(V); }

}
#endif
