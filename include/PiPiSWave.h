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

#ifndef yap_PiPiSWaveAuMorganPennington_h
#define yap_PiPiSWaveAuMorganPennington_h

#include "MassShape.h"
#include "Matrix.h"

#include <complex>

namespace yap {

/// \class PiPiSWaveAuMorganPennington
/// \brief pi pi S-wave implementation copied from rootPWA
/// \author Johannes Rauch
/// \ingroup MassShapes
class PiPiSWaveAuMorganPennington : public MassShape
{
public:

    /// Constructor
    PiPiSWaveAuMorganPennington();

    using MassShape::calculate;

    /// Calculate dynamic amplitude T for particular particle combination and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const override;

protected:

    SquareMatrix<std::complex<double>, 2> T_;
    std::vector<SquareMatrix<std::complex<double>, 2>> a_;
    std::vector<SquareMatrix<std::complex<double>, 2>> c_;
    Matrix<double, 1, 2> sP_;
    int vesSheet_;

    double piChargedMass_;
    double piNeutralMass_;
    double kaonChargedMass_;
    double kaonNeutralMass_;
    double kaonMeanMass_;

private:

    /// cached dynamic amplitude
    std::shared_ptr<ComplexCachedValue> CachedMassShape_;

};

class PiPiSWaveAuMorganPenningtonKachaev : public PiPiSWaveAuMorganPennington
{
    /// Constructor
    PiPiSWaveAuMorganPenningtonKachaev();

    using PiPiSWaveAuMorganPennington::calculate;
};

}

#endif
