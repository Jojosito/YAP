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

#ifndef yap_GounarisSakurai_h
#define yap_GounarisSakurai_h

#include "ConstantWidthBreitWigner.h"

namespace yap {

/// \class GounarisSakurai
/// \brief Class for Gounaris-Sakurai resonance shape
/// \author Johannes Rauch
/// \ingroup MassShapes
class GounarisSakurai : public ConstantWidthBreitWigner
{
public:

    /// Constructor
    /// \param m Mass of resonance [GeV]
    /// \param w Width of resonance [GeV]
    GounarisSakurai(double m, double w) : ConstantWidthBreitWigner(m, w) {}

    /// Constructor
    /// \param pde ParticleTableEntry to take mass and width from
    GounarisSakurai(const ParticleTableEntry& pde) : ConstantWidthBreitWigner(pde) {}

    using ConstantWidthBreitWigner::calculate;

    /// Calculate dynamic amplitude T for and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const override;

private:
    double h(double x, double s0) const;
    double h_prime(double x, double s0) const;
};

}

#endif
