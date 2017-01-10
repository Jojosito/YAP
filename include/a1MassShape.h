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

#ifndef yap_a1MassShape_h
#define yap_a1MassShape_h

#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"
#include "fwd/ParticleTable.h"
#include "fwd/StatusManager.h"

#include "BreitWigner.h"
#include "DecayTreeVectorIntegral.h"

#include <complex>
#include <memory>

namespace yap {

/// \class BreitWigner
/// \brief Class for a_1(1260) mass shape. Its a Breit-Wigner with mass dependent width
/// \author Daniel Greenwald
/// \ingroup MassShapes
///
/// Amplitude is 1 / (mass^2 - s - i*mass*width)\n\n
class a1MassShape : public BreitWigner
{
public:

    /// Constructor
    /// \param mass Mass of resonance [GeV]
    /// \param width Width of resonance [GeV]
    a1MassShape(double mass, double width) :
        BreitWigner(mass, width), Integral_(ownersDecayTrees()) {}

    /// Constructor
    /// \param pde ParticleTableEntry to get mass and width from
    a1MassShape(const ParticleTableEntry& pde) :
        BreitWigner(pde), Integral_(ownersDecayTrees()) {}

    using BreitWigner::calculate;
    
    /// Calculate dynamic amplitude T for particular particle combination and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const override;

private:

    /// Integral
    DecayTreeVectorIntegral Integral_;

};

}

#endif
