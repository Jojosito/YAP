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

#ifndef yap_BowlerMassShape_h
#define yap_BowlerMassShape_h

#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"
#include "fwd/ParticleTable.h"
#include "fwd/StatusManager.h"

#include "BreitWigner.h"
#include "Flatte.h"
#include "DecayTreeVectorIntegral.h"

#include <complex>
#include <memory>
#include <mutex>

namespace yap {

/// \class BowlerMassShape
/// \brief Class for mass shape with Bowler parameterization. Its a Flatte with mass dependent width
/// \author Johannes Rauch
/// \ingroup MassShapes
///
class BowlerMassShape : public BreitWigner
{
public:

    /// Constructor
    /// \param mass Mass of resonance [GeV]
    /// \param width Width of resonance [GeV]
    BowlerMassShape(double mass, double width) :
        BreitWigner(mass, width) { }

    /// Constructor
    /// \param pde ParticleTableEntry to get mass and width from
    BowlerMassShape(const ParticleTableEntry& pde) :
        BreitWigner(pde) { }

    using BreitWigner::calculate;
    
    /// Calculate dynamic amplitude T for particular particle combination and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const override;

private:

    void calculateMassDependentWidth() const;

    double massDependentWidth(double m2) const;

    /// mass dependend width (vs m2)
    mutable std::map<double, double> MassDependentWidth_;
    mutable std::mutex CacheMutex_;

};

}

#endif
