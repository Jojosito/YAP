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

#ifndef yap_BlattWeisskopf_h
#define yap_BlattWeisskopf_h

#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/DecayingParticle.h"
#include "fwd/CachedValue.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/StatusManager.h"

#include "AmplitudeComponent.h"

#include <complex>
#include <memory>
#include <string>

namespace yap {

/// \class BlattWeisskopf
/// \brief Class implementing BlattWeisskopf barrier factors
/// \author Johannes Rauch, Daniel Greenwald

class BlattWeisskopf : public RecalculableAmplitudeComponent
{
public:

    /// Constructor
    /// \param L angular momentum of Blatt-Weisskopf barrier factor
    /// \param dp raw pointer to owning DecayingParticle
    BlattWeisskopf(unsigned L, DecayingParticle* dp);

    /// \return angular momentum
    unsigned L() const
    { return L_; }

    /// \return Blatt-Weisskopf barrier factor for data point and particle combination
    /// \param d DataPoint
    /// \param pc shared_ptr to ParticleCombination
    virtual const std::complex<double> value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const override;

    /// Calculate barrier factors for and store into each data point in a data partition
    /// \param D DataPartition to calculate over
    virtual void calculate(DataPartition& D) const override;

    /// update the calculationStatus for a DataPartition
    virtual void updateCalculationStatus(StatusManager& D) const override;

    /// \return raw pointer to Model through owning DecayingParticle
    const Model* model() const override;

    const DecayingParticle* decayingParticle() const
    { return DecayingParticle_; }

    /// grant friend status to DecayingParticle to call addParticleCombination
    friend class DecayingParticle;

protected:

    /// override to throw on adding non-two-body PC
    void addParticleCombination(const ParticleCombination& pc) override;

private:

    /// raw pointer to owning DecayingParticle
    DecayingParticle* DecayingParticle_;

    /// angular momentum
    unsigned L_;

    /// Blatt-Weisskopf barrier factor
    std::shared_ptr<RealCachedValue> BarrierFactor_;

};

/// squared Blatt-Weisskopf barrier factor;
/// \param l Orbital angular momentum
/// \param z squared breakup momentum / squared radius
/// \warning approximate result for l >= 8 is not properly tested, nor streamlined!
const double squared_barrier_factor(unsigned l, double z);
}

#endif
