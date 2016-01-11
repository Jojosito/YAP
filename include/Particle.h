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

#ifndef yap_Particle_h
#define yap_Particle_h

#include "AmplitudeComponent.h"
#include "Parameter.h"
#include "QuantumNumbers.h"
#include "ReportsInitialStateParticle.h"
#include "ReportsParticleCombinations.h"

#include <string>

namespace yap {

class ParticleCombination;
class FinalStateParticle;

/// \class Particle
/// \brief Abstract Particle base class.
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup Particle Particle-related classes

// keyword virtual is needed to solve diamond problem in DecayingParticle
class Particle :
    public virtual AmplitudeComponent,
    public virtual ReportsInitialStateParticle,
    public virtual ReportsParticleCombinations
{
public:

    /// Constructor
    /// \param q Quantum numbers of particle
    /// \param m Mass of particle
    /// \param name Name of particle
    Particle(const QuantumNumbers& q, double m, std::string name);

    /// Calculate complex amplitude
    /// must be overrided in derived classes
    /// \param d DataPoint to calculate with
    /// \param pc (shared_ptr to) ParticleCombination to calculate for
    /// \param two_m 2 * the spin projection to calculate for
    /// \param dataPartitionIndex partition index for parallelization
    virtual std::complex<double> amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc,
                                           int two_m, unsigned dataPartitionIndex) const = 0;

    /// Check consitency of object
    virtual bool consistent() const override;

    /// \name Getters
    /// @{

    /// const access QuantumNumbers object
    const QuantumNumbers& quantumNumbers() const
    { return QuantumNumbers_; }

    /// \todo Do we need non-const access to the quantum numbers?
    /// \return quantum numbers
    QuantumNumbers& quantumNumbers()
    { return QuantumNumbers_; }

    /// Get mass [GeV]
    std::shared_ptr<RealParameter> mass() const
    { return Mass_; }

    /// Get name (const)
    const std::string& name() const
    { return Name_; }

    /// Get name
    std::string& name()
    { return Name_; }

    /// @}

    /// give friend status to DecayChannel to call setSymmetrizationIndexParents
    friend class DecayChannel;

protected:

    /// set symmetrization index parents
    virtual void setSymmetrizationIndexParents() = 0;

private:

    /// Quantum numbers of particle
    QuantumNumbers QuantumNumbers_;

    /// Mass [GeV]
    /// \todo share with massShape
    std::shared_ptr<RealParameter> Mass_;

    /// Name of particle
    std::string Name_;

};

/// \typedef ParticleVector
/// \ingroup Particle
using ParticleVector = std::vector<std::shared_ptr<Particle> >;

/// convert to string
inline std::string to_string(const Particle& p)
{ return p.name() + "(" + to_string(p.quantumNumbers()) + "), mass = " + std::to_string(p.mass()->value()); }

/// streamer
inline std::ostream& operator<<(std::ostream& os, const Particle& p)
{ os << to_string(p); return os; }

}

#endif
