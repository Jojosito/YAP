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

#ifndef yap_DecayingParticle_h
#define yap_DecayingParticle_h

#include "DecayChannel.h"
#include "Particle.h"

#include <memory>
#include <vector>

namespace yap {

class FinalStateParticle;

/// \class DecayingParticle
/// \brief Class for a particle that will decay
/// \authors Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class DecayingParticle : public Particle
{
public:

    /// Constructor
    DecayingParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize);

    /// \return Ampltiude for particle
    /// \param d DataPoint to evaluate on
    virtual Amp amplitude(DataPoint& d) override;

    /// Check consistency of object
    virtual bool consistent() const override;

    /// \return vector of pointers to final state particles of channel (recursively checked)
    /// \param channel Channel to return final state particles for
    const std::vector<const FinalStateParticle*> finalStateParticles(unsigned int channel = 0) const;

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param c DecayingParticle takes ownership of c
    void addChannel(DecayChannel* c);

    /// \name Getters
    /// @{

    /// \return Number of decay channels for this object
    unsigned int nChannels() const
    { return Channels_.size(); }

    /// Return Channel i
    const DecayChannel* channel(unsigned i) const
    { return Channels_.at(i).get(); }

    /// \return Radial size [GeV^-1]
    double radialSize() const
    { return RadialSize_; }

    /// @}

    /// \name Setters
    /// @{

    /// Set radial size [GeV^-1]
    void setRadialSize(double r)
    { RadialSize_ = r; }

    /// @}

    /// Print complete decay chain
    void printDecayChain() const
    { printDecayChainLevel(0); }

private:

    void printDecayChainLevel(int level) const;

    /// vector of decay channel objects
    std::vector< std::unique_ptr<yap::DecayChannel> > Channels_;

    /// Radial size parameter [GeV^-1]
    double RadialSize_;
};

}

#endif
