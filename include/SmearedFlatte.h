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

#ifndef yap_SmearedFlatte_h
#define yap_SmearedFlatte_h

#include "Flatte.h"

namespace yap {

/// \class SmearedFlatte
/// \brief Class for Flatte resonance shape smeared with Gaussian
/// \author Johannes Rauch
/// \ingroup MassShapes
///
class SmearedFlatte : public Flatte
{
public:

    /// Constructor
    /// \param m mass [GeV]
    /// \param sigma width of the gaussian [GeV]
    /// \param m_min min mass that will occur
    /// \param m_max max mass that will occur
    SmearedFlatte(double m, double sigma, double m_min, double m_max):
        Flatte(m), sigma_(sigma), m_min_(m_min), m_max_(m_max)
    { ; }

    /// Constructor
    /// \param pde ParticleTableEntry to take mass from
    SmearedFlatte(const ParticleTableEntry& pde, double sigma, double m_min, double m_max):
        Flatte(pde), sigma_(sigma), m_min_(m_min), m_max_(m_max)
    { ; }

    using Flatte::calculate;
    
    /// Calculate dynamic amplitude T for and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const override;


private:
    
    void gaussianiir1d(double* data, long length, double sigma, unsigned numsteps) const;
    
    double interpolate(double m, std::vector<double>& shape, double m_min, double m_step) const;

    /// width of the gaussian [GeV]
    double sigma_;
    double m_min_;
    double m_max_;

};

}

#endif
