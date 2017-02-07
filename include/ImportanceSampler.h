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

#ifndef yap_ImportanceSampler_h
#define yap_ImportanceSampler_h

#include "Integrator.h"

#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/DecayTreeVectorIntegral.h"
#include "fwd/FourVector.h"
#include "fwd/Model.h"
#include "fwd/ModelIntegral.h"

#include <functional>
#include <random>

namespace yap {

/// \class ImportanceSampler
/// \brief Calculates DecayTreeVectorIntegral using importance sampling
/// \author Daniel Greenwald
/// \ingroup Integration
class ImportanceSampler : public Integrator
{

public:

    /// Update calculation of ModelIntegral
    /// \param I ModelIntegral to calculate
    /// \param DPV vector of DataPartitions to calculate with
    static void calculate(ModelIntegral& I, DataPartitionVector& DPV);

    /// Update calculation of ModelIntegral
    /// \param I ModelIntegral to calculate
    /// \param D DataPartition to calculate with
    static void calculate(ModelIntegral& I, DataPartition& D);

    /// \typedef Generator
    /// function for generating new points for integration
    using Generator = std::function<std::vector<FourVector<double> >()>;

    /// Update calculation of ModelIntegral
    /// \param I ModelIntegral to calculate
    /// \param g Generator to generate new data points (as vectors of FourVector)
    /// \param N number of points to generate
    /// \param n batch size of points to generate
    /// \param t number of threads to use while integrating
    static void calculate(ModelIntegral& I, Generator g, unsigned N, unsigned n, unsigned t = 1);

    /// calculate amplitudes
    static void calculate(std::vector<std::complex<double> >& A, const DecayTreeVectorIntegral& I, const DataPoint& d);
    
    /// update DecayTreeVectorIntegral using amplitudes
    static void update(const std::vector<std::complex<double> >& A, DecayTreeVectorIntegral& I, unsigned n);

    /// \return integral_sub_map for all changed trees
    static std::vector<DecayTreeVectorIntegral*> select_changed(ModelIntegral& I);

    /// perform calculation for one data partition
    static unsigned calculate_partition(std::vector<DecayTreeVectorIntegral*>& J, DataPartition& D);

    /// \param N number of points to generate
    /// \param n batch size of points to generate
    static unsigned calculate_subset(std::vector<DecayTreeVectorIntegral*>& J, Generator& g, unsigned N, unsigned n);

};

class ImportanceSamplerGenerator : private ImportanceSampler
{
public:
    ImportanceSamplerGenerator(const Model& m, unsigned n_threads = 1, unsigned seed = 52350863);

    // generate and return integral
    double operator()(double isp_mass, unsigned n_integrationPoints = 1e4, unsigned n_batchSize = 1e4);

private:
    // number of threads
    const unsigned N_threads_;

    const Model* M_;

    // one generator per thread
    std::vector<std::mt19937> Rnd_;

};

}

#endif
