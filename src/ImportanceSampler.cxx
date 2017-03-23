#include "ImportanceSampler.h"

#include "AmplitudeBasis.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "DataSet.h"
#include "DecayTree.h"
#include "DecayTreeVectorIntegral.h"
#include "Exceptions.h"
#include "FourVector.h"
#include "IntegralElement.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "PHSP.h"
#include "MassRange.h"
#include "Model.h"
#include "ModelIntegral.h"
#include "VariableStatus.h"

#include "logging.h"

#include <assert.h>
#include <future>
#include <vector>

namespace yap {

//-------------------------
void ImportanceSampler::calculate(std::vector<std::complex<double> >& A, const DecayTreeVectorIntegral& I, const DataPoint& d)
{
    if (I.decayTrees().size() != A.size())
        throw exceptions::Exception("size mismatch", "calculate");
    // calculate the amplitudes of all decay trees
    std::transform(I.decayTrees().begin(), I.decayTrees().end(), A.begin(),
                   [&d](const std::shared_ptr<DecayTree>& dt) {return dt->dataDependentAmplitude(d);});
}

//-------------------------
void ImportanceSampler::update(const std::vector<std::complex<double> >& A, DecayTreeVectorIntegral& I, unsigned n)
{
    // increase number of points
    ++n;

    for (size_t i = 0; i < A.size(); ++i) {

        // calculate difference from mean
        double delta_diag = norm(A[i]) - I.diagonals()[i].value();
        // update mean
        diagonals(I)[i].value() += delta_diag / n;
        
        for (size_t j = i + 1; j < A.size(); ++j) {
            // calculate difference from mean
            auto delta_offdiag = conj(A[i]) * A[j] - I.offDiagonals()[i][j - i - 1].value();
            // update mean
            offDiagonals(I)[i][j - i - 1].value() += delta_offdiag / static_cast<double>(n);
        }
    }
}


//-------------------------
unsigned ImportanceSampler::calculate_partition(std::vector<DecayTreeVectorIntegral*>& J, DataPartition& D)
{
    if (J.empty())
        throw exceptions::Exception("vector is empty", "ImportanceSampler::calculate_partition");

    if (!J[0]->model())
        throw exceptions::Exception("Model is nullptr", "ImportanceSampler::calculate_partition");

    // calculate on data partition
    J[0]->model()->calculate(D);

    unsigned n = 0;
    for (auto& j : J) {
        n = 0;
        std::vector<std::complex<double> > A(j->decayTrees().size());
        // loop over data points
        for (const auto& d : D) {
            calculate(A, *j, d);
            update(A, *j, n++);
        }
    }
    return n;
}

//-------------------------
std::vector<DecayTreeVectorIntegral*> ImportanceSampler::select_changed(ModelIntegral& I)
{
    std::vector<DecayTreeVectorIntegral*> C;
    C.reserve(integrals(I).size());

    for (auto& mci : integrals(I))
        if (std::any_of(mci.Integral.decayTrees().begin(), mci.Integral.decayTrees().end(), &has_changed))
            C.push_back(&mci.Integral);

    return C;
}

//-------------------------
void ImportanceSampler::calculate(ModelIntegral& I, Generator g, unsigned N, unsigned n, unsigned t)
{
    // get DecayTreeVectorIntegral's for DecayTree's that need to be calculated
    auto J = select_changed(I);

    // if nothing requires recalculation, return
    if (J.empty())
        return;

    // reset those to be recalculated
    for (auto& j : J)
        reset(*j);

    if (t <= 1) {

        calculate_subset(J, g, N, n);

    } else {

        // create copies for running with in each partition
        std::vector<std::unique_ptr<DecayTreeVectorIntegral> > m_deleter; // RAII
        m_deleter.reserve(t*J.size());
        std::vector<std::vector<DecayTreeVectorIntegral*> > m(t);
        for (auto& m_j : m) {
            m_j.reserve(J.size());
            for (const auto& j : J) {
                m_deleter.push_back(std::make_unique<DecayTreeVectorIntegral>(*j));
                m_j.push_back(m_deleter.back().get());
            }
        }

        // run over each subset storing number of events used in each calculation
        std::vector<std::future<unsigned> > n_sub;
        n_sub.reserve(m.size());
        // create thread for each partial calculation
        unsigned NN = N;
        for (size_t i = 0; i < m.size(); ++i) {
            int nn = NN / (m.size() - i);
            // copy generator
            Generator gCopy(g);
            n_sub.push_back(std::async(std::launch::async, &ImportanceSampler::calculate_subset,
                                       std::ref(m[i]), std::ref(gCopy), nn, n / t));
            NN -= nn;
        }

        // calculate data fractions:
        // also waits for threads to finish calculating
        std::vector<double> f;
        f.reserve(n_sub.size());
        std::transform(n_sub.begin(), n_sub.end(), std::back_inserter(f), std::mem_fn(&std::future<unsigned>::get));
        double N = std::accumulate(f.begin(), f.end(), 0.);
        std::transform(f.begin(), f.end(), f.begin(), std::bind(std::divides<double>(), std::placeholders::_1, N));
        
        for (size_t i = 0; i < m.size(); ++i)
            for (size_t j = 0; j < J.size(); ++j)
                *J[j] += (*m[i][j] *= f[i]);
    }
}

//-------------------------
unsigned ImportanceSampler::calculate_subset(std::vector<DecayTreeVectorIntegral*>& J, Generator& g, unsigned N, unsigned n)
{
    if (!J[0]->model())
        throw exceptions::Exception("Model is nullptr", "ImportanceSampler::partially_calculate");
                
    // create vectors for amplitude calculation
    std::vector<std::vector<std::complex<double> > > A;
    std::transform(J.begin(), J.end(), std::back_inserter(A),
                   [](const DecayTreeVectorIntegral* j)
                   {return std::vector<std::complex<double> >(j->decayTrees().size());});

    // create DataSet and DataPoint
    auto data = const_cast<Model*>(J[0]->model())->createDataSet();

    // calculate
    for (unsigned k = 0; k < N;) {
        data.clear();
        data.reserve(std::min(n, N - k));
        data.setAll(VariableStatus::changed);
        data.setAll(CalculationStatus::uncalculated);
        std::generate_n(std::back_inserter(data), std::min(n, N - k), g);
        J[0]->model()->calculate(data);
        for (const auto& d : data) {
            for (size_t i = 0; i < J.size(); ++i) {
                calculate(A[i], *J[i], d);
                update(A[i], *J[i], k);
            }
            ++k;
        }
    }
    return N;
}

//-------------------------
void ImportanceSampler::calculate(ModelIntegral& I, DataPartition& D)
{
    // get DecayTreeVectorIntegral's for DecayTree's that need to be calculated
    auto J = select_changed(I);

    // if nothing requires recalculation, return
    if (J.empty())
        return;

    // reset those to be recalculated
    for (auto& j : J)
        reset(*j);

    // calculate it
    calculate_partition(J, D);

}

//-------------------------
void ImportanceSampler::calculate(ModelIntegral& I, DataPartitionVector& DPV)
{
    // if only one data partition, don't thread:
    if (DPV.size() == 1) {
        calculate(I, *DPV[0]);
        return;
    }

    // get DecayTreeVectorIntegral's for DecayTree's that need to be calculated
    auto J = select_changed(I);

    // if nothing requires recalculation, return
    if (J.empty())
        return;

    // reset those to be recalculated
    for (auto& j : J)
        reset(*j);

    // create copies for running with in each partition
    std::vector<std::vector<DecayTreeVectorIntegral*> > m(DPV.size());
    for (auto& m_j : m) {
        m_j.reserve(J.size());
        for (const auto& j : J)
            m_j.push_back(new DecayTreeVectorIntegral(*j));
    }

    // run over each partition storing number of events used in each calculation
    std::vector<std::future<unsigned> > n;
    n.reserve(m.size());
    // create thread for each partial calculation
    for (size_t i = 0; i < m.size(); ++i)
        n.push_back(std::async(std::launch::async, &ImportanceSampler::calculate_partition,
                               std::ref(m[i]), std::ref(*DPV[i])));

    // calculate data fractions:
    // also waits for threads to finish calculating
    std::vector<double> f;
    f.reserve(n.size());
    std::transform(n.begin(), n.end(), std::back_inserter(f), std::mem_fn(&std::future<unsigned>::get));
    double N = std::accumulate(f.begin(), f.end(), 0.);
    std::transform(f.begin(), f.end(), f.begin(), std::bind(std::divides<double>(), std::placeholders::_1, N));

    for (size_t i = 0; i < m.size(); ++i)
        for (size_t j = 0; j < J.size(); ++j)
            *J[j] += (*m[i][j] *= f[i]);

}

//-------------------------
ImportanceSamplerGenerator::ImportanceSamplerGenerator(const Model& m, unsigned n_threads, unsigned seed) :
        N_threads_(n_threads), M_(&m)
{
    for (unsigned i = 0; i < N_threads_; ++i) {
        Rnd_.push_back(std::mt19937(seed));
        Rnd_[i].discard(i * 1e7); // need to scramble, otherwise we get weird artefacts
    }
}

//-------------------------
double ImportanceSamplerGenerator::operator()(double isp_mass, unsigned n_integrationPoints, unsigned n_batchSize)
{
    // mass^2 range
    auto m2r = yap::squared(mass_range(isp_mass, M_->massAxes(), M_->finalStateParticles()));

    // generators
    using Generator = std::function<std::vector<yap::FourVector<double> >()>;
    std::vector<Generator> generators(N_threads_);

    for (unsigned i = 0; i < N_threads_; ++i) {
        generators[i] = std::bind(yap::phsp<std::mt19937>, std::cref(*M_), isp_mass, M_->massAxes(), m2r,
                                  Rnd_[i], std::numeric_limits<unsigned>::max());
    }

    // Model integral
    yap::ModelIntegral Integral_(*M_);

    // get DecayTreeVectorIntegral's for DecayTree's that need to be calculated
    auto J = select_changed(Integral_);

    assert(J[0]->model() == M_);

    // reset those to be recalculated
    for (auto& j : J)
        Integrator::reset(*j);

    if (N_threads_ <= 1)
        calculate_subset(J, generators[0], n_integrationPoints, n_batchSize);
    else {
        // create copies for running with in each partition
        std::vector<std::unique_ptr<DecayTreeVectorIntegral> > m_deleter; // RAII
        m_deleter.reserve(N_threads_ * J.size());
        std::vector<std::vector<DecayTreeVectorIntegral*> > m(N_threads_);
        for (auto& m_j : m) {
            m_j.reserve(J.size());
            for (const auto& j : J) {
                m_deleter.push_back(std::make_unique<DecayTreeVectorIntegral>(*j));
                m_j.push_back(m_deleter.back().get());
            }
        }

        // run over each subset storing number of events used in each calculation
        std::vector<std::future<unsigned> > n_sub;
        n_sub.reserve(N_threads_);
        // create thread for each partial calculation
        for (size_t i = 0; i < N_threads_; ++i) {
            n_sub.push_back(std::async(std::launch::async, &calculate_subset,
                                       std::ref(m[i]), std::ref(generators[i]),
                                       n_integrationPoints/N_threads_, n_batchSize/N_threads_));
        }

        // calculate data fractions:
        // also waits for threads to finish calculating
        std::vector<double> f;
        f.reserve(n_sub.size());
        std::transform(n_sub.begin(), n_sub.end(), std::back_inserter(f), std::mem_fn(&std::future<unsigned>::get));
        double N = std::accumulate(f.begin(), f.end(), 0.);
        std::transform(f.begin(), f.end(), f.begin(), std::bind(std::divides<double>(), std::placeholders::_1, N));

        for (size_t i = 0; i < m.size(); ++i)
            for (size_t j = 0; j < J.size(); ++j)
                *J[j] += (*m[i][j] *= f[i]);
    }

    return integral(Integral_).value();
}

}
