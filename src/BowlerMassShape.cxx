#include "BowlerMassShape.h"

#include "AmplitudeBasis.h"
#include "Attributes.h"
#include "BreitWigner.h"
#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DalitzPhspIntegral.h"
#include "DataPartition.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "DecayTree.h"
#include "FinalStateParticle.h"
#include "Flatte.h"
#include "FourMomenta.h"
#include "FourVector.h"
#include "FreeAmplitude.h"
#include "HelicityFormalism.h"
#include "ImportanceSampler.h"
#include "logging.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "MassRange.h"
#include "MathUtilities.h"
#include "Model.h"
#include "ModelIntegral.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleTable.h"
#include "PDL.h"
#include "PHSP.h"
#include "QuantumNumbers.h"
#include "RelativisticBreitWigner.h"
#include "SpinAmplitudeCache.h"

#include <algorithm>
#include <assert.h>
#include <random>

namespace yap {

//-------------------------
void BowlerMassShape::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    /// \todo also recalculate if mass or width have changed
    if (MassDependentWidth_.empty())
        calculateMassDependentWidth();

    // if no calculation necessary, exit
    if (D.status(*T(), si) != CalculationStatus::uncalculated)
        return;

    for (auto& d : D) {
        auto m2 = model()->fourMomenta()->m2(d, pc);
        auto M = mass()->value();
        // T := 1 / (M^2 - m^2 - i * M * Gamma)
        T()->setValue(1. / (M*M - m2 - 1_i * M * massDependentWidth(m2)),
                d, si, D);
    }

    D.status(*T(), si) = CalculationStatus::calculated;
}

//-------------------------
void BowlerMassShape::calculateMassDependentWidth() const
{
    std::lock_guard<std::mutex> guard(CacheMutex_);

    if (not MassDependentWidth_.empty())
        return;

    // \todo: not hardcode
    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus  = FinalStateParticle::create(T[211]);
    auto piMinus = FinalStateParticle::create(T[-211]);

    auto M = std::make_unique<yap::Model>(std::make_unique<yap::HelicityFormalism>());
    M->setFinalState(piPlus, piMinus, piPlus);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto a_1 = DecayingParticle::create(T["a_1+"], radialSize);

    // rho
    auto rho = DecayingParticle::create(T["rho0"], radialSize, std::make_shared<RelativisticBreitWigner>(T["rho0"]));
    rho->addStrongDecay(piPlus, piMinus);

    // sigma / f_0(500)
    auto sigma = DecayingParticle::create(T["f_0(500)"], radialSize, std::make_shared<BreitWigner>(T["f_0(500)"]));
    sigma->addStrongDecay(piPlus, piMinus);


    // a_1 -> sigma pi
    a_1->addStrongDecay(sigma, piPlus);
    *free_amplitude(*a_1, to(sigma)) = std::polar(1.0739853519 * 0.439, rad(193.));

    // a_1 -> rho pi
    a_1->addStrongDecay(rho,   piPlus);
    // S wave
    *free_amplitude(*a_1, to(rho), l_equals(0)) = 1;
    // D wave
    *free_amplitude(*a_1, to(rho), l_equals(2)) = std::polar(1.0739853519 * 0.241, rad(82.));

    M->lock();

    //-------------------------
    const unsigned n_integrationPoints = 1e6;
    const unsigned n_threads = 4;
    const unsigned nBins = 150;

    const double m_pi = T["pi+"].mass();
    const double low_m = 3.*m_pi;
    const double hi_m = T["D0"].mass() - m_pi;



    /// stores integral result
    yap::ModelIntegral Integral_(*M);

    // Get normalizing width
    // get FSP mass ranges
    const double a1_mass = mass()->value(); // nominal mass
    auto m2r = yap::squared(mass_range(a1_mass, M->massAxes(), M->finalStateParticles()));

    /// function for generating new points for integration
    std::mt19937 rnd(164852419);
    using Generator = std::function<std::vector<yap::FourVector<double> >()>;
    Generator g = std::bind(yap::phsp<std::mt19937>, std::cref(*M), a1_mass, M->massAxes(), m2r,
                            rnd, std::numeric_limits<unsigned>::max());
    yap::ImportanceSampler::calculate(Integral_, g, n_integrationPoints, n_integrationPoints, n_threads);

    double norm_width = dalitz_phasespace_volume(a1_mass, M->finalStateParticles())
                        * integral(Integral_).value()
                        / pow(a1_mass, 3. / 2);

    if (isnan(norm_width) or norm_width == 0)
        LOG(ERROR) << "norm_width invalid";
    LOG(INFO) << "norm_width = " << norm_width;

    LOG(INFO) << "m2; width";
    for (unsigned i = 0; i <= nBins; ++i) {

        static const double scaling = 1.5; // gives higher density of samples at lower m2, where the width is changing more rapidly

        const double m2 = low_m*low_m + (hi_m*hi_m - low_m*low_m) * pow(i, scaling)/pow(nBins-1, scaling);
        const double mass = sqrt(m2);

        // get FSP mass ranges
        m2r = yap::squared(mass_range(mass, M->massAxes(), M->finalStateParticles()));
        g = std::bind(yap::phsp<std::mt19937>, std::cref(*M), mass, M->massAxes(), m2r,
                      rnd, std::numeric_limits<unsigned>::max());
        yap::ImportanceSampler::calculate(Integral_, g, n_integrationPoints, n_integrationPoints, n_threads);

        const double phsp = dalitz_phasespace_volume(mass, M->finalStateParticles());
        const double density = integral(Integral_).value();

        double w = phsp * density / pow(mass, 3./2) * width()->value() / norm_width;
        if (std::isnan(w))
            w = 0;
        w = std::max(0., w);

        LOG(INFO) << m2 << "; " << w;

        if (i > 0)
            assert(w > 0);

        MassDependentWidth_[m2] = w;
    }
}

//-------------------------
double BowlerMassShape::massDependentWidth(double m2) const
{
    // linear interpolation
    auto next = MassDependentWidth_.upper_bound(m2);
    assert(next->first >= m2);
    if (next == MassDependentWidth_.begin())
        return next->second;

    auto prev = next;
    --prev;
    double range = next->first - prev->first;
    double dist = m2 - prev->first;
    double rel = dist/range;

    return rel * prev->second + (1. - rel) * next->second;
}

}




