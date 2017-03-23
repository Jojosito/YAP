#include "BowlerMassShape.h"

#include "../../examples/bat_gen/models/d4pi_scales.h"
#include "../../data/set_parities.h"
#include "../../data/deduce_parities.h"

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
#include "MathUtilities.h"
#include "Model.h"
#include "ModelIntegral.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleTable.h"
#include "PDL.h"

#include "QuantumNumbers.h"
#include "RelativisticBreitWigner.h"
#include "SpinAmplitudeCache.h"

#include <algorithm>
#include <assert.h>
#include <thread>


namespace yap {

//-------------------------
BowlerMassShape::BowlerMassShape(double mass, double width) :
    BreitWigner(mass, width),
    KsKCoupling_(std::make_shared<PositiveRealParameter>(0.06)),
    CalculatedForMass_(0.), CalculatedForWidth_(0.)
{}

//-------------------------
BowlerMassShape::BowlerMassShape(const ParticleTableEntry& pde) :
    BreitWigner(pde),
    KsKCoupling_(std::make_shared<PositiveRealParameter>(0.06)),
    CalculatedForMass_(0.), CalculatedForWidth_(0.)
{}

//-------------------------
void BowlerMassShape::lock()
{
    // make the model

    // \todo: not hardcode
    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");
    deduce_meson_parities(T);
    set_parities(T);

    // final state particles
    auto piPlus  = FinalStateParticle::create(T[211]);
    auto piMinus = FinalStateParticle::create(T[-211]);

    auto M = std::make_unique<yap::Model>(std::make_unique<yap::HelicityFormalism>());
    M->setFinalState(piPlus, piMinus, piPlus);

    // use common radial size for all resonances
    double r = 1.2; // [GeV^-1]

    // initial state particle
    auto a_1 = DecayingParticle::create(T["a_1+"], r);

    // rho
    auto rho = DecayingParticle::create(T["rho0"], r, std::make_shared<BreitWigner>(T["rho0"]));
    rho->addStrongDecay(piPlus, piMinus);

    // sigma / f_0(500)
    auto sigma = DecayingParticle::create(T["f_0(500)"], r, std::make_shared<BreitWigner>(T["f_0(500)"]));
    sigma->addStrongDecay(piPlus, piMinus);

    a_1->addStrongDecay(rho,   piPlus);
    a_1->addStrongDecay(sigma, piPlus);

    M->lock();

    amp_a1_rho_S_ = free_amplitude(*a_1, to(rho), l_equals(0));
    amp_a1_rho_D_ = free_amplitude(*a_1, to(rho), l_equals(2));
    amp_a1_sigma_ = free_amplitude(*a_1, to(sigma));

    auto model_a_1 = std::static_pointer_cast<DecayingParticle>(particle(*model(), is_named("a_1+")));
    auto model_rho = std::static_pointer_cast<DecayingParticle>(particle(*model(), is_named("rho0")));
    auto model_sigma = std::static_pointer_cast<DecayingParticle>(particle(*model(), is_named("f_0(500)")));
    amp_model_a1_rho_S_ = free_amplitude(*model_a_1, to(model_rho), l_equals(0));
    amp_model_a1_rho_D_ = free_amplitude(*model_a_1, to(model_rho), l_equals(2));
    amp_model_a1_sigma_ = free_amplitude(*model_a_1, to(model_sigma));

    Model_.swap(M);
}

//-------------------------
bool BowlerMassShape::amplitudeChanged() const
{
    if (amp_model_a1_rho_S_->variableStatus() == VariableStatus::changed or
            amp_model_a1_rho_D_->variableStatus() == VariableStatus::changed or
            amp_model_a1_sigma_->variableStatus() == VariableStatus::changed)
        return true;

    return false;
}

//-------------------------
void BowlerMassShape::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    if (MassDependentWidth_.empty() or
        mass() ->value() != CalculatedForMass_ or
        width()->value() != CalculatedForWidth_  or
        amplitudeChanged()) {

        calculateMassDependentWidth();
    }

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

    if (not MassDependentWidth_.empty() and
        mass() ->value() == CalculatedForMass_ and
        width()->value() == CalculatedForWidth_ and
        amp_a1_rho_S_->value() == amp_model_a1_rho_S_->value() and
        amp_a1_rho_D_->value() == amp_model_a1_rho_D_->value() and
        amp_a1_sigma_->value() == amp_model_a1_sigma_->value())
        return;

    LOG(INFO) << "BowlerMassShape::calculateMassDependentWidth() - recalculate";

    amp_a1_rho_S_->setValue(amp_model_a1_rho_S_->value());
    amp_a1_rho_D_->setValue(amp_model_a1_rho_D_->value());
    amp_a1_sigma_->setValue(amp_model_a1_sigma_->value());

    LOG(INFO) << amp_a1_rho_D_->value();
    LOG(INFO) << amp_a1_sigma_->value();

    //-------------------------
    static const unsigned n_integrationPoints = 1e4;
    static const unsigned n_threads = 4;//std::max(1u, std::thread::hardware_concurrency());
    static const unsigned nBins = 150;
    // get FSP mass ranges
    static auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");
    static const double m_pi = T["pi+"].mass();
    static const double low_m = 3.*m_pi;
    static const double hi_m = 1.01 * (T["D0"].mass() - m_pi);
    static ImportanceSamplerGenerator impSampGen(*Model_, n_threads);

    // Get normalizing width
    const double a1_mass = mass()->value(); // nominal mass
    const double phsp = dalitz_phasespace_volume(a1_mass, Model_->finalStateParticles());
    assert(phsp > 0.);
    const double density = impSampGen(a1_mass, phsp * n_integrationPoints);
    assert(density > 0.);

    double norm_width = phsp * density / pow(a1_mass, 3);
    assert(not isnan(norm_width) and norm_width > 0.);

    DEBUG("norm_width = " << norm_width);

    assert(width()->value() > 0.);

    DEBUG("m2; width");
    for (unsigned i = 0; i <= nBins; ++i) {
        // first bin is 0
        if (i == 0) {
            MassDependentWidth_[low_m*low_m] = 0.;
            continue;
        }
        static const double scaling = 1.5; // gives higher density of samples at lower m2, where the width is changing more rapidly

        const double m2 = low_m*low_m + (hi_m*hi_m - low_m*low_m) * pow(i, scaling)/pow(nBins-1, scaling);
        const double mass = sqrt(m2);

        const double phsp = dalitz_phasespace_volume(mass, Model_->finalStateParticles());
        assert(phsp > 0.);

        const double density = impSampGen(mass, phsp * n_integrationPoints);
        assert(density > 0.);

        double w = phsp * density / pow(mass, 3) * width()->value() / norm_width;
        assert(w > 0.);

        DEBUG(m2 << "; " << w);

        MassDependentWidth_[m2] = w;
    }

    CalculatedForMass_  = mass() ->value();
    CalculatedForWidth_ = width()->value();
}

//-------------------------
double BowlerMassShape::massDependentWidth(double m2) const
{
    if (MassDependentWidth_.empty())
        throw exceptions::Exception("Cache is empty", "BowlerMassShape::massDependentWidth");

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

    double w = rel * prev->second + (1. - rel) * next->second;

    // K*K threshold
    double mKK2 = m2 - pow(8.9166000e-01 + 4.9367700e-01, 2);
    if (mKK2 > 0)
        w += KsKCoupling_->value() * sqrt(mKK2);

    return w;
}

}




