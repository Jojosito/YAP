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
    KsKCoupling_(std::make_shared<PositiveRealParameter>(0.06))
{}

//-------------------------
BowlerMassShape::BowlerMassShape(const ParticleTableEntry& pde) :
    BreitWigner(pde),
    KsKCoupling_(std::make_shared<PositiveRealParameter>(0.06))
{}

//-------------------------
void BowlerMassShape::updateCalculationStatus(StatusManager& D) const
{
    BreitWigner::updateCalculationStatus(D);
    if (amp_model_a1_rho_S_->variableStatus() == VariableStatus::changed or
        amp_model_a1_rho_D_->variableStatus() == VariableStatus::changed or
        amp_model_a1_sigma_->variableStatus() == VariableStatus::changed)
        D.set(*T(), CalculationStatus::uncalculated);

    // also set amplitudes for next iteration
    amp_a1_rho_S_->setValue(amp_model_a1_rho_S_->value());
    amp_a1_rho_D_->setValue(amp_model_a1_rho_D_->value());
    amp_a1_sigma_->setValue(amp_model_a1_sigma_->value());
}

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

    fillCache();
}

//-------------------------
void BowlerMassShape::fillCache()
{
    LOG(INFO) << "Fill Bowler cache ...";
    static const unsigned n_integrationPoints = 1e5;
    static const unsigned n_threads = 4;//std::max(1u, std::thread::hardware_concurrency());
    static const unsigned nBins = 150;
    // get FSP mass ranges
    static auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");
    static const double m_pi = T["pi+"].mass();
    static const double low_m = 3.*m_pi;
    static const double hi_m = 1.01 * (T["D0"].mass() - m_pi);
    static ImportanceSamplerGenerator impSampGen(*Model_, n_threads);

    for (unsigned i = 0; i <= nBins; ++i) {
        static const double scaling = 1.5; // gives higher density of samples at lower m2, where the width is changing more rapidly

        const double m2 = low_m*low_m + (hi_m*hi_m - low_m*low_m) * pow(i, scaling)/pow(nBins-1, scaling);
        const double mass = sqrt(m2);

        double phsp = dalitz_phasespace_volume(mass, Model_->finalStateParticles());
        if (i==0)
            phsp = 0;
        else
            assert(phsp > 0.);

        const int nPoints = std::max(unsigned(100), unsigned(phsp * n_integrationPoints));
        Cache_[m2].Phsp_ = phsp;
        Cache_[m2].MI_ = impSampGen.modelIntegral(mass, nPoints, nPoints / 2 + 1);
    }
    LOG(INFO) << "Done.";
}

//-------------------------
void BowlerMassShape::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    // if no calculation necessary, exit
    if (D.status(*T(), si) != CalculationStatus::uncalculated)
        return;

    for (auto& d : D) {
        auto m2 = model()->fourMomenta()->m2(d, pc);
        auto M = mass()->value();
        // T := 1 / (M^2 - m^2 - i * M * Gamma)
        T()->setValue(1. / (M*M - m2 - 1_i * M * massDependentWidth(D, m2)),
                d, si, D);
    }

    D.status(*T(), si) = CalculationStatus::calculated;
}

//-------------------------
double BowlerMassShape::massDependentWidth(DataPartition& D, double m2) const
{
    if (Cache_.empty())
        throw exceptions::Exception("Cache is empty", "BowlerMassShape::massDependentWidth");

    // make sure amplitudes have been set correctly
    assert(amp_a1_rho_D_->value() ==  amp_model_a1_rho_D_->value());

    // norm_width
    auto next = Cache_.upper_bound(pow(mass()->value(), 2));
    auto prev = next;
    if (next != Cache_.begin())
        --prev;
    double rel = (pow(mass()->value(), 2) - prev->first)/(next->first - prev->first);

    double phsp = (1. - rel) * prev->second.Phsp_ + rel * next->second.Phsp_;
    double density = (1. - rel) * integral(prev->second.MI_).value() + rel * integral(next->second.MI_).value();
    const double norm_width = phsp * density / pow(mass()->value(), 3);
    assert(norm_width > 0.);

    // linear interpolation
    next = Cache_.upper_bound(m2);
    prev = next;
    if (next != Cache_.begin())
        --prev;
    rel = (m2 - prev->first)/(next->first - prev->first);

    phsp = (1. - rel) * prev->second.Phsp_ + rel * next->second.Phsp_;
    density = (phsp > 0.) ? (1. - rel) * integral(prev->second.MI_).value() + rel * integral(next->second.MI_).value() : 0.;

    double w = phsp * density / pow(m2, 3./2.) * width()->value() / norm_width;

    // K*K threshold
    // \todo dont hardcode masses
    double mKK2 = m2 - pow(8.9166000e-01 + 4.9367700e-01, 2);
    if (mKK2 > 0)
        w += KsKCoupling_->value() * sqrt(mKK2);

    return w;
}

}




