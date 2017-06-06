// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D4PI__H
#define __BAT__D4PI__H

#include "../../../data/set_parities.h"
#include "../../../data/deduce_parities.h"

#include "d4pi_scales.h"

#include <AmplitudeBasis.h>
#include <Attributes.h>
#include <BowlerMassShape.h>
#include <BreitWigner.h>
#include <CompensatedSum.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FourMomenta.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <make_unique.h>
#include <MassRange.h>
#include <MathUtilities.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <PHSP.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <SpinAmplitudeCache.h>

#include <assert.h>
#include <complex>
#include <future>
#include <memory>

using namespace yap;

//-------------------------
inline std::unique_ptr<Model> d4pi_phsp()
{
    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus  = FinalStateParticle::create(T[211]);
    auto piMinus = FinalStateParticle::create(T[-211]);

    auto M = std::make_unique<yap::Model>(std::make_unique<yap::HelicityFormalism>());
    M->setFinalState(piPlus, piMinus, piPlus, piMinus);

    // use common radial size for all resonances
    double r = 1.2; // [GeV^-1]

    // initial state particle
    auto D = DecayingParticle::create(T["D0"], r);

    // pi+ pi- flat
    D->addWeakDecay(piPlus, piMinus, piPlus, piMinus);

    M->lock();

    return M;
}

//-------------------------
inline std::unique_ptr<Model> d4pi()
{
    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");
    try {
        deduce_meson_parities(T);
    }
    catch (yap::exceptions::Exception& e) {
        std::cerr << e.what();
    }
    set_parities(T);

    // final state particles
    auto piPlus  = FinalStateParticle::create(T[211]);
    auto piMinus = FinalStateParticle::create(T[-211]);

    auto M = std::make_unique<yap::Model>(std::make_unique<yap::HelicityFormalism>());
    M->setFinalState(piPlus, piMinus, piPlus, piMinus);
    //M->setFinalState(piMinus, piPlus, piMinus, piPlus);

    // use common radial size for all resonances
    double r = 1.2; // [GeV^-1]

    // initial state particle
    auto D = DecayingParticle::create(T["D0"], r);
    
    //
    // resonances
    //
    // a_1
    auto a_1_shape = a1_bowler ? std::make_shared<BowlerMassShape>(T["a_1+"]) : std::make_shared<BreitWigner>(T["a_1+"]);
    if (a1_bowler) {
        //a_1_shape->width()->setValue(0.560);
        a_1_shape->width()->setValue(0.430);
    }

    auto a_1_plus  = DecayingParticle::create(T["a_1+"], r, a_1_shape);
    auto a_1_minus = DecayingParticle::create(T["a_1-"], r, a_1_shape);

    // rho
    auto rho = DecayingParticle::create(T["rho0"], r, std::make_shared<BreitWigner>(T["rho0"]));
    rho->addStrongDecay(piPlus, piMinus);

    // omega
    auto omega = DecayingParticle::create(T["omega"], r, std::make_shared<BreitWigner>(T["omega"]));
    omega->addStrongDecay({piPlus, piMinus});

    // sigma / f_0(500)
    auto sigma = DecayingParticle::create(T["f_0(500)"], r, std::make_shared<BreitWigner>(T["f_0(500)"]));
    sigma->addStrongDecay(piPlus, piMinus);

    // f_0(980) (as Flatte)
    auto f_0_980_flatte = std::make_shared<Flatte>(T["f_0"]);
    f_0_980_flatte->add(FlatteChannel(0.20, *piPlus, *piMinus));
    f_0_980_flatte->add(FlatteChannel(0.50, T[321], T[-321])); // K+K-
    auto f_0_980 = DecayingParticle::create(T["f_0"], r, f_0_980_flatte);
    f_0_980->addStrongDecay(piPlus, piMinus);
                                            
    // f_2(1270)
    auto f_2 = DecayingParticle::create(T["f_2"], r, std::make_shared<BreitWigner>(T["f_2"]));
    f_2->addStrongDecay(piPlus, piMinus); 

    // f_0(1370)
    ParticleTableEntry pdl_f_0_1370(10221, "f_0(1370)", QuantumNumbers(0, 0, 1, 0, 1, 1),
            1.370, {0.350}); // mass, width
    auto f_0_1370 = DecayingParticle::create(pdl_f_0_1370, r, std::make_shared<BreitWigner>(pdl_f_0_1370));
    f_0_1370->addStrongDecay(piPlus, piMinus);

    auto pipi = DecayingParticle::create("pipi", QuantumNumbers(0, 0), r);
    pipi->addWeakDecay(piPlus, piMinus);

    //
    // decays and amplitudes
    //
    if (f_0_pipi) {
        D->addWeakDecay(f_0_980, piPlus, piMinus);
        *free_amplitude(*D, to(f_0_980, piPlus, piMinus)) = amp_f_0_pipi;
    }

    if (f_2_pipi) {
        D->addWeakDecay(f_2, pipi);
        *free_amplitude(*D, to(f_2, pipi)) = amp_f_2_pipi;
    }

    if (sigma_pipi) {
        D->addWeakDecay(sigma, piPlus, piMinus);
        *free_amplitude(*D, to(sigma, piPlus, piMinus)) = amp_sigma_pipi;
    }

    if (sigma_f_0_1370) {
        D->addWeakDecay(sigma, f_0_1370);
        *free_amplitude(*D, to(sigma, f_0_1370)) = amp_sigma_f_0_1370;
    }

    if (flat_4pi) {
        D->addWeakDecay(piPlus, piMinus, piPlus, piMinus);
        *free_amplitude(*D, to(piPlus, piMinus, piPlus, piMinus)) = amp_pipiFlat;

        //D->addWeakDecay(pipi, pipi);
        //*free_amplitude(*D, to(pipi, pipi)) = amp_pipiFlat;
    }

    if (rho_rho) {
        D->addWeakDecay(rho, rho);
        assert(free_amplitudes(*D, to(rho, rho)).size() == 3);

        for (auto& fa : free_amplitudes(*D, to(rho, rho)))
            *fa = static_cast<std::complex<double> >(amp_rho_rho[fa->spinAmplitude()->L()]);

        if (omega_omega) {
            D->addWeakDecay(omega, omega);
            assert(free_amplitudes(*D, to(omega, omega)).size() == 3);

            for (auto& fa : free_amplitudes(*D, to(omega, omega)))
                *fa = static_cast<std::complex<double> >(amp_omega_omega[fa->spinAmplitude()->L()]);
        }
    }

    // a_1
    if (a_rho_pi_S or a_rho_pi_D) {
        a_1_plus->addStrongDecay(rho, piPlus);
        a_1_minus->addStrongDecay(rho, piMinus);

        if (a_rho_pi_S) {
            *free_amplitude(*a_1_plus, to(rho), l_equals(0)) = 1.; // reference amplitude
            free_amplitude(*a_1_plus, to(rho), l_equals(0))->variableStatus() = VariableStatus::fixed;

            if (a1_shared)
                free_amplitude(*a_1_minus, to(rho), l_equals(0))->shareFreeAmplitude(*free_amplitude(*a_1_plus, to(rho), l_equals(0)));
            else
                *free_amplitude(*a_1_minus, to(rho), l_equals(0)) = 1.;
        }
        else {
            *free_amplitude(*a_1_plus, to(rho), l_equals(0)) = 0.;
            free_amplitude(*a_1_plus, to(rho), l_equals(0))->variableStatus() = VariableStatus::fixed;

            *free_amplitude(*a_1_minus, to(rho), l_equals(0)) = 0.;
            free_amplitude(*a_1_minus, to(rho), l_equals(0))->variableStatus() = VariableStatus::fixed;
        }

        assert(free_amplitudes(*a_1_plus, to(rho), l_equals(1)).empty());
        assert(free_amplitudes(*a_1_minus, to(rho), l_equals(1)).empty());

        if (a_rho_pi_D) {
            *free_amplitude(*a_1_plus, to(rho), l_equals(2)) = amp_a_rho_pi_D;

            if (a1_shared)
                free_amplitude(*a_1_minus, to(rho), l_equals(2))->shareFreeAmplitude(*free_amplitude(*a_1_plus, to(rho), l_equals(2)));
            else
                *free_amplitude(*a_1_minus, to(rho), l_equals(2)) = amp_a_rho_pi_D;
        }
        else {
            *free_amplitude(*a_1_plus, to(rho), l_equals(2)) = 0.;
            free_amplitude(*a_1_plus, to(rho), l_equals(2))->variableStatus() = VariableStatus::fixed;

            *free_amplitude(*a_1_minus, to(rho), l_equals(2)) = 0.;
            free_amplitude(*a_1_minus, to(rho), l_equals(2))->variableStatus() = VariableStatus::fixed;
        }

        if (omega_omega) {
            a_1_plus->addStrongDecay(omega, piPlus);
            a_1_minus->addStrongDecay(omega, piPlus);

            if (a_rho_pi_S) {
                *free_amplitude(*a_1_plus, to(omega), l_equals(0)) = amp_a_omega_pi_S;

                if (a1_shared)
                    free_amplitude(*a_1_minus, to(omega), l_equals(0))->shareFreeAmplitude(*free_amplitude(*a_1_plus, to(omega), l_equals(0)));
                else
                    *free_amplitude(*a_1_minus, to(omega), l_equals(0)) = amp_a_omega_pi_S;
            }
            else {
                *free_amplitude(*a_1_plus, to(omega), l_equals(0)) = 0.;
                free_amplitude(*a_1_plus, to(omega), l_equals(0))->variableStatus() = VariableStatus::fixed;

                *free_amplitude(*a_1_minus, to(omega), l_equals(0)) = 0.;
                free_amplitude(*a_1_minus, to(omega), l_equals(0))->variableStatus() = VariableStatus::fixed;
            }

            assert(free_amplitudes(*a_1_plus, to(omega), l_equals(1)).empty());
            assert(free_amplitudes(*a_1_minus, to(omega), l_equals(1)).empty());

            if (a_rho_pi_D) {
                *free_amplitude(*a_1_plus, to(omega), l_equals(2)) = amp_a_rho_pi_D;

                if (a1_shared)
                    free_amplitude(*a_1_minus, to(omega), l_equals(2))->shareFreeAmplitude(*free_amplitude(*a_1_plus, to(omega), l_equals(2)));
                else
                    *free_amplitude(*a_1_minus, to(omega), l_equals(2)) = amp_a_rho_pi_D;
            }
            else {
                *free_amplitude(*a_1_plus, to(omega), l_equals(2)) = 0.;
                free_amplitude(*a_1_plus, to(omega), l_equals(2))->variableStatus() = VariableStatus::fixed;

                *free_amplitude(*a_1_minus, to(omega), l_equals(2)) = 0.;
                free_amplitude(*a_1_minus, to(omega), l_equals(2))->variableStatus() = VariableStatus::fixed;
            }
        }
    }
    if (a_sigma_pi) {
        a_1_plus->addStrongDecay(sigma, piPlus);
        *free_amplitude(*a_1_plus, to(sigma)) = amp_a_sigma_pi;

        a_1_minus->addStrongDecay(sigma, piMinus);

        if (a1_shared)
            free_amplitude(*a_1_minus, to(sigma))->shareFreeAmplitude(*free_amplitude(*a_1_plus, to(sigma)));
        else
            *free_amplitude(*a_1_minus, to(sigma)) = amp_a_sigma_pi;
    }
    
    if (not free_amplitudes(*a_1_plus).empty()) {
        D->addWeakDecay(a_1_plus, piMinus);
        if (a1_plus)
            *free_amplitude(*D, to(a_1_plus)) = 1.;
        else
            *free_amplitude(*D, to(a_1_plus)) = 0.;
        free_amplitude(*D, to(a_1_plus))->variableStatus() = VariableStatus::fixed;

        D->addWeakDecay(a_1_minus, piPlus);
        if (a1_minus)
            *free_amplitude(*D, to(a_1_minus)) = 1.;
        else
            *free_amplitude(*D, to(a_1_minus)) = 0.;
        free_amplitude(*D, to(a_1_minus))->variableStatus() = VariableStatus::fixed;

    }

    //
    // background
    //
    if (bg_flat_4pi) {
        auto flat_4pi = DecayingParticle::create("flat_4pi", QuantumNumbers(0, 0), r);
        flat_4pi->addWeakDecay(piPlus, piMinus, piPlus, piMinus);
        //flat_4pi->addWeakDecay(pipi, pipi); // produces wrong result since integral would have to be initialized with 4 due to symmetrization
        M->addInitialState(flat_4pi);
    }

    if (bg_rho) {
        auto rho_2pi = DecayingParticle::create("rho_2pi", QuantumNumbers(0, 0), r);
        rho_2pi->addWeakDecay(rho, pipi);
        M->addInitialState(rho_2pi);
        //M->addInitialState(rho);
    }

    if (bg_a1) {
        M->addInitialState(a_1_plus);
        M->addInitialState(a_1_minus);
    }


    //
    // finish and print summary
    //
    M->lock();

    LOG(INFO) << "D Decay trees:";
    LOG(INFO) << to_string(D->decayTrees());

    LOG(INFO) << "\nFree amplitudes: ";
    for (const auto& fa : free_amplitudes(*M, yap::is_not_fixed()))
        LOG(INFO) << yap::to_string(*fa) << "  \t (mag, phase) = (" << abs(fa->value()) << ", " << deg(arg(fa->value())) << "°)"
            << "  \t (real, imag) = (" << real(fa->value()) << ", " << imag(fa->value()) << ")"
            << "  \t" << fa.get()
            << "  \t" << fa->freeAmplitude().get();

    // loop over admixtures
    for (auto& comp : M->components()) {
        if (comp.particle() == &*D) {
            if (bg_only)
                comp.admixture()->setValue(0);
            comp.admixture()->variableStatus() = yap::VariableStatus::fixed;
            assert( comp.admixture()->variableStatus() == yap::VariableStatus::fixed );
            LOG(INFO) << "fixed D0 admixture";
        }
        else if (comp.admixture()->variableStatus() != yap::VariableStatus::fixed)
            *comp.admixture() = 3.0;

        LOG(INFO) << "admixture " << comp.particle()->name()
                << "; m = " << spin_to_string(comp.decayTrees()[0]->initialTwoM())
                << "; variableStatus = " << to_string(comp.admixture()->variableStatus())
                << "; value = " << comp.admixture()->value();
    }

    return M;
}

#endif
