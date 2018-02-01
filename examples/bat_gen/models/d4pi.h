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
#include <GounarisSakurai.h>
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
#include <PiPiSWave.h>
#include <PoleMass.h>
#include <PDL.h>
#include <PHSP.h>
#include <QuantumNumbers.h>
#include <SmearedFlatte.h>
#include <SpinAmplitudeCache.h>

#include <assert.h>
#include <complex>
#include <future>
#include <memory>

using namespace yap;

//-------------------------
inline void add_3pi_decays(std::shared_ptr<DecayingParticle>& isp,
        std::shared_ptr<FinalStateParticle>& piPlus, std::shared_ptr<FinalStateParticle>& piMinus,
        std::shared_ptr<DecayingParticle>& part_plus, std::shared_ptr<DecayingParticle>& part_minus,
        const std::complex<double>& amp_plus, const std::complex<double>& amp_minus,
        std::vector<std::shared_ptr<DecayingParticle> > daughters,
        std::vector<std::vector<std::complex<double> > > amps_daughters)
{
    // check
    assert(daughters.size() == amps_daughters.size());
    assert(isp->name() == "D0");
    assert(piPlus->name() == "pi+");
    assert(piMinus->name() == "pi-");
    assert(not daughters.empty());

    // add
    for (unsigned i_daughter = 0; i_daughter < daughters.size(); ++i_daughter) {
        part_plus->addStrongDecay<d4pi_max_L>(daughters[i_daughter], piPlus);
        part_minus->addStrongDecay<d4pi_max_L>(daughters[i_daughter], piMinus);

        unsigned L_min(999999);
        unsigned L_max(0);
        for (auto fa : free_amplitudes(*part_plus, to(daughters[i_daughter]))) {
            unsigned L = fa->spinAmplitude()->L();
            if (amps_daughters.at(i_daughter).size() <= L) {
                throw yap::exceptions::Exception(part_plus->name() + " to " + daughters[i_daughter]->name() + ": no amplitude for L = " + std::to_string(L) + " given", "add_3pi_decays");
            }
            if (amps_daughters.at(i_daughter).at(L) == std::complex<double>(0)) {
                LOG(INFO) << part_plus->name() + " to " + daughters[i_daughter]->name() + ": zero amplitude for L = " + std::to_string(L) + " given";
            }
            *fa = amps_daughters.at(i_daughter).at(L);

            L_min = std::min(L_min, L);
            L_max = std::max(L_max, L);

            // share amplitudes
            free_amplitude(*part_minus, to(daughters[i_daughter]), l_equals(L))->shareFreeAmplitude(*fa);
        }

        if (amps_daughters.at(i_daughter).size() > L_max + 1) {
            LOG(INFO) << part_plus->name() + " to " + daughters[i_daughter]->name() + ": max amplitude for L = " + std::to_string(L_max);
        }

        // fix first amplitude
        if (i_daughter == 0) {
            free_amplitude(*part_plus, to(daughters[i_daughter]), l_equals(L_min))->variableStatus() = VariableStatus::fixed;
        }
    }

    isp->addWeakDecay<d4pi_max_L>(part_plus, piMinus);
    isp->addWeakDecay<d4pi_max_L>(part_minus, piPlus);

    *free_amplitude(*isp, to(part_plus)) = amp_plus;

    if (d4pi_pm_shared)
        free_amplitude(*isp, to(part_plus))->shareFreeAmplitude(*free_amplitude(*isp, to(part_minus)));
    else
        *free_amplitude(*isp, to(part_minus)) = amp_plus;
}

//-------------------------
void fix_free_amplitudes(Model& M)
{
    for (const auto& fa : free_amplitudes(M, yap::is_not_fixed()))
        fa->freeAmplitude()->variableStatus() = yap::VariableStatus::fixed;

    for (auto& comp : M.components())
        comp.admixture()->variableStatus() = yap::VariableStatus::fixed;
}

//-------------------------
inline std::unique_ptr<Model> d4pi_phsp()
{
    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/d4pi.pdl");

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
    D->addWeakDecay<d4pi_max_L>(piPlus, piMinus, piPlus, piMinus);

    M->lock();

    return M;
}

//-------------------------
inline void printAmplitudes(const Model& M)
{

    LOG(INFO) << "\nFixed amplitudes: ";
    for (const auto& fa : free_amplitudes(M, yap::is_fixed()))
        LOG(INFO) << yap::to_string(*fa) << "  \t (mag, phase) = (" << abs(fa->value()) << ", " << deg(arg(fa->value())) << "°)"
            << "  \t (real, imag) = (" << real(fa->value()) << ", " << imag(fa->value()) << ")"
            << "  \t" << fa.get()
            << "  \t" << fa->freeAmplitude().get();

    LOG(INFO) << "\nFree amplitudes: ";
    for (const auto& fa : free_amplitudes(M, yap::is_not_fixed()))
        LOG(INFO) << yap::to_string(*fa) << "  \t (mag, phase) = (" << abs(fa->value()) << ", " << deg(arg(fa->value())) << "°)"
            << "  \t (real, imag) = (" << real(fa->value()) << ", " << imag(fa->value()) << ")"
            << "  \t" << fa.get()
            << "  \t" << fa->freeAmplitude().get();

    for (auto& comp : M.components()) {
        LOG(INFO) << "admixture " << comp.particle()->name()
                << "; m = " << spin_to_string(comp.decayTrees()[0]->initialTwoM())
                << "; variableStatus = " << to_string(comp.admixture()->variableStatus())
                << "; value = " << comp.admixture()->value();
    }
}

//-------------------------
inline std::unique_ptr<Model> d4pi()
{
    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/d4pi.pdl");
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
    auto a_1_shape = d4pi_a1_bowler ? std::make_shared<BowlerMassShape>(T["a_1+"]) : std::make_shared<ConstantWidthBreitWigner>(T["a_1+"]);
    if (d4pi_a1_bowler) {
        //a_1_shape->width()->setValue(0.560);
        a_1_shape->width()->setValue(0.459);
    }

    auto a_1_plus  = DecayingParticle::create(T["a_1+"], r, a_1_shape);
    auto a_1_minus = DecayingParticle::create(T["a_1-"], r, a_1_shape);


    // pi(1300)
    auto pi_1300_shape = std::make_shared<ConstantWidthBreitWigner>(T["pi(1300)+"]);
    auto pi_1300_plus = DecayingParticle::create(T["pi(1300)+"], r, pi_1300_shape);
    auto pi_1300_minus = DecayingParticle::create(T["pi(1300)-"], r, pi_1300_shape);

    // pi(1800)
    auto pi_1800_shape = std::make_shared<ConstantWidthBreitWigner>(T["pi(1800)+"]);
    auto pi_1800_plus = DecayingParticle::create(T["pi(1800)+"], r, pi_1800_shape);
    auto pi_1800_minus = DecayingParticle::create(T["pi(1800)-"], r, pi_1800_shape);

    // a_1(1420)
    auto a_1_1420_shape = std::make_shared<ConstantWidthBreitWigner>(T["a_1(1420)+"]);
    auto a_1_1420_plus = DecayingParticle::create(T["a_1(1420)+"], r, a_1_1420_shape);
    auto a_1_1420_minus = DecayingParticle::create(T["a_1(1420)-"], r, a_1_1420_shape);

    // a_1(1640)
    auto a_1_1640_shape = std::make_shared<ConstantWidthBreitWigner>(T["a_1(1640)+"]);
    auto a_1_1640_plus = DecayingParticle::create(T["a_1(1640)+"], r, a_1_1640_shape);
    auto a_1_1640_minus = DecayingParticle::create(T["a_1(1640)-"], r, a_1_1640_shape);

    // a_2(1320)
    auto a_2_1320_shape = std::make_shared<ConstantWidthBreitWigner>(T["a_2(1320)+"]);
    auto a_2_1320_plus = DecayingParticle::create(T["a_2(1320)+"], r, a_2_1320_shape);
    auto a_2_1320_minus = DecayingParticle::create(T["a_2(1320)-"], r, a_2_1320_shape);

    // pi_2(1670)
    auto pi_2_1670_shape = std::make_shared<ConstantWidthBreitWigner>(T["pi_2(1670)+"]);
    auto pi_2_1670_plus = DecayingParticle::create(T["pi_2(1670)+"], r, pi_2_1670_shape);
    auto pi_2_1670_minus = DecayingParticle::create(T["pi_2(1670)-"], r, pi_2_1670_shape);

    // rho
    auto rho = DecayingParticle::create(T["rho0"], r, std::make_shared<GounarisSakurai>(T["rho0"]));
    rho->addStrongDecay(piPlus, piMinus);

    // omega
    auto omega = DecayingParticle::create(T["omega"], r, std::make_shared<BreitWigner>(T["omega"]));
    omega->addStrongDecay({piPlus, piMinus});

    // sigma / f_0(500)
    //auto sigma = DecayingParticle::create(T["f_0(500)"], r, std::make_shared<BreitWigner>(T["f_0(500)"]));
    //sigma->addStrongDecay(piPlus, piMinus);

    // (pi pi)S wave
    auto pipiS = DecayingParticle::create("pipiS", QuantumNumbers(0, 0), r, std::make_shared<PiPiSWaveAuMorganPenningtonKachaev>());
    pipiS->addStrongDecay<d4pi_max_L>(piPlus, piMinus);

    auto pipi = DecayingParticle::create("pipi", QuantumNumbers(0, 0), r);
    pipi->addStrongDecay<d4pi_max_L>(piPlus, piMinus);

    // f_0(980) (as Flatte)
    //auto f_0_980_flatte = std::make_shared<Flatte>(T["f_0"]);
    auto f_0_980_flatte = std::make_shared<SmearedFlatte>(T["f_0"], 0.00397333297611, 0.2, 1.6); // sigma from K resolution
    f_0_980_flatte->add(FlatteChannel(0.406 * 0.406,     *piPlus, *piMinus));
    f_0_980_flatte->add(FlatteChannel(0.406 * 0.406 * 4, T["K+"], T["K-"]));
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

    // f_0(1500)
    auto f_0_1500 = DecayingParticle::create(T["f_0(1500)"], r, std::make_shared<BreitWigner>(T["f_0(1500)"]));
    f_0_1500->addStrongDecay({piPlus, piMinus});

    //
    // decays and amplitudes
    //
    if (d4pi_pipiS_pipiS) {
        D->addWeakDecay<d4pi_max_L>(pipiS, pipiS);
        *free_amplitude(*D, to(pipiS, pipiS)) = d4pi_amp_pipiS_pipiS;
    }

    if (d4pi_pipiS_f_0) {
        D->addWeakDecay<d4pi_max_L>(pipiS, f_0_980);
        *free_amplitude(*D, to(pipiS, f_0_980)) = d4pi_amp_pipiS_f_0;
    }

    if (d4pi_pipiS_rho) {
        D->addWeakDecay<d4pi_max_L>(pipiS, rho);
        *free_amplitude(*D, to(pipiS, rho)) = d4pi_amp_pipiS_rho;
    }

    if (d4pi_pipiS_omega) {
        D->addWeakDecay<d4pi_max_L>(pipiS, omega);
        *free_amplitude(*D, to(pipiS, omega)) = d4pi_amp_pipiS_omega;
    }

    if (d4pi_pipiS_f_2) {
        D->addWeakDecay<d4pi_max_L>(pipiS, f_2);
        *free_amplitude(*D, to(pipiS, f_2)) = d4pi_amp_pipiS_f_2;
    }


    if (d4pi_f_0_f_0) {
        D->addWeakDecay<d4pi_max_L>(f_0_980, f_0_980);
        *free_amplitude(*D, to(f_0_980, f_0_980)) = d4pi_amp_f_0_f_0;
    }

    if (d4pi_f_0_rho) {
        D->addWeakDecay<d4pi_max_L>(f_0_980, rho);
        *free_amplitude(*D, to(f_0_980, rho)) = d4pi_amp_f_0_rho;
    }

    if (d4pi_f_0_omega) {
        D->addWeakDecay<d4pi_max_L>(f_0_980, omega);
        *free_amplitude(*D, to(f_0_980, omega)) = d4pi_amp_f_0_omega;
    }

    if (d4pi_f_0_f_2) {
        D->addWeakDecay<d4pi_max_L>(f_0_980, f_2);
        *free_amplitude(*D, to(f_0_980, f_2)) = d4pi_amp_f_0_f_2;
    }


    if (d4pi_rho_rho) {
        D->addWeakDecay<d4pi_max_L>(rho, rho);
        assert(free_amplitudes(*D, to(rho, rho)).size() == 3);

        for (auto& fa : free_amplitudes(*D, to(rho, rho)))
            *fa = static_cast<std::complex<double> >(d4pi_amp_rho_rho[fa->spinAmplitude()->L()]);
    }

    if (d4pi_rho_omega) {
        D->addWeakDecay<d4pi_max_L>(rho, omega);
        //assert(free_amplitudes(*D, to(rho, omega)).size() == 3);

        for (auto& fa : free_amplitudes(*D, to(rho, omega)))
            *fa = static_cast<std::complex<double> >(d4pi_amp_rho_omega[fa->spinAmplitude()->L()]);
    }

    if (d4pi_rho_f_2) {
        D->addWeakDecay<d4pi_max_L>(rho, f_2);
        //assert(free_amplitudes(*D, to(rho, f_2)).size() == 3);

        for (auto& fa : free_amplitudes(*D, to(rho, f_2)))
            *fa = static_cast<std::complex<double> >(d4pi_amp_rho_f_2[fa->spinAmplitude()->L()]);
    }


    if (d4pi_omega_omega) {
        D->addWeakDecay<d4pi_max_L>(omega, omega);
        //assert(free_amplitudes(*D, to(omega, omega)).size() == 3);

        for (auto& fa : free_amplitudes(*D, to(omega, omega)))
            *fa = static_cast<std::complex<double> >(d4pi_amp_omega_omega[fa->spinAmplitude()->L()]);
    }

    if (d4pi_omega_f_2) {
        D->addWeakDecay<d4pi_max_L>(omega, f_2);
        //assert(free_amplitudes(*D, to(omega, f_2)).size() == 3);

        for (auto& fa : free_amplitudes(*D, to(omega, f_2)))
            *fa = static_cast<std::complex<double> >(d4pi_amp_omega_f_2[fa->spinAmplitude()->L()]);
    }


    if (d4pi_f_2_f_2) {
        D->addWeakDecay<0>(f_2, f_2);
        for (auto& fa : free_amplitudes(*D, to(f_2, f_2)))
            *fa = static_cast<std::complex<double> >(d4pi_amp_f_2_f_2[fa->spinAmplitude()->L()]);
    }


    if (d4pi_f_0_1370_pipiS) {
        D->addWeakDecay<d4pi_max_L>(f_0_1370, pipiS);
        *free_amplitude(*D, to(f_0_1370, pipiS)) = d4pi_amp_f_0_1370_pipiS;
    }


    // a_1
    /*if (d4pi_a_rho_pi_S or d4pi_a_rho_pi_D) {
        a_1_plus->addStrongDecay(rho, piPlus);
        a_1_minus->addStrongDecay(rho, piMinus);

        if (d4pi_a_rho_pi_S) {
            *free_amplitude(*a_1_plus, to(rho), l_equals(0)) = 1.; // reference amplitude
            if (d4pi_fix_a_rho_pi_S)
                free_amplitude(*a_1_plus, to(rho), l_equals(0))->variableStatus() = VariableStatus::fixed;

            free_amplitude(*a_1_minus, to(rho), l_equals(0))->shareFreeAmplitude(*free_amplitude(*a_1_plus, to(rho), l_equals(0)));
        }
        else {
            *free_amplitude(*a_1_plus, to(rho), l_equals(0)) = 0.;
            free_amplitude(*a_1_plus, to(rho), l_equals(0))->variableStatus() = VariableStatus::fixed;

            *free_amplitude(*a_1_minus, to(rho), l_equals(0)) = 0.;
            free_amplitude(*a_1_minus, to(rho), l_equals(0))->variableStatus() = VariableStatus::fixed;
        }

        assert(free_amplitudes(*a_1_plus, to(rho), l_equals(1)).empty());
        assert(free_amplitudes(*a_1_minus, to(rho), l_equals(1)).empty());

        if (d4pi_a_rho_pi_D) {
            *free_amplitude(*a_1_plus, to(rho), l_equals(2)) = d4pi_amp_a_rho_pi_D;

            free_amplitude(*a_1_minus, to(rho), l_equals(2))->shareFreeAmplitude(*free_amplitude(*a_1_plus, to(rho), l_equals(2)));
        }
        else {
            *free_amplitude(*a_1_plus, to(rho), l_equals(2)) = 0.;
            free_amplitude(*a_1_plus, to(rho), l_equals(2))->variableStatus() = VariableStatus::fixed;

            *free_amplitude(*a_1_minus, to(rho), l_equals(2)) = 0.;
            free_amplitude(*a_1_minus, to(rho), l_equals(2))->variableStatus() = VariableStatus::fixed;
        }

    }
    if (d4pi_a_pipiS_pi) {
        a_1_plus->addStrongDecay(pipiS, piPlus);
        *free_amplitude(*a_1_plus, to(pipiS)) = d4pi_amp_a_pipiS_pi;

        a_1_minus->addStrongDecay(pipiS, piMinus);

        free_amplitude(*a_1_minus, to(pipiS))->shareFreeAmplitude(*free_amplitude(*a_1_plus, to(pipiS)));
    }
    
    if (not free_amplitudes(*a_1_plus).empty()) {
        D->addWeakDecay<d4pi_max_L>(a_1_plus, piMinus);
        if (d4pi_a1_plus)
            *free_amplitude(*D, to(a_1_plus)) = 1.;
        else
            *free_amplitude(*D, to(a_1_plus)) = 0.;
        free_amplitude(*D, to(a_1_plus))->variableStatus() = VariableStatus::fixed;

        if (d4pi_a1_minus) {
            D->addWeakDecay<d4pi_max_L>(a_1_minus, piPlus);
            if (d4pi_pm_shared)
                free_amplitude(*D, to(a_1_minus))->shareFreeAmplitude(*free_amplitude(*D, to(a_1_plus)));
            else
                *free_amplitude(*D, to(a_1_minus)) = d4pi_amp_a_minus;
        }
    }*/

    if (d4pi_a_1) {
        add_3pi_decays(D, piPlus, piMinus,
                a_1_plus, a_1_minus,
                1, d4pi_amp_a_minus,
                {rho, pipiS, f_2}, // kein f0
                {{1, 0, d4pi_amp_a_rho_pi_D}, {0, d4pi_amp_a_pipiS_pi}, {0, 1}});
    }
    free_amplitude(*a_1_plus, to(rho), l_equals(0))->variableStatus() = VariableStatus::fixed;
    free_amplitude(*D, to(a_1_plus))->variableStatus() = VariableStatus::fixed;


    if (d4pi_pi1300) {
        add_3pi_decays(D, piPlus, piMinus,
                pi_1300_plus, pi_1300_minus,
                d4pi_amp_pi1300_plus, d4pi_amp_pi1300_minus,
                {pipiS, rho},
                {{1}, {0, d4pi_amp_pi1300_rho}});
    }

    if (d4pi_pi1800) {
        add_3pi_decays(D, piPlus, piMinus,
                pi_1800_plus, pi_1800_minus,
                d4pi_amp_pi_1800_plus, d4pi_amp_pi_1800_minus,
                {/*pipiS, f_0_980,*/ f_0_1500}, // pipiS, f0 too much overlap with pipiS pipiS
                {{1}/*, {d4pi_amp_pi_1800_f_0}, {d4pi_amp_pi_1800_f_0_1500}*/});
    }

    if (d4pi_a_1_1420) {
        add_3pi_decays(D, piPlus, piMinus,
                a_1_1420_plus, a_1_1420_minus,
                d4pi_amp_a_1_1420_plus, d4pi_amp_a_1_1420_minus,
                {f_0_980},
                {{0, 1}});
    }

    if (d4pi_a_1_1640) {
        add_3pi_decays(D, piPlus, piMinus,
                a_1_1640_plus, a_1_1640_minus,
                d4pi_amp_a_1_1640_plus, d4pi_amp_a_1_1640_minus,
                {/*pipiS,*/ rho, f_2}, // no f0; no f_2 in data
                {/*{0, 1},*/ {1, 0, 1}, {0, 1}});
    }

    if (d4pi_a_2_1320) {
        add_3pi_decays(D, piPlus, piMinus,
                a_2_1320_plus, a_2_1320_minus,
                d4pi_amp_a_2_1320_plus, d4pi_amp_a_2_1320_minus,
                {rho}, // f2 vernachlaessigbar
                {{0, 0, 1}});
    }

    if (d4pi_pi_2_1670) {
        add_3pi_decays(D, piPlus, piMinus,
                pi_2_1670_plus, pi_2_1670_minus,
                d4pi_amp_pi_2_1670_plus, d4pi_amp_pi_2_1670_minus,
                {f_2, rho, pipiS},
                {{1, 0, d4pi_amp_f_2_D}, {0, d4pi_amp_rho}, {0, 0, d4pi_amp_pipiS}});
    }

    //
    // background
    //
    std::vector<double> admixtures;

    if (d4pi_bg_flat_4pi) {
        auto flat_4pi = DecayingParticle::create("flat_4pi", QuantumNumbers(0, 0), r);
        flat_4pi->addWeakDecay<d4pi_max_L>(piPlus, piMinus, piPlus, piMinus);
        M->addInitialState(flat_4pi);
        admixtures.push_back(d4pi_admixture_bg_flat_4pi);
    }

    /*if (d4pi_bg_pipiS) {
        M->addInitialState(pipiS);
        admixtures.push_back(1.);
    }*/

    if (d4pi_bg_rho) {
        // nonresonant, just want the mass shape
        auto rho_bg = DecayingParticle::create("rho_bg", QuantumNumbers(0, 0), r, rho->massShape());
        rho_bg->addStrongDecay<d4pi_max_L>(piPlus, piMinus);

        M->addInitialState(rho_bg);
        admixtures.push_back(d4pi_admixture_bg_rho);
    }


    /*if (d4pi_bg_omega) {
        // nonresonant, just want the mass shape
        auto omega_bg = DecayingParticle::create("omega_bg", QuantumNumbers(0, 0), r, omega->massShape());
        omega_bg->addStrongDecay<d4pi_max_L>(piPlus, piMinus);

        M->addInitialState(omega_bg);
        admixtures.push_back(1.);
    }*/

    if (d4pi_bg_f_0) {
        M->addInitialState(f_0_980);
        admixtures.push_back(1.);
    }

    if (d4pi_bg_a_1) {
        M->addInitialState(a_1_plus);
        M->addInitialState(a_1_minus);
        admixtures.push_back(1.);
        admixtures.push_back(1.);
    }

    if (d4pi_bg_K) {
        auto K0 = DecayingParticle::create(T["K0"], r, std::make_shared<BreitWigner>(T["K0"].mass(), 0.1));
        K0->addStrongDecay<d4pi_max_L>(piPlus, piMinus);
        M->addInitialState(K0);
        admixtures.push_back(1.);
    }

    if (d4pi_bg_sigma) {
        auto sigma = DecayingParticle::create(T["f_0(500)"], r, std::make_shared<PoleMass>(std::complex<double>(0.26012, -0.116053)));
        sigma->addStrongDecay<d4pi_max_L>(piPlus, piMinus);
        M->addInitialState(sigma);
        admixtures.push_back(d4pi_admixture_bg_sigma);
    }

    // lock before dealing with admixtures
    M->lock();

    // loop over admixtures
    unsigned i_adm(0);
    bool D0component(false);
    bool fixed(false);
    for (auto& comp : M->components()) {
        if (comp.particle() == &*D) {
            D0component = true;
            if (d4pi_bg_only) {
                comp.admixture()->setValue(0);
                comp.admixture()->variableStatus() = yap::VariableStatus::fixed;
                LOG(INFO) << "fixed D0 admixture to 0";
            }
            if (not d4pi_fix_bg) {
                comp.admixture()->variableStatus() = yap::VariableStatus::fixed;
                LOG(INFO) << "fixed D0 admixture";
            }
        }
        else if (comp.admixture()->variableStatus() != yap::VariableStatus::fixed) {

            //if (pars.empty())
            *comp.admixture() = admixtures[i_adm++];

            if (d4pi_fix_bg) {
                comp.admixture()->variableStatus() = yap::VariableStatus::fixed;
                LOG(INFO) << "fixed bg admixture";
            }
        }

        if (comp.admixture()->variableStatus() == yap::VariableStatus::fixed)
            fixed = true;
    }

    if (not D0component and not fixed and not d4pi_fix_amplitudes) {
        //M->components()[0].admixture()->setValue(1.);
        M->components()[0].admixture()->variableStatus() = yap::VariableStatus::fixed;
    }

    if (d4pi_fix_amplitudes)
        fix_free_amplitudes(*M);

    // print summary
    LOG(INFO) << "D Decay trees:";
    LOG(INFO) << to_string(D->decayTrees());

    printAmplitudes(*M);

    return M;
}

#endif
