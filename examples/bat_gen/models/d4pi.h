// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D4PI__H
#define __BAT__D4PI__H

#include "../bat_fit.h"
#include "../fit_fitFraction.h"
#include "../tools.h"
#include "../../../data/set_parities.h"
#include "../../../data/deduce_parities.h"

#include "d4pi_scales.h"

#include <AmplitudeBasis.h>
#include <Attributes.h>
#include <BowlerMassShape.h>
#include <BreitWigner.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MathUtilities.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <SpinAmplitudeCache.h>

#include <DTo4piStructs.h>

#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCConstantPrior.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <TDirectory.h>
#include <TEventList.h>
#include <TTree.h>

#include <assert.h>
#include <complex>
#include <memory>

using namespace yap;

//-------------------------
inline size_t load_data_4pi(yap::DataSet& data, TTree& t, int N, const double BDT_cut,
                     const double K0_cut = 0, // cut away events with pi+pi- mass within +-K0_cut around K0 mass
                     const bool mc_phasespace_only = false) // only load phasespace MC events
{
    const auto old_size = data.size();
    const auto n_entries = t.GetEntries();

    // get event list (apply BDT cut), speeds up loading
    t.Draw(">>eventList", ("BDT >= " + std::to_string(BDT_cut)).c_str(), "");
    auto eventList = (TEventList*)gDirectory->Get("eventList");
    LOG(INFO) << eventList->GetN() << " events after BDT cut";
    // apply phasespace cut
    if (mc_phasespace_only) {
        t.SetEventList(eventList);
        t.Draw(">>eventListPHSP", "mc_good_D0 && mc_good_DS && mc_direct_4pi", "");
        eventList = (TEventList*)gDirectory->Get("eventListPHSP");
        LOG(INFO) << eventList->GetN() << " events after phasespace cut";
    }
    auto eventListArray = eventList->GetList();

    // set branch addresses
    EVENT* E = nullptr;
    t.SetBranchAddress("E", &E);
    t.SetBranchStatus("*", 0);
    t.SetBranchStatus("mom_uc_pi*", 1);
    t.SetBranchStatus("mom_pi*", 1);

    if (N <= 0) // attempt to load all data
        N = eventList->GetN();
    else
        N = std::min(N, eventList->GetN());

    data.reserve(old_size + N);

    int n_loaded = 0;

    std::vector<yap::FourVector<double>> P; // pi+ pi- pi+ pi-
    std::vector<double> masses;

    // get K0 mass
    const double K0_mass = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl")["K0"].mass();

    for (long i=0; i<eventList->GetN(); ++i) {
        unsigned long long n = eventListArray[i];

        if (t.GetEntry(n) < 1)
            continue;

        // K0 cut with UNconstrained 4 momenta
        if (K0_cut > 0) {
            masses.clear();
            masses.push_back((E->mom_uc_piPlus1 + E->mom_uc_piMinus1).M());
            masses.push_back((E->mom_uc_piPlus1 + E->mom_uc_piMinus2).M());
            masses.push_back((E->mom_uc_piPlus2 + E->mom_uc_piMinus1).M());
            masses.push_back((E->mom_uc_piPlus2 + E->mom_uc_piMinus2).M());

            bool cut(false);
            for (auto m : masses)
                if (fabs(m - K0_mass) < K0_cut) {
                    cut = true;
                    break;
                }

            if (cut)
                continue;
        }

        // Fill CONSTRAINED 4-momenta
        P.clear();
        // \todo ATTENTION: I had a bug in my analysis,
        // which means piPlus2 and piMinus1 (constrained) are swapped
        P.push_back(convert(E->mom_piPlus1));
        P.push_back(convert(E->mom_piPlus2)); // actually mom_piMinus1
        P.push_back(convert(E->mom_piMinus1));// actually mom_piPlus2
        P.push_back(convert(E->mom_piMinus2));

        data.push_back(P);

        if (++n_loaded >= N)
            break;

        if (n_loaded%10 == 0)
            std::cout << "\r" << 100. * n_loaded/N << "%                           ";
    }
    std::cout << std::endl;

    if (data.size() == old_size)
        LOG(INFO) << "No data loaded.";
    else {
        LOG(INFO) << "Loaded " << data.size() - old_size << " data points (" << ((data.size() - old_size) * data[0].bytes() * 1.e-6) << " MB)"
                << " from a tree of size " << n_entries;
        if (old_size)
            LOG(INFO) << "Total data size now " << data.size() << " points (" << (data.bytes() * 1.e-6) << " MB)";
    }

    if (int(data.size() - old_size) < N)
        LOG(WARNING) << "could not load as many data points as requested.";

    return data.size() - old_size;
}

//-------------------------
inline std::unique_ptr<Model> d4pi()
{
    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");
    deduce_meson_parities(T);
    set_parities(T);

    // final state particles
    auto piPlus  = FinalStateParticle::create(T[211]);
    auto piMinus = FinalStateParticle::create(T[-211]);

    auto M = std::make_unique<yap::Model>(std::make_unique<yap::HelicityFormalism>());
    M->setFinalState(piPlus, piMinus, piPlus, piMinus);

    // use common radial size for all resonances
    double r = 1.2; // [GeV^-1]

    // initial state particle
    auto D = DecayingParticle::create(T["D0"], r);
    
    //
    // resonances
    //
    // a_1
    auto a_1 = a1_bowler ? DecayingParticle::create(T["a_1+"], r, std::make_shared<BowlerMassShape>(T["a_1+"])) :
                           DecayingParticle::create(T["a_1+"], r, std::make_shared<BreitWigner>(T["a_1+"]));
    if (a1_bowler)
        //std::dynamic_pointer_cast<BowlerMassShape>(a_1->massShape())->width()->setValue(0.560);
        std::dynamic_pointer_cast<BowlerMassShape>(a_1->massShape())->width()->setValue(0.430);

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

    // pi+ pi- flat
    auto pipiFlat = DecayingParticle::create("pipiFlat", QuantumNumbers(0, 0), r);
    pipiFlat->addStrongDecay(piPlus, piMinus);

    //
    // decays and amplitudes
    //
    if (f_0_pipi) {
        D->addWeakDecay(f_0_980, piPlus, piMinus);
        *free_amplitude(*D, to(f_0_980, piPlus, piMinus)) = amp_f_0_pipi;
    }

    if (f_2_pipi) {
        D->addWeakDecay(f_2, pipiFlat);
        *free_amplitude(*D, to(f_2, pipiFlat)) = amp_f_2_pipi;
    }

    if (sigma_pipi) {
        D->addWeakDecay(sigma, piPlus, piMinus);
        *free_amplitude(*D, to(sigma, piPlus, piMinus)) = amp_sigma_pipi;
    }

    if (flat_4pi) {
        D->addWeakDecay(pipiFlat, pipiFlat);
        *free_amplitude(*D, to(pipiFlat, pipiFlat)) = amp_pipiFlat;
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
        a_1->addStrongDecay(rho, piPlus);

        if (a_rho_pi_S) {
            *free_amplitude(*a_1, to(rho), l_equals(0)) = 1.; // reference amplitude
            free_amplitude(*a_1, to(rho), l_equals(0))->variableStatus() = VariableStatus::fixed;
        }
        else {
            *free_amplitude(*a_1, to(rho), l_equals(0)) = 0.;
            free_amplitude(*a_1, to(rho), l_equals(0))->variableStatus() = VariableStatus::fixed;
        }

        assert(free_amplitudes(*a_1, to(rho), l_equals(1)).empty());

        if (a_rho_pi_D) {
            *free_amplitude(*a_1, to(rho), l_equals(2)) = amp_a_rho_pi_D;
        }
        else {
            *free_amplitude(*a_1, to(rho), l_equals(2)) = 0.;
            free_amplitude(*a_1, to(rho), l_equals(2))->variableStatus() = VariableStatus::fixed;
        }

        if (omega_omega) {
            a_1->addStrongDecay(omega, piPlus);

            if (a_rho_pi_S) {
                *free_amplitude(*a_1, to(omega), l_equals(0)) = amp_a_omega_pi_S;
            }
            else {
                *free_amplitude(*a_1, to(omega), l_equals(0)) = 0.;
                free_amplitude(*a_1, to(omega), l_equals(0))->variableStatus() = VariableStatus::fixed;
            }

            assert(free_amplitudes(*a_1, to(omega), l_equals(1)).empty());

            if (a_rho_pi_D)
                *free_amplitude(*a_1, to(omega), l_equals(2)) = amp_a_rho_pi_D;
            else {
                *free_amplitude(*a_1, to(omega), l_equals(2)) = 0.;
                free_amplitude(*a_1, to(omega), l_equals(2))->variableStatus() = VariableStatus::fixed;
            }
        }
    }
    if (a_sigma_pi) {
        a_1->addStrongDecay(sigma, piPlus);
        *free_amplitude(*a_1, to(sigma)) = amp_a_sigma_pi;
    }
    
    if (not free_amplitudes(*a_1).empty()) {
        D->addWeakDecay(a_1, piMinus);
        *free_amplitude(*D, to(a_1)) = 1.;
        free_amplitude(*D, to(a_1))->variableStatus() = VariableStatus::fixed;
    }

    
    //
    // background
    //
    if (bg_flat_4pi) {
        auto flat_4pi = DecayingParticle::create("flat_4pi", QuantumNumbers(0, 0), r);
        flat_4pi->addWeakDecay(pipiFlat, pipiFlat);
        M->addInitialState(flat_4pi);
    }

    if (bg_rho_2pi) {
        auto rho_2pi = DecayingParticle::create("rho_2pi", QuantumNumbers(0, 0), r);
        rho_2pi->addWeakDecay(rho, pipiFlat);
        M->addInitialState(rho_2pi);
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
            << "  \t" << fa.get();

    // loop over admixtures
    for (auto& comp : M->components()) {
        LOG(INFO) << "admixture " << comp.particle()->name()
                << "; m = " << spin_to_string(comp.decayTrees()[0]->initialTwoM())
                << "; variableStatus = " << to_string(comp.admixture()->variableStatus());
        if (comp.particle() == &*D) {
            comp.admixture()->variableStatus() = yap::VariableStatus::fixed;
            assert( comp.admixture()->variableStatus() == yap::VariableStatus::fixed );
            LOG(INFO) << "fixed D0 admixture";
        }
        else if (comp.admixture()->variableStatus() != yap::VariableStatus::fixed)
            *comp.admixture() = 0.1;
    }

    return M;
}

//-------------------------
inline bat_fit d4pi_fit(std::string name, std::vector<std::vector<unsigned> > pcs = {})
{
    bat_fit m(name, d4pi(), pcs);

    LOG(INFO) << "setting priors";
    // set priors: complete range in real and imag
    /*for (const auto& fa : m.freeAmplitudes()) {
        double re = real(fa->value());
        double im = imag(fa->value());

        double rangeHi = 2.0;

        m.setRealImagRanges(fa, -rangeHi*fabs(re), rangeHi*fabs(re), -rangeHi*fabs(im), rangeHi*fabs(im));
    }*/

    // set priors: range around expected value
    for (const auto& fa : m.freeAmplitudes()) {
        double re = real(fa->value());
        double im = imag(fa->value());

        double ab = abs(fa->value());
        double ar = deg(arg(fa->value()));
        double rangeLo = -100.;
        double rangeHi = 102.;
        double rangeMin = 0.4;
        m.setPriors(fa, new BCConstantPrior(std::max(0., rangeLo*ab), rangeHi*ab),
                new BCConstantPrior(ar - 180, ar + 180));
        m.setRealImagRanges(fa, std::min(re-rangeMin, std::min(rangeLo*re, rangeHi*re)), std::max(re+rangeMin, std::max(rangeLo*re, rangeHi*re)),
                                std::min(im-rangeMin, std::min(rangeLo*im, rangeHi*im)), std::max(im+rangeMin, std::max(rangeLo*im, rangeHi*im)));
        m.setAbsArgRanges(fa, std::max(0., rangeLo*ab), rangeHi*ab,
                          ar - 180, ar + 180);
    }




    // free a_1 width
    auto a_1   = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("a_1+")));
    //m.addParameter("width(a_1)", std::dynamic_pointer_cast<BreitWigner>(a_1->massShape())->width(), 0, 1.4);
    //m.GetParameters().Back().SetPriorConstant();

    if (a1_bowler) {
        //m.addParameter("K*K_coupling", std::dynamic_pointer_cast<BowlerMassShape>(a_1->massShape())->coupling(), 0, 1.);
        //m.GetParameters().Back().SetPriorConstant();
    }

    return m;
}


//-------------------------
inline void d4pi_printFitFractions(bat_fit& m)
{
    LOG(INFO) << "\nFit fractions for single decay trees:";
    double sum(0);
    for (const auto& mci : m.modelIntegral().integrals()) {
        auto ff = fit_fractions(mci.Integral);
        for (size_t i = 0; i < ff.size(); ++i) {
            LOG(INFO) << to_string(*mci.Integral.decayTrees()[i]) << "\t" << ff[i].value()*100. << " %";
            sum += ff[i].value();
        }
    }
    LOG(INFO) << "Sum = " << sum*100 << " %";

    LOG(INFO) << "\nFit fractions for grouped decay trees:";
    sum = 0;
    for (const auto& mci : m.modelIntegral().integrals()) {

        // particles
        auto rho   = rho_rho     ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("rho0"))) : nullptr;
        auto omega = omega_omega ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("omega"))) : nullptr;
        auto sigma = (a_sigma_pi or sigma_pipi) ?  std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(500)"))) : nullptr;
        auto a_1   = (a_rho_pi_S or a_rho_pi_D or a_sigma_pi) ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("a_1+"))) : nullptr;
        auto f_0   = f_0_pipi   ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0"))) : nullptr;
        auto f_2   = f_2_pipi   ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_2"))) : nullptr;
        auto pipiFlat = flat_4pi ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("pipiFlat"))) : nullptr;

        auto decayTrees = mci.Integral.decayTrees();
        std::vector<yap::DecayTreeVector> groupedDecayTrees;


        // a_rho_pi_S
        if (a_rho_pi_S) {
            auto blub = yap::filter(decayTrees, yap::to(a_1));
            blub.erase(std::remove_if(blub.begin(), blub.end(),
                    [&](const std::shared_ptr<yap::DecayTree>& dt){return (filter(dt->daughterDecayTreeVector(), to(rho), l_equals(0))).empty() and (filter(dt->daughterDecayTreeVector(), to(omega), l_equals(0))).empty();}),
                    blub.end());
            groupedDecayTrees.push_back(blub);

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }

        // a_rho_pi_D
        if (a_rho_pi_D) {
            auto blub = yap::filter(decayTrees, yap::to(a_1));
            blub.erase(std::remove_if(blub.begin(), blub.end(),
                    [&](const std::shared_ptr<yap::DecayTree>& dt){return (filter(dt->daughterDecayTreeVector(), to(rho), l_equals(2))).empty() and (filter(dt->daughterDecayTreeVector(), to(omega), l_equals(2))).empty();}),
                    blub.end());
            groupedDecayTrees.push_back(blub);

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }

        // a_sigma_pi
        if (a_sigma_pi) {
            auto blub = yap::filter(decayTrees, yap::to(a_1));
            blub.erase(std::remove_if(blub.begin(), blub.end(),
                    [&](const std::shared_ptr<yap::DecayTree>& dt){return (filter(dt->daughterDecayTreeVector(), yap::to(sigma))).empty();}),
                    blub.end());
            groupedDecayTrees.push_back(blub);

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }

        // rho_rho
        // omega_omega
        if (rho_rho or omega_omega) {
            DecayTreeVector rhos = yap::filter(decayTrees, yap::to(rho));
            if (omega_omega) {
                DecayTreeVector omegas = yap::filter(decayTrees, yap::to(omega));
                rhos.insert(rhos.end(), omegas.begin(), omegas.end());
            }
            groupedDecayTrees.push_back(rhos);

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }

        // f_0_pipi
        if (f_0_pipi) {
            groupedDecayTrees.push_back(yap::filter(decayTrees, yap::to(f_0)));

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }


        // f_2_pipi
        if (f_2_pipi) {
            groupedDecayTrees.push_back(yap::filter(decayTrees, yap::to(f_2)));

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }

        // sigma_pipi
        if (sigma_pipi) {
            groupedDecayTrees.push_back(yap::filter(decayTrees, yap::to(sigma)));

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }

        if (flat_4pi)
            groupedDecayTrees.push_back(yap::filter(decayTrees, yap::to(pipiFlat)));


        auto ff = fit_fractions(mci.Integral, groupedDecayTrees);
        for (size_t i = 0; i < ff.size(); ++i) {

            LOG(INFO) << "" ;
            LOG(INFO) << yap::daughter_name()(groupedDecayTrees[i][0])
                      << "   L = " << yap::orbital_angular_momentum()(groupedDecayTrees[i][0])
                      << " :: " << std::to_string(groupedDecayTrees[i].size()) << " amplitudes:";
            for (const auto& fa : groupedDecayTrees[i])
                LOG(INFO) << yap::to_string(*fa);

            LOG(INFO) << "\t" << ff[i].value()*100. << " %";
            sum += ff[i].value();
        }
        LOG(INFO) << "Sum = " << sum*100 << " %\n";


        // loop over admixtures
        for (auto& comp : m.model()->components()) {
            LOG(INFO) << "admixture " << comp.particle()->name()
                    << "; m = " << spin_to_string(comp.decayTrees()[0]->initialTwoM())
                    << "; variableStatus = " << to_string(comp.admixture()->variableStatus())
                    << "; \t value = " << comp.admixture()->value();
        }


        if (a_rho_pi_S and a_rho_pi_D and a_sigma_pi and rho_rho and f_0_pipi and f_2_pipi and sigma_pipi) {
            // scaling to match FOCUS
            LOG(INFO) << "\nnew scales to match FOCUS model";
            std::vector<double> focusFractions{43.3, 2.5, 8.3, 24.5, 2.4, 4.9, 8.2};
            std::vector<double> currentScales{1.0, scale_a_rho_pi_D, scale_a_sigma_pi, scale_rho_rho, scale_f_0_pipi, scale_f_2_pipi, scale_sigma_pipi};
            double fix = sqrt(focusFractions[0] / ff[0].value()); // a_rho_pi_S
            //assert(ff.size() == focusFractions.size() and ff.size() == currentScales.size());
            for (size_t i = 1; i < ff.size(); ++i) {
                LOG(INFO) << currentScales[i] * sqrt(focusFractions[i] / ff[i].value()) / fix;
            }
        }

    }
}

#endif
