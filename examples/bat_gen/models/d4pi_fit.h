// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D4PI_FIT__H
#define __BAT__D4PI_FIT__H

#include "../bat_fit.h"
#include "../fit_fitFraction.h"
#include "../tools.h"

#include "d4pi.h"

#include <DTo4piStructs.h>

//#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCConstantPrior.h>
//#include <BAT/BCLog.h>
//#include <BAT/BCParameterSet.h>

#include <TDirectory.h>
#include <TEventList.h>
#include <TFile.h>
#include <TTree.h>
#include <TKDTreeBinning.h>
#include <TH2D.h>

#include <assert.h>
#include <complex>
#include <future>
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
    if (BDT_cut >= 0)
        t.Draw(">>eventList", ("BDT >= " + std::to_string(BDT_cut)).c_str(), "");
    else
        t.Draw(">>eventList", ("BDT < " + std::to_string(BDT_cut)).c_str(), "");
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
    const double K0_mass = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/d4pi.pdl")["K0"].mass();
    LOG(INFO) << "K0 mass = " << K0_mass;

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
inline void generate_fit_fraction_data(bat_fit& m, unsigned nPoints, unsigned nPartitions)
{
    auto T = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    //assert(m.model()->initialStates()[0]->name() == "D0");
    auto isp_mass = T["D0"].mass();

    auto A = m.model()->massAxes();
    auto m2r = yap::squared(yap::mass_range(isp_mass, A, m.model()->finalStateParticles()));

    std::mt19937 g(0);
    // fill data set with nPoints points
    std::generate_n(std::back_inserter(m.fitFractionData()), nPoints,
            std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), isp_mass, A, m2r, g, std::numeric_limits<unsigned>::max()));
    m.fitFractionPartitions() = yap::DataPartitionBlock::create(m.fitFractionData(), nPartitions);

    yap::ImportanceSampler::calculate(m.fitFractionIntegral(), m.fitFractionPartitions(), true);
}

//-------------------------
inline bat_fit d4pi_fit(std::string name, std::vector<std::vector<unsigned> > pcs = {})
{
    bat_fit m(name, d4pi(), pcs);

    //LOG(INFO) << "setting priors";
    // set priors: complete range in real and imag
    /*for (const auto& fa : m.freeAmplitudes()) {
        double re = real(fa->value());
        double im = imag(fa->value());

        double rangeHi = 2.0;

        m.setRealImagRanges(fa, -rangeHi*fabs(re), rangeHi*fabs(re), -rangeHi*fabs(im), rangeHi*fabs(im));
    }*/

    // set priors: range around expected value
    for (const auto& fa : m.freeAmplitudes()) {
        //double re = real(fa->value());
        //double im = imag(fa->value());

        double ab = abs(fa->value());
        double ar = deg(arg(fa->value()));
        double range = 2.;
        m.setPriors(fa, new BCConstantPrior(0, range*ab),
                new BCConstantPrior(ar - 180, ar + 180));
        /*m.setRealImagRanges(fa, std::min(re-rangeMin, std::min(rangeLo*re, rangeHi*re)), std::max(re+rangeMin, std::max(rangeLo*re, rangeHi*re)),
                                std::min(im-rangeMin, std::min(rangeLo*im, rangeHi*im)), std::max(im+rangeMin, std::max(rangeLo*im, rangeHi*im)));*/
        m.setRealImagRanges(fa, -range*ab, range*ab, -range*ab, range*ab);
        m.setAbsArgRanges(fa, 0, range*ab,
                          ar - 180, ar + 180);
    }




    // free a_1 width
    //auto a_1   = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("a_1+")));
    //m.addParameter("width(a_1)", std::dynamic_pointer_cast<BreitWigner>(a_1->massShape())->width(), 0, 1.4);
    //m.GetParameters().Back().SetPriorConstant();

    if (d4pi_free_parameters) {
        /*if (a_sigma_pi or sigma_f_0_1370) {
            auto sigma  = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(500)")));
            auto width = std::dynamic_pointer_cast<BreitWigner>(sigma->massShape())->width();
            m.addParameter("width(f_0(500))", width, 0.5*width->value(), 1.5*width->value());
            m.GetParameters().Back().SetPriorConstant();
        }*/

        if (d4pi_f_0_1370_pipiS) {
            auto f_0_1370  = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(1370)")));
            auto width = std::dynamic_pointer_cast<BreitWigner>(f_0_1370->massShape())->width();
            m.addParameter("width(f_0(1370))", width, 0.5*width->value(), 1.5*width->value());
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.35, 0.15));
        }

        if (d4pi_a1_bowler) {
            auto a_1   = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("a_1+")));
            auto shape = std::dynamic_pointer_cast<BowlerMassShape>(a_1->massShape());

            double mass = shape->mass()->value();
            m.addParameter("mass(a_1)", shape->mass(), 0.9*mass, 1.1*mass);
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(1.23, 0.04));

            double width = shape->width()->value();
            m.addParameter("width(a_1)", shape->width(), 0.7*width, 1.3*width);
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.425, 0.175));

            double coupling = shape->coupling()->value();
            m.addParameter("K*K_coupling", shape->coupling(), 0.5*coupling, 2.*coupling);
            m.GetParameters().Back().SetPriorConstant();
        }

        if (d4pi_pi1300) {
            auto pi1300 = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("pi(1300)+")));
            auto shape  = std::dynamic_pointer_cast<ConstantWidthBreitWigner>(pi1300->massShape());

            m.addParameter("mass(pi(1300)+)", shape->mass(), 1.1, 1.5);
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(1.3, 0.1));

            m.addParameter("width(pi(1300)+)", shape->width(), 0.2, 0.6);
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.4, 0.2));
        }

        // f_0(980)
        if (not particles(*m.model(), yap::is_named("f_0")).empty()) {
            auto f_0 = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0")));
            auto f_0_980_flatte = std::dynamic_pointer_cast<Flatte>(f_0->massShape());

            m.addParameter("mass(f_0)", f_0_980_flatte->mass(), 0.97, 1.01);
            for (auto ch : f_0_980_flatte->channels()) {
                double coupling = ch.Coupling->value();
                m.addParameter("f_0_coupling_" + ch.Particles[0]->name(), ch.Coupling, 0.5*coupling, 2*coupling);
            }
        }
    }

    if (d4pi_bg_K) {
        auto K0   = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("K0")));
        auto shape = std::dynamic_pointer_cast<BreitWigner>(K0->massShape());

        m.addParameter("width(K0)", shape->width(), 0, 0.2);
        m.GetParameters().Back().SetPriorConstant();
    }

    if (d4pi_free_parameters and d4pi_bg_sigma) {
        auto sigma   = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(500)")));
        auto shape = std::dynamic_pointer_cast<PoleMass>(sigma->massShape());

        // mass is complex
        m.addParameter("mass(f_0(500))", shape->mass(), std::complex<double>(0.2, -0.7), std::complex<double>(0.3, -0.16));

    }

    return m;
}

//-------------------------
void draw_overlap_integrals(const DecayTreeVectorIntegral& dtvi)
{
    auto& diagonals = dtvi.diagonals();
    auto offDiagonals = dtvi.offDiagonals();

    TFile file("output/overlapIntegrals.root", "RECREATE");
    TH2D histo("overlapIntegrals", "overlapIntegrals", diagonals.size(), 0, diagonals.size(), diagonals.size(), 0, diagonals.size());

    for (unsigned i = 0; i < diagonals.size(); ++i) {
        auto str = to_string(*dtvi.decayTrees().at(i)->freeAmplitude()->decayChannel());
        str.erase(0, 7);
        /*for (const auto& d_dt : dtvi.decayTrees().at(i)->daughterDecayTrees()) {
            if

        }*/

        histo.GetXaxis()->SetBinLabel(i+1, str.c_str());
        histo.GetYaxis()->SetBinLabel(diagonals.size()-i, str.c_str());
    }

    for (unsigned i = 0; i < diagonals.size(); ++i) {
        for (unsigned j = i; j < diagonals.size(); ++j) {
            double val = 0.;
            if (i == j) {
                val = 1.;
            }
            else {
                double norm = 1. / sqrt(diagonals.at(i).value() * diagonals.at(j).value());
                val = norm * abs(offDiagonals.at(i).at(j-i-1).value());
            }
            histo.Fill(j+0.5, diagonals.size()-i-0.5, val);
        }
    }

    histo.Write();
    file.Close();
}

//-------------------------
inline void d4pi_normalizeFitFractions(bat_fit& m)
{
    LOG(INFO) << "\nd4pi_normalizeFitFractions";
    auto groupedDecayTrees = group_by_free_amplitudes(m.fitFractionIntegral().integrals().at(0).Integral.decayTrees());
    auto ff = fit_fractions(m.fitFractionIntegral().integrals().at(0).Integral, groupedDecayTrees);
    for (size_t i = 0; i < ff.size(); ++i) {

        if (groupedDecayTrees[i].empty())
            continue;
        if (!groupedDecayTrees[i][0])
            continue;

        ff = fit_fractions(m.fitFractionIntegral().integrals().at(0).Integral, groupedDecayTrees);

        double frac = ff[i].value();

        std::shared_ptr<FreeAmplitude> fa;
        for (auto ddtv : groupedDecayTrees[i][0]->daughterDecayTreeVector()) {
            fa = ddtv->freeAmplitude();
            if (fa->variableStatus() != yap::VariableStatus::fixed)
                break;
        }

        if (not fa or fa->variableStatus() == yap::VariableStatus::fixed)
            fa =  groupedDecayTrees[i][0]->freeAmplitude();

        if (fa and fa->variableStatus() != yap::VariableStatus::fixed)
            fa->setValue(fa->value() / sqrt(frac*100));


        /*LOG(INFO) << "" ;
        LOG(INFO) << yap::daughter_name()(*(groupedDecayTrees[i][0]))
                  << "   L = " << yap::orbital_angular_momentum()(*(groupedDecayTrees[i][0]))
                  << " :: " << std::to_string(groupedDecayTrees[i].size()) << " amplitudes:";
        for (const auto& dt : groupedDecayTrees[i])
            LOG(INFO) << yap::to_string(*dt);*/

        LOG(INFO) << to_string(*fa) << " = " << fa->value() << "; \t ff =" << frac*100. << " %";

    }
}

//-------------------------
inline void d4pi_printFitFractions(bat_fit& m)
{
    // sum of integrals
    double sumIntegrals(0);
    // loop over admixtures
    for (const auto& mci : m.fitFractionIntegral().integrals()) {
        sumIntegrals += mci.Admixture->value() * integral(mci.Integral).value();
    }
    LOG(INFO) << "\nsum of integrals = " << sumIntegrals;

    double sum(0);

    LOG(INFO) << "\nFit fractions for single decay trees:";
    for (const auto& mci : m.fitFractionIntegral().integrals()) {
        double sum_admixture(0);
        auto ff = fit_fractions(mci.Integral);
        for (size_t i = 0; i < ff.size(); ++i) {
            double fit_frac = mci.Admixture->value()  * integral(mci.Integral).value() / sumIntegrals * ff[i].value();
            LOG(INFO) << to_string(*mci.Integral.decayTrees()[i]) << "\t" << fit_frac*100. << " %; fit fraction in admixture = " << ff[i].value()*100. << " %";
            sum += fit_frac;
            sum_admixture += ff[i].value();
        }
        LOG(INFO) << "Sum for admixture = " << sum_admixture*100 << " %";
        LOG(INFO) << "Deviation from 100% = " << (sum_admixture-1.)*100 << " %";
    }
    LOG(INFO) << "Sum = " << sum*100 << " %";
    LOG(INFO) << "Deviation from 100% = " << (sum-1.)*100 << " %";

    draw_overlap_integrals(m.fitFractionIntegral().integrals()[0].Integral);



    LOG(INFO) << "\nFit fractions grouped by free amplitudes:";
    sum = 0;
    std::multimap<double, unsigned> ff_map;
    auto groupedDecayTrees = group_by_free_amplitudes(m.fitFractionIntegral().integrals().at(0).Integral.decayTrees());
    auto ff = fit_fractions(m.fitFractionIntegral().integrals().at(0).Integral, groupedDecayTrees);
    for (size_t i = 0; i < ff.size(); ++i) {

        if (groupedDecayTrees[i].empty())
            continue;
        if (!groupedDecayTrees[i][0])
            continue;

        LOG(INFO) << "" ;
        LOG(INFO) << yap::daughter_name()(*(groupedDecayTrees[i][0]))
                  << "   L = " << yap::orbital_angular_momentum()(*(groupedDecayTrees[i][0]))
                  << " :: " << std::to_string(groupedDecayTrees[i].size()) << " amplitudes:";
        for (const auto& dt : groupedDecayTrees[i])
            LOG(INFO) << yap::to_string(*dt);

        double frac = ff[i].value();

        LOG(INFO) << "\t" << frac*100. << " %";
        sum += frac;

        ff_map.insert( std::pair<double, unsigned>(-frac, i) );
    }
    LOG(INFO) << "Sum = " << sum*100 << " %\n";


    LOG(INFO) << "\nFit fractions grouped by free amplitudes, sorted by decreasing fit fraction:";
    for (auto kv : ff_map) {
        unsigned i = kv.second;

        LOG(INFO) << "" ;
        LOG(INFO) << yap::daughter_name()(*(groupedDecayTrees[i][0]))
                  << "   L = " << yap::orbital_angular_momentum()(*(groupedDecayTrees[i][0]))
                  << " :: " << std::to_string(groupedDecayTrees[i].size()) << " amplitudes:";
        for (const auto& dt : groupedDecayTrees[i])
            LOG(INFO) << yap::to_string(*dt);
        LOG(INFO) << "\t" << -kv.first*100. << " %";
    }

/*
    LOG(INFO) << "\nFit fractions for grouped decay trees:";
    sum = 0;
    for (const auto& mci : m.fitFractionIntegral().integrals()) {

        // particles
        auto piPlus = std::static_pointer_cast<FinalStateParticle>(particle(*m.model(), is_named("pi+")));
        auto piMinus = std::static_pointer_cast<FinalStateParticle>(particle(*m.model(), is_named("pi-")));

        auto a_1_plus   = (d4pi_a_rho_pi_S or d4pi_a_rho_pi_D or d4pi_a_pipiS_pi) ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("a_1+"))) : nullptr;
        auto a_1_minus  = (a_1_plus and d4pi_a1_minus) ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("a_1-"))) : nullptr;

        auto pi_1300_plus  = particles(*m.model(), is_named("pi+(1300)")).empty() ? nullptr : std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("pi+(1300)")));
        auto pi_1300_minus = particles(*m.model(), is_named("pi-(1300)")).empty() ? nullptr : std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("pi-(1300)")));

        auto pipiS = particles(*m.model(), is_named("pipiS")).empty() ? nullptr : std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("pipiS")));
        auto f_0 = particles(*m.model(), is_named("f_0")).empty() ? nullptr : std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0")));
        auto rho = particles(*m.model(), is_named("rho0")).empty() ? nullptr : std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("rho0")));
        auto omega = particles(*m.model(), is_named("omega")).empty() ? nullptr : std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("omega")));
        auto f_2 = particles(*m.model(), is_named("f_2")).empty() ? nullptr : std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_2")));
        auto f_0_1370 = particles(*m.model(), is_named("f_0(1370)")).empty() ? nullptr : std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(1370)")));


        std::vector<std::shared_ptr<yap::Particle>> twoPiResonances;
        if (pipiS)
            twoPiResonances.push_back(pipiS);
        if (f_0)
            twoPiResonances.push_back(f_0);
        if (rho)
            twoPiResonances.push_back(rho);
        if (omega)
            twoPiResonances.push_back(omega);
        if (f_2)
            twoPiResonances.push_back(f_2);


        auto decayTrees = mci.Integral.decayTrees();
        std::vector<yap::DecayTreeVector> groupedDecayTrees;


        // a_rho_pi_S
        if (d4pi_a_rho_pi_S) {
            auto blub = yap::filter(decayTrees, yap::to(a_1_plus));
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
        if (d4pi_a_rho_pi_D) {
            auto blub = yap::filter(decayTrees, yap::to(a_1_plus));
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
        if (d4pi_a_pipiS_pi) {
            auto blub = yap::filter(decayTrees, yap::to(a_1_plus));
            blub.erase(std::remove_if(blub.begin(), blub.end(),
                    [&](const std::shared_ptr<yap::DecayTree>& dt){return (filter(dt->daughterDecayTreeVector(), yap::to(pipiS))).empty();}),
                    blub.end());
            groupedDecayTrees.push_back(blub);

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }

        // a_rho_pi_S
        if (d4pi_a_rho_pi_S) {
            auto blub = yap::filter(decayTrees, yap::to(a_1_minus));
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
        if (d4pi_a_rho_pi_D) {
            auto blub = yap::filter(decayTrees, yap::to(a_1_minus));
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
        if (d4pi_a_pipiS_pi) {
            auto blub = yap::filter(decayTrees, yap::to(a_1_minus));
            blub.erase(std::remove_if(blub.begin(), blub.end(),
                    [&](const std::shared_ptr<yap::DecayTree>& dt){return (filter(dt->daughterDecayTreeVector(), yap::to(pipiS))).empty();}),
                    blub.end());
            groupedDecayTrees.push_back(blub);

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }

        // pi 1300
        if (pi_1300_plus) {
            DecayTreeVector dtv = yap::filter(decayTrees, yap::to(pi_1300_plus));
            groupedDecayTrees.push_back(dtv);

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }
        if (pi_1300_minus) {
            DecayTreeVector dtv = yap::filter(decayTrees, yap::to(pi_1300_minus));
            groupedDecayTrees.push_back(dtv);

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }

        // 2 pi resonances
        for (unsigned i = 0; i < twoPiResonances.size(); ++i)
            for (unsigned j = i; j < twoPiResonances.size(); ++j) {
                DecayTreeVector dtv = yap::filter(decayTrees, yap::to(twoPiResonances[i], twoPiResonances[j]));
                groupedDecayTrees.push_back(dtv);

                for(auto dt : groupedDecayTrees.back())  {
                    auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                    if(iter != decayTrees.end())
                        decayTrees.erase(iter);
                }
            }

        // rest
        groupedDecayTrees.push_back(decayTrees);

        if (groupedDecayTrees.empty())
            return;

        auto ff = fit_fractions(mci.Integral, groupedDecayTrees);
        for (size_t i = 0; i < ff.size(); ++i) {

            if (groupedDecayTrees[i].empty())
                continue;
            if (!groupedDecayTrees[i][0])
                continue;

            LOG(INFO) << "" ;
            LOG(INFO) << yap::daughter_name()(*(groupedDecayTrees[i][0]))
                      << "   L = " << yap::orbital_angular_momentum()(*(groupedDecayTrees[i][0]))
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

    }*/
}

//-------------------------
inline double chi2_adaptiveBinning(bat_fit& batFit, unsigned EvtsPerBin = 30)
{
    auto massAxes = batFit.model()->massAxes({{0, 1}, {0, 3}, {1, 2}, {2, 3}, {0, 2}});

    unsigned dataDim = 5;
    unsigned dataSize_meas = batFit.fitData().size();
    std::vector<double> data_meas; // stride=1 for the points and = N (the dataSize) for the dimension
    data_meas.reserve(dataDim * dataSize_meas);

    // fill data_meas
    for (auto& indices : massAxes) { // must first loop over mass axes
        for (auto& dp : batFit.fitData()) { // then over data points
            FourVector<double> s;
            for (auto fsp : indices->indices())
                s += batFit.model()->fourMomenta()->finalStateMomenta(dp)[fsp];

            data_meas.push_back(s*s);
        }
    }

    // adaptive binning
    TKDTreeBinning kdBins(dataSize_meas, dataDim, &data_meas[0], dataSize_meas/EvtsPerBin);

    /*std::cout << "bin contents (measured):\n";
    for (unsigned i = 0; i < kdBins.GetNBins(); ++i) {
        std::cout << kdBins.GetBinContent(i) << "  ";
    }
    std::cout << "\n";*/


    // bin integralData
    std::vector<double> histo_exp(kdBins.GetNBins(), 0.);
    std::vector<unsigned> histo_exp_N(kdBins.GetNBins(), 0.);
    std::vector<double> ss(5, 0.); // invariant mass squares

    for (auto& dp : batFit.integralData()) {
        for (unsigned i = 0; i < massAxes.size(); ++i) {
            FourVector<double> s;
            for (auto fsp : massAxes[i]->indices())
                s += batFit.model()->fourMomenta()->finalStateMomenta(dp)[fsp];

            ss[i] = s*s;
        }

        // fill bin with intensity
        histo_exp.at(kdBins.FindBin(&ss[0])) += intensity(*batFit.model(), dp);
        histo_exp_N.at(kdBins.FindBin(&ss[0])) += 1;
    }

    /*std::cout << "bin contents (expected):\n";
    for (unsigned i = 0; i < histo_exp_N.size(); ++i) {
        std::cout << histo_exp_N.at(i) << "  ";
    }
    std::cout << "\n";

    std::cout << "bin contents (expected) * intensity:\n";
    for (unsigned i = 0; i < histo_exp.size(); ++i) {
        std::cout << histo_exp.at(i) << "  ";
    }
    std::cout << "\n";*/


    // calculate chi^2
    double n_exp = std::accumulate(histo_exp.begin(), histo_exp.end(), 0.0);

    double norm_exp = dataSize_meas/n_exp;

    double chi2(0.);
    double NDF = -double(batFit.GetNFreeParameters());

    for (unsigned i = 0; i < histo_exp.size(); ++i) {
        // at least this much events
        if (histo_exp_N.at(i) < 5.) {
            //std::cout<< "x  ";
            continue;
        }

        double NbExp = histo_exp.at(i) * norm_exp;
        double Nb = kdBins.GetBinContent(i);

        double chi2inc = pow(Nb - NbExp, 2) / NbExp;
        //std::cout<< chi2inc << "  ";
        chi2 += chi2inc;
        NDF += 1.;
    }
    std::cout << "\n";

    std::cout << "total chi2 with adaptive binning = " << chi2 << "\n";
    std::cout << "NDF = " << NDF << "\n";
    std::cout << "reduced chi2 = " << chi2/NDF << "\n";

    return chi2;

}

//-------------------------
inline double chi2(bat_fit& batFit, std::vector<double> sbinUpperBounds = {0.237108, 0.549527, 0.734061})
{
    unsigned nBins = sbinUpperBounds.size() + 1;

    //batFit.model()->calculate(batFit.integralData());
    //batFit.model()->calculate(batFit.fitData());

    auto massAxes = batFit.model()->massAxes({{0, 1}, {0, 3}, {1, 2}, {2, 3}, {0, 2}});
    //auto massAxes = batFit.model()->massAxes();

    std::vector<double> ss(5, 0.); // invariant mass squares
    std::vector<unsigned> pointCoords(5, 0); // in 5D space

    std::vector<unsigned> histo_exp_N(pow(nBins, 5), 0.); // predicted, number of events in bins
    std::vector<double> histo_exp(pow(nBins, 5), 0.); // predicted
    std::vector<double> histo_meas(pow(nBins, 5), 0.); // measured

    // phsp data -> predicted events
    for (auto& dp : batFit.integralData()) {
        for (unsigned i = 0; i < massAxes.size(); ++i) {
            FourVector<double> s;
            for (auto fsp : massAxes[i]->indices())
                s += batFit.model()->fourMomenta()->finalStateMomenta(dp)[fsp];

            ss[i] = s*s;
        }

        // determine bin number
        unsigned iCoord(0);
        for (auto& c : ss) {
          for (unsigned i = 0; i < sbinUpperBounds.size(); ++i) {
            if (c < sbinUpperBounds.at(i)) {
              pointCoords.at(iCoord) = i;
              break;
            }
            pointCoords.at(iCoord) = sbinUpperBounds.size();
          }
          ++iCoord;
        }

        unsigned iBin(0);
        for (unsigned i=0; i<pointCoords.size(); ++i) {
          iBin += pointCoords.at(i) * pow(nBins, i);
        }

        // fill bin with intensity
        histo_exp.at(iBin) += intensity(*batFit.model(), dp);
        histo_exp_N.at(iBin) += 1;
    }


    /*std::cout << "histo_exp\n";
    for (auto c : histo_exp) {
      std::cout<< c << "  ";
    }
    std::cout << "\n";*/


    // data -> measured events
    for (auto& dp : batFit.fitData()) {
        for (unsigned i = 0; i < massAxes.size(); ++i) {
            FourVector<double> s;
            for (auto fsp : massAxes[i]->indices())
                s += batFit.model()->fourMomenta()->finalStateMomenta(dp)[fsp];

            ss[i] = s*s;
        }

        // determine bin number
        unsigned iCoord(0);
        for (auto& c : ss) {
          for (unsigned i = 0; i < sbinUpperBounds.size(); ++i) {
            if (c < sbinUpperBounds.at(i)) {
              pointCoords.at(iCoord) = i;
              break;
            }
            pointCoords.at(iCoord) = sbinUpperBounds.size();
          }
          ++iCoord;
        }

        unsigned iBin(0);
        for (unsigned i=0; i<pointCoords.size(); ++i) {
          iBin += pointCoords.at(i) * pow(nBins, i);
        }

        // fill bin with 1
        histo_meas.at(iBin) += 1;
    }


    /*std::cout << "histo_meas\n";
    for (auto c : histo_meas) {
      std::cout<< c << "  ";
    }
    std::cout << "\n";


    std::cout << "chi2 increments\n";*/

    // calculate chi^2
    double n_exp = std::accumulate(histo_exp.begin(), histo_exp.end(), 0.0);
    double n_meas = std::accumulate(histo_meas.begin(), histo_meas.end(), 0.0);

    double norm_exp = n_meas/n_exp;

    double chi2(0.);
    double NDF = -double(batFit.GetNFreeParameters());

    // see arXiv:1611.09253v2 [hep-ex], Eq. 6
    for (unsigned i = 0; i < histo_exp.size(); ++i) {
        double NbExp = histo_exp.at(i) * norm_exp;

        // at least this much events
        if (histo_exp_N.at(i) < 5.) {
            //std::cout<< "x  ";
            continue;
        }

        double Nb = histo_meas.at(i);

        double chi2inc = pow(Nb - NbExp, 2) / NbExp;
        //std::cout<< chi2inc << "  ";
        chi2 += chi2inc;
        NDF += 1.;
    }
    std::cout << "\n";

    std::cout << "total chi2 = " << chi2 << "\n";
    std::cout << "NDF = " << NDF << "\n";
    std::cout << "reduced chi2 = " << chi2/NDF << "\n";

    return chi2;
}


#endif
