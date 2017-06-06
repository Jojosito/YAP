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
//#include <BAT/BCGaussianPrior.h>
#include <BAT/BCConstantPrior.h>
//#include <BAT/BCLog.h>
//#include <BAT/BCParameterSet.h>

#include <TDirectory.h>
#include <TEventList.h>
#include <TTree.h>
#include <TKDTreeBinning.h>

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
    const double K0_mass = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl")["K0"].mass();
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
    //auto a_1   = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("a_1+")));
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
    // generate integration data
    yap::DataSet data(m.model()->createDataSet());

    {
        auto T = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");
        auto isp_mass = T[m.model()->initialStates()[0]->name()].mass();

        auto A = m.model()->massAxes();
        auto m2r = yap::squared(yap::mass_range(isp_mass, A, m.model()->finalStateParticles()));

        std::mt19937 g(0);
        // fill data set with nPoints points
        std::generate_n(std::back_inserter(data), 500000,
                        std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), isp_mass, A, m2r, g, std::numeric_limits<unsigned>::max()));
    }
    auto partitions = yap::DataPartitionBlock::create(data, 4);

    yap::ModelIntegral mi(*m.model());
    yap::ImportanceSampler::calculate(mi, partitions);

    // sum of integrals
    double sumIntegrals(0);
    // loop over admixtures
    for (const auto& mci : mi.integrals()) {
        sumIntegrals += mci.Admixture->value() * integral(mci.Integral).value();
    }
    LOG(INFO) << "\nsum of integrals = " << sumIntegrals;


    LOG(INFO) << "\nFit fractions for single decay trees:";
    double sum(0);
    for (const auto& mci : mi.integrals()) {
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

    LOG(INFO) << "\nFit fractions for grouped decay trees:";
    sum = 0;
    for (const auto& mci : mi.integrals()) {

        // particles
        auto piPlus = std::static_pointer_cast<FinalStateParticle>(particle(*m.model(), is_named("pi+")));
        auto piMinus = std::static_pointer_cast<FinalStateParticle>(particle(*m.model(), is_named("pi-")));
        auto rho   = rho_rho     ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("rho0"))) : nullptr;
        auto omega = omega_omega ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("omega"))) : nullptr;
        auto sigma = (a_sigma_pi or sigma_pipi) ?  std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(500)"))) : nullptr;
        auto a_1_plus   = (a_rho_pi_S or a_rho_pi_D or a_sigma_pi) ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("a_1+"))) : nullptr;
        auto a_1_minus  = (a_rho_pi_S or a_rho_pi_D or a_sigma_pi) ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("a_1-"))) : nullptr;
        auto f_0   = f_0_pipi   ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0"))) : nullptr;
        auto f_2   = f_2_pipi   ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_2"))) : nullptr;
        auto f_0_1370 = sigma_f_0_1370 ? std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(1370)"))) : nullptr;

        auto decayTrees = mci.Integral.decayTrees();
        std::vector<yap::DecayTreeVector> groupedDecayTrees;


        // a_rho_pi_S
        if (a_rho_pi_S) {
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
        if (a_rho_pi_D) {
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
        if (a_sigma_pi) {
            auto blub = yap::filter(decayTrees, yap::to(a_1_plus));
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

        // a_rho_pi_S
        if (a_rho_pi_S) {
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
        if (a_rho_pi_D) {
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
        if (a_sigma_pi) {
            auto blub = yap::filter(decayTrees, yap::to(a_1_minus));
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

        // sigma_f_0_1370
        if (f_0_1370) {
            groupedDecayTrees.push_back(yap::filter(decayTrees, yap::to(f_0_1370)));

            for(auto dt : groupedDecayTrees.back())  {
                auto iter = std::find(decayTrees.begin(), decayTrees.end(), dt);
                if(iter != decayTrees.end())
                    decayTrees.erase(iter);
            }
        }

        if (flat_4pi)
            groupedDecayTrees.push_back(yap::filter(decayTrees, yap::to(piPlus, piMinus, piPlus, piMinus)));


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


        /*if (a_rho_pi_S and a_rho_pi_D and a_sigma_pi and rho_rho and f_0_pipi and f_2_pipi and sigma_pipi) {
            // scaling to match FOCUS
            LOG(INFO) << "\nnew scales to match FOCUS model";
            std::vector<double> focusFractions{43.3, 2.5, 8.3, 24.5, 2.4, 4.9, 8.2};
            std::vector<double> currentScales{1.0, scale_a_rho_pi_D, scale_a_sigma_pi, scale_rho_rho, scale_f_0_pipi, scale_f_2_pipi, scale_sigma_pipi};
            double fix = sqrt(focusFractions[0] / ff[0].value()); // a_rho_pi_S
            //assert(ff.size() == focusFractions.size() and ff.size() == currentScales.size());
            for (size_t i = 1; i < ff.size(); ++i) {
                LOG(INFO) << currentScales[i] * sqrt(focusFractions[i] / ff[i].value()) / fix;
            }
        }*/

    }
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
