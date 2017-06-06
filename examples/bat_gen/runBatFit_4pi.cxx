// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_fit.h"
#include "bat_gen.h"
#include "models/d4pi_fit.h"
#include "tools.h"

#include <Exceptions.h>
#include <logging.h>
#include <make_unique.h>

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <TChain.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <random>

#include <fenv.h>

int main()
{
    yap::plainLogs(el::Level::Global);
    //feenableexcept (FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
    //feenableexcept (FE_INVALID|FE_OVERFLOW);

    std::string model_name = "D4PI_data";;

    const double BDT_cut = 0.15; // 0.15 // negative for BG fit
    const double K0_cut = 3. * 0.00397333297611; // sigma from constrained masses

    // create model
    bat_fit m(d4pi_fit(model_name + "_fit"));

    const unsigned nData = 200000; // max number of Data points we want
    const unsigned nThreads = 4;


    std::string dir("/home/ne53mad/CopiedFromKEK/");

    // real data
    //std::string dir = "/nfs/hicran/scratch/user/jrauch/CopiedFromKEK/";
    /*{
        TChain t("t");
        t.Add((dir + "DataSkim_0?_analysis.root").c_str());
        t.AddFriend("t", (dir + "DataSkim_analysis_TMVA_weights.root").c_str());
        LOG(INFO) << "Load data";
        load_data_4pi(m.fitData(), t, nData, BDT_cut, K0_cut, false);
        m.fitPartitions() = yap::DataPartitionBlock::create(m.fitData(), 4);
    }
    // MC data, load from several streams if needed
    for (int i = 0; i < 6; ++i) {
        auto nLoaded = m.integralData().size();
        if (nLoaded >= nData)
            break;
        TChain t_mcmc("t");
        t_mcmc.Add(dir + TString::Format("FourPionsSkim_s%d_analysis.root", i));
        t_mcmc.AddFriend("t", dir + TString::Format("FourPionsSkim_analysis_s%d_TMVA_weights.root", i));
        LOG(INFO) << "Load integration (MC) data";
        load_data_4pi(m.integralData(), t_mcmc, nData-nLoaded, BDT_cut, K0_cut, true);
    }*/

    // real data, pre-filtered
    {
        TChain t("t");
        if (BDT_cut >= 0) {
            t.Add((dir + "DataSkim_analysis_bdt_gt_0.025.root").c_str());
            t.AddFriend("t", (dir + "DataSkim_analysis_bdt_gt_0.025_TMVA_weights.root").c_str());
        }
        else { // BG
            if (BDT_cut > -0.2) {
                t.Add((dir + "DataSkim_analysis_bdt_lt_-0.1.root").c_str());
                t.AddFriend("t", (dir + "DataSkim_analysis_bdt_lt_-0.1_TMVA_weights.root").c_str());
            }
            else {
                t.Add((dir + "DataSkim_analysis_bdt_lt_-0.2.root").c_str());
                t.AddFriend("t", (dir + "DataSkim_analysis_bdt_lt_-0.2_TMVA_weights.root").c_str());
            }
        }
        LOG(INFO) << "Load data";
        load_data_4pi(m.fitData(), t, nData, BDT_cut, K0_cut, false);
        m.fitPartitions() = yap::DataPartitionBlock::create(m.fitData(), nThreads);
    }
    { // MC data, pre-filtered
        LOG(INFO) << "Load integration (MC) data";
        dir = "/home/ne53mad/CopiedFromKEK/";
        TChain t_mcmc("t");
        if (BDT_cut > 0.) {
            t_mcmc.Add((dir + "FourPionsSkim_analysis_phsp_bdt_gt_0.025.root").c_str());
            t_mcmc.AddFriend("t", (dir + "FourPionsSkim_analysis_phsp_bdt_gt_0.025_TMVA_weights.root").c_str());
            load_data_4pi(m.integralData(), t_mcmc, nData, BDT_cut, K0_cut, false); // phsp cut was already applied
        }
        else {
            t_mcmc.Add((dir + "FourPionsSkim_analysis_phsp_bdt_lt_0.0.root").c_str());
            t_mcmc.AddFriend("t", (dir + "FourPionsSkim_analysis_phsp_bdt_lt_0.0_TMVA_weights.root").c_str());
            load_data_4pi(m.integralData(), t_mcmc, nData, -0.1, K0_cut, false); // phsp cut was already applied
        }
    }


    // MC generated data
    /*dir = "/home/ne53mad/YAPfork/buildRelease/output/";
    const double D0_mass = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl")["D0"].mass();
    {
        //TChain t("D4pi_rho_pi_S_flatBG_mcmc");
        //t.Add((dir + "D4pi_rho_pi_S_flatBG_mcmc.root").c_str());
        //TChain t("D4pi_rho_pi_S_rho2piBG_mcmc");
        //t.Add((dir + "D4pi_rho_pi_S_rho2piBG_mcmc.root").c_str());
        TChain t("D4pi_flat4pi_rho2piBG_mcmc");
        t.Add((dir + "D4pi_flat4pi_rho2piBG_mcmc.root").c_str());

        LOG(INFO) << "Load data";
        load_data(m.fitData(), *m.model(), m.axes(), D0_mass, t, nData, 10);
        m.fitPartitions() = yap::DataPartitionBlock::create(m.fitData(), nThreads);
    }
    { // MC data, pre-filtered
        TChain t_mcmc("D4pi_phsp_mcmc");
        t_mcmc.Add((dir + "D4pi_phsp_mcmc.root").c_str());
        LOG(INFO) << "Load integration (MC) data";
        load_data(m.integralData(), *m.model(), m.axes(), D0_mass, t_mcmc, nData, 25);
    }*/

    // partition integral data
    m.integralPartitions() = yap::DataPartitionBlock::create(m.integralData(), nThreads);


    // open log file
    BCLog::OpenLog("output/" + m.GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNIterationsPreRunMax(1e5);
    m.SetNChains(nThreads); //
    m.SetMinimumEfficiency(0.15);
    m.SetMaximumEfficiency(0.35);


    m.SetNIterationsRun(static_cast<int>(1e5 / m.GetNChains()));

    LOG(INFO) << "Initialization scheme: " << m.GetInitialPositionScheme();

    // start timing:
    const auto start = std::chrono::steady_clock::now();

    chi2(m);
    chi2_adaptiveBinning(m);

    LOG(INFO) << "MCMCUserInitialize";
    m.MCMCUserInitialize();

    // run MCMC, marginalizing posterior
    bool mcmc = false;
    if (mcmc) {
        m.setUseJacobian(true);
        m.MarginalizeAll(BCIntegrate::kMargMetropolis);
    }
    else {
        m.setUseJacobian(false);
        LOG(INFO) << "FindMode";
        m.FindMode(m.getInitialPositions());
        //for (unsigned i = 0; i<10; ++i)
        //    m.FindMode(m.getRandomInitialPositions());
    }


    // end timing
    const auto end = std::chrono::steady_clock::now();

    auto bestFitPars = m.GetBestFitParameters();
    LOG(INFO) << "Best fit parameters:";
    for (auto par : bestFitPars)
        LOG(INFO) << "\t" << par;

    if (mcmc) {
        LOG(INFO) << "\nFind global mode";
        m.setUseJacobian(false);
        auto globalMode = m.FindMode(bestFitPars);

        double llGlobMode = std::numeric_limits<double>::min();
        if (globalMode.empty())
            LOG(ERROR) << "global mode is empty";
        else
            llGlobMode = m.LogLikelihood(globalMode);

        LOG(INFO) << "\nLogLikelihood of global mode = " << llGlobMode;
    }

    chi2(m);
    chi2_adaptiveBinning(m);

    // delete data
    m.fitData().clear();
    m.integralData().clear();

    d4pi_printFitFractions(m);

    LOG(INFO) << "Print Summary";
    m.PrintSummary();
    LOG(INFO) << "Generate plots";
    m.PrintAllMarginalized("output/" + m.GetSafeName() + "_plots.pdf", 2, 2);


    //LOG(INFO) << "LogLikelihood of global mode = " << llGlobMode;

    // timing:
    const auto diff = end - start;
    const auto ms = std::chrono::duration<double, std::micro>(diff).count();
    const auto nevents = (m.GetNIterationsPreRun() + m.GetNIterationsRun()) * m.GetNChains();
    BCLog::OutSummary(std::string("Seconds = ") + std::to_string(ms / 1.e6) + " for " + std::to_string(nevents) + " iterations, " + std::to_string(m.likelihoodCalls()) + " calls");
    BCLog::OutSummary(std::to_string(ms / nevents) + " microsec / iteration");
    BCLog::OutSummary(std::to_string(ms / m.likelihoodCalls()) + " microsec / call");

    // // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();


    // generate with fitted amplitudes
    if (true) {
        LOG(INFO) << "Generate data with fitted amplitudes";

        auto model = std::unique_ptr<yap::Model>(nullptr);
        model.swap(m.model());

        LOG(INFO) << "\nAmplitudes used for YAP generation: ";
        for (const auto& fa : free_amplitudes(*model))
            LOG(INFO) << yap::to_string(*fa) << "  \t (mag, phase) = (" << abs(fa->value()) << ", " << deg(arg(fa->value())) << "°)"
                << "  \t (real, imag) = (" << real(fa->value()) << ", " << imag(fa->value()) << ")";

        const double D0_mass = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl")["D0"].mass();
        bat_gen m_gen("D4pi_fitted_params", std::move(model), D0_mass);

        // open log file
        BCLog::OpenLog("output/" + m_gen.GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

        // set precision
        m_gen.SetPrecision(BCEngineMCMC::kMedium);
        m_gen.SetNChains(4);
        m_gen.SetNIterationsPreRunMax(1e6);
        m_gen.SetMinimumEfficiency(0.7);
        m_gen.SetMaximumEfficiency(0.9);
        m_gen.SetInitialPositionAttemptLimit(1e5);

        m_gen.SetNIterationsRun(static_cast<int>(1e8 / m_gen.GetNChains()));

        m_gen.WriteMarkovChain("output/" + m_gen.GetSafeName() + "_mcmc.root", "RECREATE", true, false);

        // start timing:
        auto start = std::chrono::steady_clock::now();

        // run MCMC, marginalizing posterior
        m_gen.MarginalizeAll(BCIntegrate::kMargMetropolis);

        // end timing
        auto end = std::chrono::steady_clock::now();

        // timing:
        auto diff = end - start;
        auto ms = std::chrono::duration<double, std::micro>(diff).count();
        auto nevents = (m_gen.GetNIterationsPreRun() + m_gen.GetNIterationsRun()) * m_gen.GetNChains();
        BCLog::OutSummary(std::string("Seconds = ") + std::to_string(ms / 1.e6) + " for " + std::to_string(nevents) + " iterations, " + std::to_string(m_gen.likelihoodCalls()) + " calls");
        BCLog::OutSummary(std::to_string(ms / nevents) + " microsec / iteration");
        BCLog::OutSummary(std::to_string(ms / m_gen.likelihoodCalls()) + " microsec / call");

        // close log file
        BCLog::OutSummary("Exiting");
        BCLog::CloseLog();
    }

    return 0;
}
