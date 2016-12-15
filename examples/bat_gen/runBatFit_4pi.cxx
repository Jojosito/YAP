// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_fit.h"
#include "models/d4pi.h"
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

int main()
{
    yap::plainLogs(el::Level::Info);
    yap::plainLogs(el::Level::Debug);

    std::string model_name = "D4PI_data";;

    const double BDT_cut = 0.2;
    const double K0_cut = 3. * 0.00397333297611; // sigma from constrained masses

    // create model
    bat_fit m(d4pi_fit(model_name + "_fit"));

    // real data
    {
        TChain t("t");
        t.Add("/nfs/hicran/scratch/user/jrauch/CopiedFromKEK/DataSkim_0?_analysis.root");
        t.AddFriend("t", "/nfs/hicran/scratch/user/jrauch/CopiedFromKEK/DataSkim_analysis_TMVA_weights.root");
        LOG(INFO) << "Load data";
        load_data_4pi(m.fitData(), t, 300, BDT_cut, K0_cut, false);
        m.fitPartitions() = yap::DataPartitionBlock::create(m.fitData(), 4);
    }

    // MC data
    {
        TChain t_mcmc("t");
        t_mcmc.Add("/nfs/hicran/scratch/user/jrauch/CopiedFromKEK/FourPionsSkim_s0_analysis.root");
        t_mcmc.AddFriend("t", "/nfs/hicran/scratch/user/jrauch/CopiedFromKEK/FourPionsSkim_analysis_s0_TMVA_weights.root");
        LOG(INFO) << "Load integration (MC) data";
        load_data_4pi(m.integralData(), t_mcmc, 300, BDT_cut, K0_cut, true);
        m.integralPartitions() = yap::DataPartitionBlock::create(m.integralData(), 4);
    }

    // open log file
    BCLog::OpenLog("output/" + m.GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNIterationsPreRunMax(1e5);
    m.SetNChains(4);
    // m.SetMinimumEfficiency(0.85);
    // m.SetMaximumEfficiency(0.99);

    m.SetNIterationsRun(static_cast<int>(1e5 / m.GetNChains()));

    // start timing:
    const auto start = std::chrono::steady_clock::now();

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // end timing
    const auto end = std::chrono::steady_clock::now();

    m.FindMode(m.GetBestFitParameters());

    m.PrintSummary();
    m.PrintAllMarginalized("output/" + m.GetSafeName() + "_plots.pdf", 2, 2);

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

    return 0;
}
