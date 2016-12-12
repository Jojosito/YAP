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
#include <HelicityFormalism.h>
#include <MassRange.h>
#include <PHSP.h>
#include <RelativisticBreitWigner.h>
#include <ZemachFormalism.h>
#include <logging.h>
#include <make_unique.h>

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <TCanvas.h>
#include <TChain.h>
#include <TEventList.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <random>

int main()
{
    yap::plainLogs(el::Level::Info);

    std::string model_name = "D4PI_data";;

    // real data
    TChain* t = new TChain("t");
    t->Add("/nfs/hicran/scratch/user/jrauch/CopiedFromKEK/DataSkim_0?_analysis.root");
    t->AddFriend("t", "/nfs/hicran/scratch/user/jrauch/CopiedFromKEK/DataSkim_analysis_TMVA_weights.root");

    // MC data
    TChain* t_mcmc = new TChain("t");
    // stream 1
    t_mcmc->Add("/nfs/hicran/scratch/user/jrauch/CopiedFromKEK/FourPionsSkim_s01_??_analysis.root");
    t_mcmc->AddFriend("t", "/nfs/hicran/scratch/user/jrauch/CopiedFromKEK/");

    double BDT_cut = 0.2;

    // create model
    bat_fit* m = new bat_fit(d4pi_fit(model_name + "_fit"));
    double D_mass = 1.8648400; // D0

    // load fit data and partition it
    LOG(INFO) << "Load data";
    load_data_4pi(m->fitData(), *t, 200000, BDT_cut);
    m->fitPartitions() = yap::DataPartitionBlock::create(m->fitData(), 4);


    // load integration (MC) data
    load_data_4pi(m->integralData(), *t_mcmc, 200000, BDT_cut, true);
    m->integralPartitions() = yap::DataPartitionBlock::create(m->integralData(), 4);

    // open log file
    BCLog::OpenLog("output/" + m->GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m->SetPrecision(BCEngineMCMC::kMedium);
    m->SetNIterationsPreRunMax(1e6);
    m->SetNChains(4);
    // m->SetMinimumEfficiency(0.85);
    // m->SetMaximumEfficiency(0.99);

    m->SetNIterationsRun(static_cast<int>(1e5 / m->GetNChains()));

    // start timing:
    auto start = std::chrono::steady_clock::now();

    // run MCMC, marginalizing posterior
    m->MarginalizeAll(BCIntegrate::kMargMetropolis);

    // end timing
    auto end = std::chrono::steady_clock::now();

    m->FindMode(m->GetBestFitParameters());

    m->PrintSummary();
    m->PrintAllMarginalized("output/" + m->GetSafeName() + "_plots.pdf", 2, 2);

    // timing:
    auto diff = end - start;
    auto ms = std::chrono::duration<double, std::micro>(diff).count();
    auto nevents = (m->GetNIterationsPreRun() + m->GetNIterationsRun()) * m->GetNChains();
    BCLog::OutSummary(std::string("Seconds = ") + std::to_string(ms / 1.e6) + " for " + std::to_string(nevents) + " iterations, " + std::to_string(m->likelihoodCalls()) + " calls");
    BCLog::OutSummary(std::to_string(ms / nevents) + " microsec / iteration");
    BCLog::OutSummary(std::to_string(ms / m->likelihoodCalls()) + " microsec / call");

    // // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    delete m;

    return 0;
}
