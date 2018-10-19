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
#include <fstream>
#include <random>

#include <fenv.h>

int main()
{
    yap::plainLogs(el::Level::Global);

    const std::string model_name = "D4PI_data";;

    const double BDT_cut = 0.12; // 0.12 // negative for BG fit
    const double K0_cut = 3. * 0.00397333297611; // sigma from constrained masses

    // create model
    bat_fit_fractions m(d4pi_fit_fractions(model_name + "_fit_fractions",
            std::string(
"Summary :   0)  Parameter \"real(0 D0 --> a_1(1420)+, pi- L = 1 S = 1)\"    : -0.6243 +- 0.1098\n\
Summary :   1)  Parameter \"imag(0 D0 --> a_1(1420)+, pi- L = 1 S = 1)\"    : -0.9217 +- 0.1062\n\
Summary :   2)  Parameter \"real(1 D0 --> a_1(1640)+, pi- L = 1 S = 1)\"    : -1.604 +- 0.04952\n\
Summary :   3)  Parameter \"imag(1 D0 --> a_1(1640)+, pi- L = 1 S = 1)\"    : -2.932 +- 0.04551\n\
Summary :   4)  Parameter \"real(2 D0 --> f_0, f_0 L = 0 S = 0)\"           : -30.05 +- 0.5703\n\
Summary :   5)  Parameter \"imag(2 D0 --> f_0, f_0 L = 0 S = 0)\"           : 2.086 +- 0.5612\n\
Summary :   6)  Parameter \"real(3 D0 --> f_2, f_2 L = 0 S = 0)\"           : 9119 +- 1219\n\
Summary :   7)  Parameter \"imag(3 D0 --> f_2, f_2 L = 0 S = 0)\"           : 4.187e+04 +- 1051\n\
Summary :   8)  Parameter \"real(4 D0 --> pi(1300)+, pi- L = 0 S = 0)\"     : 0.9394 +- 0.01492\n\
Summary :   9)  Parameter \"imag(4 D0 --> pi(1300)+, pi- L = 0 S = 0)\"     : -0.1199 +- 0.01473\n\
Summary :  10)  Parameter \"real(5 D0 --> pi(1800)+, pi- L = 0 S = 0)\"     : -1.872 +- 0.08975\n\
Summary :  11)  Parameter \"imag(5 D0 --> pi(1800)+, pi- L = 0 S = 0)\"     : -1.042 +- 0.0893\n\
Summary :  12)  Parameter \"real(6 D0 --> pipiS, f_0 L = 0 S = 0)\"         : 6.341 +- 0.1313\n\
Summary :  13)  Parameter \"imag(6 D0 --> pipiS, f_0 L = 0 S = 0)\"         : -0.6367 +- 0.1302\n\
Summary :  14)  Parameter \"real(7 D0 --> pipiS, f_2 L = 2 S = 2)\"         : -8.045 +- 0.6819\n\
Summary :  15)  Parameter \"imag(7 D0 --> pipiS, f_2 L = 2 S = 2)\"         : 3.286 +- 0.8542\n\
Summary :  16)  Parameter \"real(8 D0 --> pipiS, pipiS L = 0 S = 0)\"       : -3.171 +- 0.0599\n\
Summary :  17)  Parameter \"imag(8 D0 --> pipiS, pipiS L = 0 S = 0)\"       : -1.276 +- 0.05688\n\
Summary :  18)  Parameter \"real(9 D0 --> rho0, rho0 L = 0 S = 0)\"         : -0.5887 +- 0.01974\n\
Summary :  19)  Parameter \"imag(9 D0 --> rho0, rho0 L = 0 S = 0)\"         : 0.6231 +- 0.01809\n\
Summary :  20)  Parameter \"real(10 D0 --> rho0, rho0 L = 1 S = 1)\"        : -0.01686 +- 0.02762\n\
Summary :  21)  Parameter \"imag(10 D0 --> rho0, rho0 L = 1 S = 1)\"        : 0.09199 +- 0.03549\n\
Summary :  22)  Parameter \"real(11 D0 --> rho0, rho0 L = 2 S = 2)\"        : -2.218 +- 0.0373\n\
Summary :  23)  Parameter \"imag(11 D0 --> rho0, rho0 L = 2 S = 2)\"        : -1.308 +- 0.0353\n\
Summary :  24)  Parameter \"real(12 a_1(1640)+ --> f_2, pi+ L = 1 S = 2)\"  : 0.861 +- 0.6925\n\
Summary :  25)  Parameter \"imag(12 a_1(1640)+ --> f_2, pi+ L = 1 S = 2)\"  : -0.0748 +- 0.7222\n\
Summary :  26)  Parameter \"real(13 a_1(1640)+ --> rho0, pi+ L = 2 S = 1)\" : 0.08458 +- 0.04543\n\
Summary :  27)  Parameter \"imag(13 a_1(1640)+ --> rho0, pi+ L = 2 S = 1)\" : 0.01641 +- 0.04266\n\
Summary :  28)  Parameter \"real(14 a_1+ --> f_2, pi+ L = 1 S = 2)\"        : 8.104 +- 2.33\n\
Summary :  29)  Parameter \"imag(14 a_1+ --> f_2, pi+ L = 1 S = 2)\"        : 0.6459 +- 2.426\n\
Summary :  30)  Parameter \"real(15 a_1+ --> pipiS, pi+ L = 1 S = 0)\"      : 0.2447 +- 0.03713\n\
Summary :  31)  Parameter \"imag(15 a_1+ --> pipiS, pi+ L = 1 S = 0)\"      : -2.343 +- 0.03136\n\
Summary :  32)  Parameter \"real(16 a_1+ --> rho0, pi+ L = 2 S = 1)\"       : 1.462 +- 0.07874\n\
Summary :  33)  Parameter \"imag(16 a_1+ --> rho0, pi+ L = 2 S = 1)\"       : 0.4838 +- 0.079\n\
Summary :  34)  Parameter \"real(17 pi(1300)+ --> rho0, pi+ L = 1 S = 1)\"  : 0.1879 +- 0.01835\n\
Summary :  35)  Parameter \"imag(17 pi(1300)+ --> rho0, pi+ L = 1 S = 1)\"  : 0.5084 +- 0.01664\n\
Summary :  36)  Parameter \"18 admixture_flat_4pi_0\"                       : 5.636 +- 0.1364"
            )));

    const unsigned nPreRun = 0;
    const unsigned nRun    = 10000;
    const unsigned nChains = 4;

    const unsigned nData = 300000; // 300000; max number of Data points we want
    const bool chargeBlind = false;
    const unsigned nThreads = 4;

    std::string dir("/home/ne53mad/CopiedFromKEK/");


    // real data, pre-filtered
    {
        TChain t("t");
        if (BDT_cut >= 0) {
            if (chargeBlind) {
                t.Add((dir + "DataSkim_analysis_chargeBlind_bdt_gt_0.025.root").c_str());
                t.AddFriend("t", (dir + "DataSkim_analysis_chargeBlind_bdt_gt_0.025_TMVA_weights.root").c_str());
            }
            else {
                t.Add((dir + "DataSkim_analysis_bdt_gt_0.025.root").c_str());
                t.AddFriend("t", (dir + "DataSkim_analysis_bdt_gt_0.025_TMVA_weights.root").c_str());
            }
        }
        else { // BG
            if (BDT_cut < -0.2) {
                t.Add((dir + "DataSkim_analysis_bdt_lt_-0.2.root").c_str());
                t.AddFriend("t", (dir + "DataSkim_analysis_bdt_lt_-0.2_TMVA_weights.root").c_str());
            }
            else if (BDT_cut < -0.1) {
                t.Add((dir + "DataSkim_analysis_bdt_lt_-0.1.root").c_str());
                t.AddFriend("t", (dir + "DataSkim_analysis_bdt_lt_-0.1_TMVA_weights.root").c_str());
            }
            else {
                t.Add("/nfs/hicran/scratch/user/jrauch/CopiedFromKEK/DataSkim_??_analysis.root");
                t.AddFriend("t", "/nfs/hicran/scratch/user/jrauch/CopiedFromKEK/DataSkim_analysis_TMVA_weights.root");
            }
        }
        LOG(INFO) << "Load data";
        load_data_4pi(m.fitData(), t, nData, BDT_cut, K0_cut, false);
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
            load_data_4pi(m.integralData(), t_mcmc, nData, BDT_cut, K0_cut, false); // phsp cut was already applied
        }
    }

    // fit fraction data
    { // root generated PHSP
        LOG(INFO) << "Load fitFractionData data";
        TChain t_phsp("t");
        t_phsp.Add("/home/ne53mad/YAPfork/output/D4pi_phsp_root.root");
        load_data_4pi(m.fitFractionData(), t_phsp, nData, 0, 0, false);
    }

    // partition data
    m.fitPartitions() = yap::DataPartitionBlock::create(m.fitData(), nThreads);
    m.integralPartitions() = yap::DataPartitionBlock::create(m.integralData(), nThreads);
    m.fitFractionPartitions() = yap::DataPartitionBlock::create(m.fitFractionData(), nThreads);

    // force calculate all integrals
    yap::ImportanceSampler::calculate(m.fitFractionIntegral(), m.fitFractionPartitions(), true);
    yap::ImportanceSampler::calculate(m.modelIntegral(), m.integralPartitions(), true);

    // open log file
    BCLog::OpenLog("output/" + m.GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNIterationsPreRunMax(nPreRun);
    m.SetNChains(nChains);
    m.SetMinimumEfficiency(0.15);
    m.SetMaximumEfficiency(0.9);


    m.SetNIterationsRun(static_cast<int>(nRun / m.GetNChains()));

    LOG(INFO) << "Initialization scheme: " << m.GetInitialPositionScheme();

    // start timing:
    const auto start = std::chrono::steady_clock::now();

    chi2(m);
    chi2_adaptiveBinning(m);
    chi2_adaptiveBinning(m, true);

    LOG(INFO) << "MCMCUserInitialize";
    m.MCMCUserInitialize();


    // set fit fraction observable ranges
    {
        m.setParameters(get_bat_parameter_values(m));
        m.CalculateObservables(get_bat_parameter_values(m));
        auto& observables = m.GetObservables();
        for (unsigned i = 0; i < observables.Size(); ++i)
        {
            auto& observable = observables.At(i);
            if (observable.GetName().find("fit_frac") != std::string::npos) {
                observable.SetLimits(0., 10.*observable.Value());
                LOG(INFO) << "set range for " << observable.GetName() << " to " << observable.GetUpperLimit();
            }
        }
    }



    // run MCMC, marginalizing posterior
    m.PrintSummary();
    m.setUseJacobian(true);
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    auto bestFitPars = m.GetBestFitParameters();
    std::ofstream out("output/MCMC_result.dat");
    out << forRootAsTxt(m, bestFitPars);
    out.close();


    LOG(INFO) << "amplitudes after fit:";
    printAmplitudes(*m.model());

    // end timing
    const auto end = std::chrono::steady_clock::now();

    LOG(INFO) << "Best fit parameters:";
    for (auto par : bestFitPars)
        LOG(INFO) << "\t" << par;


    // delete data
    m.fitData().clear();
    m.integralData().clear();

    d4pi_printFitFractions(m);
    m.fitFractionData().clear();

    LOG(INFO) << "Print Summary";
    m.PrintSummary();
    LOG(INFO) << "Generate plots";
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
