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
    bat_fit m(d4pi_fit(model_name + "_fit"));
    const bool adjustRangesFF = true; // adjust the real and imag ranges so that all waves have similar ff

    //bat_fit m(d4pi_fit(model_name + "_fit", {4.03148, -32.5435, -4.51938, -1.59346, 136.728, -73.0129, -699953, -125468, -0.444905, -0.93394, 10.5643, 3.10409, -30.4939, 15.9983, 52.5308, -37.9762, 5.14999, 2.77111, 0.00602175, 0.756096, 236.493, 17.6788, 111.186, -55.8603, -2.25629, 1.67204, 3.4199, -1.36056, -1.87353, 1.34619, -17.2051, -33.5261, 1.16347, -10.496, -254.044, -77.1403, 8.07964, 3.70565, -26.6026, -21.6966, -5.52496, 1.28902, 1.02664e-05}));
    //const bool adjustRangesFF = false; // adjust the real and imag ranges so that all waves have similar ff

    //m.setModelSelection(0.75); // LASSO/BCM parameter


    const bool minuit = true;
    m.SetRelativePrecision(0.1); // default is 1%, to speed things up a bit

    const bool mcmc = false;
    const unsigned nPreRun = 1500;
    const unsigned nRun    = 100;
    const unsigned nChains = 128;

    const unsigned nData = 300000; // 300000; max number of Data points we want
    const bool chargeBlind = false;
    const unsigned nThreads = 4;
    const double maxHours = 18;

    const bool generate_mc = true;


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

    if (adjustRangesFF) {
        for (unsigned i = 0; i < 10; ++i)
            d4pi_normalizeFitFractions(m);
        m.setRealImagRanges();
    }


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


    // run MCMC, marginalizing posterior
    if (mcmc) {
        m.PrintSummary();
        m.setUseJacobian(true);
        m.MarginalizeAll(BCIntegrate::kMargMetropolis);

        auto bestFitPars = m.GetBestFitParameters();
        std::ofstream out("output/MCMC_result.dat");
        out << forRootAsTxt(m, bestFitPars);
        out.close();
    }
    if (minuit) {
        m.setUseJacobian(false);

        // keep on searching for longer
        unsigned iteration(0);
        auto elapsed = std::chrono::duration_cast<std::chrono::hours>(std::chrono::steady_clock::now() - start).count();
        LOG(INFO) << "time elapsed: " << elapsed << " hours";

        auto bestFitPars = m.GetBestFitParameters();

        while (elapsed <= maxHours) {

            // get random start parameters
            auto initPars = (iteration == 0) ? m.getInitialPositions() : m.getRandomInitialPositions();

            // set ranges very high
            if (adjustRangesFF) {
                if (iteration > 0)
                    m.setParameters(bestFitPars);
                m.setRealImagRanges(1e3);
                m.setAdmixtureRanges(1e3);
            }

            m.PrintSummary();

            std::ofstream out("output/Minuit_init_" + std::to_string(iteration) + ".dat");
            out << forRootAsTxt(m, initPars);
            out.close();

            // DO THE FIT
            auto mode = m.FindMode(initPars);

            std::ofstream out2("output/Minuit_result_" + std::to_string(iteration) + ".dat");
            out2 << forRootAsTxt(m, mode);
            out2.close();

            // set parameter ranges for next random parameter choice
            bestFitPars = m.GetBestFitParameters();

            if (adjustRangesFF) {
                m.setParameters(bestFitPars);
                m.setRealImagRanges(2.);
                m.setAdmixtureRanges(2.);
            }

            // get elapsed time
            elapsed = 1 + std::chrono::duration_cast<std::chrono::hours>(std::chrono::steady_clock::now() - start).count();

            LOG(INFO) << "Best fit parameters after iteration " << iteration++ << "; time elasped: " << elapsed << "h :";
            for (auto par : mode)
                LOG(INFO) << "\t" << par;
        }
    }

    LOG(INFO) << "amplitudes after fit:";
    printAmplitudes(*m.model());

    // end timing
    const auto end = std::chrono::steady_clock::now();

    auto bestFitPars = m.GetBestFitParameters();
    LOG(INFO) << "Best fit parameters:";
    for (auto par : bestFitPars)
        LOG(INFO) << "\t" << par;

    //if (mcmc) {
        LOG(INFO) << "\nFind global mode";
        m.setUseJacobian(false);
        m.SetRelativePrecision(0.1);
        auto globalMode = m.FindMode(bestFitPars);

        double llGlobMode = std::numeric_limits<double>::min();
        if (globalMode.empty())
            LOG(ERROR) << "global mode is empty";
        else
            llGlobMode = m.LogLikelihood(globalMode);

        LOG(INFO) << "\nLogLikelihood of global mode = " << llGlobMode;
    //}


    LOG(INFO) << "amplitudes after finding global mode:";
    printAmplitudes(*m.model());

    chi2(m);
    chi2_adaptiveBinning(m);
    chi2_adaptiveBinning(m, true);

    if (generate_mc) {
        LOG(INFO) << "Generate weight file with fitted amplitudes";
        write_weighted_integral_data(m, "output/fitResults.root");
    }

    // delete data
    m.fitData().clear();
    m.integralData().clear();

    d4pi_printFitFractions(m);
    m.fitFractionData().clear();

    LOG(INFO) << "Print Summary";
    m.PrintSummary();
    LOG(INFO) << "Generate plots";
    m.PrintAllMarginalized("output/" + m.GetSafeName() + "_plots.pdf", 2, 2);
    if (mcmc and not minuit) {
        try { m.PrintParameterPlot("output/" + m.GetSafeName() + "_parameterPlots.pdf", 2, 2); }
        catch (...) { LOG(ERROR) << "Failed PrintParameterPlot"; }
        try { m.PrintCorrelationMatrix("output/" + m.GetSafeName() + "_matrix.pdf"); }
        catch (...) { LOG(ERROR) << "Failed PrintCorrelationMatrix"; }
        /*try { m.PrintCorrelationPlot("output/" + m.GetSafeName() + "_correlation_observables.pdf"); }
        catch (...) { LOG(ERROR) << "Failed PrintCorrelationPlot for obsevables"; }
        try { m.PrintCorrelationPlot("output/" + m.GetSafeName() + "_correlation.pdf", false); }
        catch (...) { LOG(ERROR) << "Failed PrintCorrelationPlot"; }*/
    }

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
