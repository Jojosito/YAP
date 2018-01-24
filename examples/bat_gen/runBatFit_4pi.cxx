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
    //bat_fit m(d4pi_fit(model_name + "_fit", {   -0.982423, -1.11708, -37.4235, -27.7717, -3.08474, -0.609421, 5.86998, 10.4598, 0.832805, 1.53719, 4.3105, 2.70736, -4.50548, 5.54519, 30.3196, 200.503, 124190, 45460.8, 4.13664, -1.0761, -4.08, -3.8142, 2.32384, 4.50658, 61.379, 1.22991, 0.424383, 0.06, 1.48762, -131.33, 46.6024, -143.421, 3.14436, -168.825, 11.9943, 60.6991, 1.74829, 61.5524, 5.09021, 32.1323, 7.14482, 129.094, 202.782, 81.401, 132249, 20.1056, 4.27431, -14.5816, 5.58521, -136.928, 5.07046, 62.7219}));
    //bat_fit m(d4pi_fit(model_name + "_fit", {-0.154518, 0.227784, 0.177584, -0.392354, -20.707, 28.1154, 7.63314, -2.7599, 0.115322, 0.00264234, -2.82195, 2.16605, 82.8522, -19.4367, -1.0357, 0.421327, -23.0284, 3.68123, -0.631946, 0.476604, -0.475645, 0.227638, 1.02482, -0.668613, -404.364, 171.518, -6.79664, 2.55212, -6.78063, -12.8324, 1.41051, 6.41238, 114.503, 57.3835, -15.6803, -74.0754, -26.5365, -3.97416, 103.276, 8.50773, -750.726, -5590.62, 515.759, -105.61, 32.8272, 10.4841, 41302, 23896.5, -0.551647, 0.00909297, 1.65934, -0.810574, -0.474255, -0.664477, -1.77448, -0.158776, -3.92576, 1.96692, 0.0899473, -0.117294, 4.20759, 0.811348, 6.63504, 2.383, -4.99893, 2.51687, -0.210127, 0.248135, -4.66397, -0.950354, 1.07129, -0.0454435, -21.448, -33.0476, -1.14807, 2.97919, -0.24193, -0.0889014, 99.5525, 23.527, -0.865321, -1.21978, -2.0214, -0.45496, -0.0743514, -0.164393, -0.0360475, -0.09143, -46.8871, -58.8546, 0.0153325}));
    m.setModelSelection(0.75); // LASSO/BCM parameter

    const bool adjustRangesFF = true; // adjust the real and imag ranges so that all waves have similar ff

    const bool minuit = true;
    m.SetRelativePrecision(0.1); // default is 1%, to speed things up a bit

    const bool mcmc = false;
    const unsigned nPreRun = 1500;
    const unsigned nRun    = 100;
    const unsigned nChains = 128;

    const unsigned nData = 300000; // 300000; max number of Data points we want
    const bool chargeBlind = false;
    const unsigned nThreads = 4;
    const double maxHours = 37;

    const bool generate_mc = true;


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

    /*{ // root generated PHSP
        LOG(INFO) << "Load integral data";
        TChain t_phsp("t");
        t_phsp.Add("/home/ne53mad/YAPfork/output/D4pi_phsp_root.root");
        load_data_4pi(m.integralData(), t_phsp, nData, 0, K0_cut, false);
    }*/


    // MC generated data
    /*dir = "/home/ne53mad/YAPfork/buildRelease/output/";
    {
        TChain t("D4pi_mcmc");
        t.Add((dir + "D4pi_mcmc.root").c_str());

        LOG(INFO) << "Load data";
        load_data(m.fitData(), *m.model(), m.axes(), D0_mass, t, nData, -1);
        m.fitPartitions() = yap::DataPartitionBlock::create(m.fitData(), nThreads);
    }
    const double D0_mass = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/d4pi.pdl")["D0"].mass();
    { // MC data, pre-filtered
        TChain t_mcmc("D4pi_phsp_mcmc");
        t_mcmc.Add("/home/ne53mad/YAPfork/buildRelease/output/D4pi_phsp_mcmc.root");
        LOG(INFO) << "Load integration (MC) data";
        load_data(m.integralData(), *m.model(), m.axes(), D0_mass, t_mcmc, nData, -1);
    }*/

    // fit fraction data

    { // root generated PHSP
        LOG(INFO) << "Load fitFractionData data";
        TChain t_phsp("t");
        t_phsp.Add("/home/ne53mad/YAPfork/output/D4pi_phsp_root.root");
        load_data_4pi(m.fitFractionData(), t_phsp, nData, 0, 0, false);
    }

    /*{ // YAP generated; seems to have some problems!
        LOG(INFO) << "Generate data to calculate fit fractions";
        generate_fit_fraction_data(m, nData, nThreads);
    }*/

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

    chi2(m);
    chi2_adaptiveBinning(m);

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
