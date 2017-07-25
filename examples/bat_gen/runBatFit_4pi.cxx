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
    m.setModelSelection(0.1); // LASSO/BCM parameter

    bool adjustRangesFF = true; // adjust the real and imag ranges so that all waves have similar ff

    bool mcmc = true;
    unsigned nPreRun = 50000;
    unsigned nRun    = 20000;

    const unsigned nData = 10000; // max number of Data points we want
    const bool chargeBlind = false;
    const unsigned nThreads = 4;
    const double maxHours = 0.1;


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
        if (chargeBlind and BDT_cut >= 0) {
            t.Add((dir + "DataSkim_analysis_chargeBlind_bdt_gt_0.025.root").c_str());
            t.AddFriend("t", (dir + "DataSkim_analysis_chargeBlind_bdt_gt_0.025_TMVA_weights.root").c_str());
        }
        else if (BDT_cut >= 0) {
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
    { // fit fraction data
        LOG(INFO) << "Generate data to calculate fit fractions";
        generate_fit_fraction_data(m, nData, nThreads);
    }


    // MC generated data
    /*dir = "/home/ne53mad/YAPfork/buildRelease/output/";
    const double D0_mass = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/d4pi.pdl")["D0"].mass();
    {
        //TChain t("D4pi_rho_pi_S_flatBG_mcmc");
        //t.Add((dir + "D4pi_rho_pi_S_flatBG_mcmc.root").c_str());
        //TChain t("D4pi_rho_pi_S_rho2piBG_mcmc");
        //t.Add((dir + "D4pi_rho_pi_S_rho2piBG_mcmc.root").c_str());
        TChain t("D4pi_fitted_params_mcmc");
        t.Add((dir + "fit_06_2017_fixed_a1_shape/D4pi_fitted_params_mcmc.root").c_str());

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

    // force calculate all integrals
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
    m.SetNChains(nThreads); //
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
    }
    else {
        m.setUseJacobian(false);
        LOG(INFO) << "FindMode";


        //m.FindMode(m.getInitialPositions());
        //m.FindMode(m.getInitialPositions(true)); // start free amplitudes with 0 phase around 0

        //m.FindMode({-9.58814, -6.50694, -8.64714, -6.52594, -1.53999, 0.0525196, 0.703545, -1.35909, 2.1226, -0.314182, -96.9996, -89.4811, 0.0839771, -0.268775, 3.60997, -0.678504, -47.912, -39.3641, -0.202606, -0.26742, -1.12792, -2.09835, -9.52495, -8.51107, 17.4215, 15.7363, -0.00556024, -0.0047041, -0.80164, 0.478782, -1.25677, 1.52691, 0.0772175, -0.0794042, 0.269826, -0.258303, 0.00347969, 0.00144554, -2.80713, 0.457574, -0.699671, 0.708493, 0.0207715, -0.0121311, 0.0444898, 0.0996983, -0.0235575, -0.0221101, 3.94205, -8.6873, 0.0593776, 0.0236461, -0.0116456, -0.00305649, 0.155337, 0.477325, 0.279316});
        m.FindMode({0.0557778, 0.0633547, -4.38364, -10.452, -2.1563, -5.85761, 0.419586, -0.294097, -0.994249, -0.95886, -1.79373, 9.71762, -0.00434213, -0.000476509, -9.32922, -1.60055, -0.0526993, -0.415825, 0.307403, -0.345281, 0.0842681, 1.53982, -0.0800081, 0.045667, 0.00114439, -0.003268, -0.00758869, -0.0102014, 0.0223363, 0.0335351, -0.602239, -0.281768, 0.00931527, 0.0426345, 0.0100023, 0.000121869, 0.00217086, 1.18898, 0.490842, 0.0427989, 1.09312, 0.337672, 0.252606, -0.25507, 0.0844096, 48.6391, 11.334, -112.754, 6.24189, -110.21, 0.512392, -35.0274, 1.38128, -136.038, 9.88178, 100.458, 0.00436819, -173.737, 9.46552, -170.265, 0.419151, -97.2228, 0.462293, -48.3214, 1.54212, 86.8676, 0.0921237, 150.283, 0.00346258, -70.7008, 0.0127144, -126.645, 0.0402929, 56.3342, 0.664894, -154.927, 0.0436403, 77.6751});

        // keep on searching for longer
        unsigned iteration(0);
        auto elapsed = std::chrono::duration_cast<std::chrono::hours>(std::chrono::steady_clock::now() - start).count();
        LOG(INFO) << "time elapsed: " << elapsed << " hours";
        while (elapsed < maxHours) {

            LOG(INFO) << "random initial parameters of size " << m.getRandomInitialPositions().size();

            auto mode = m.FindMode(m.getRandomInitialPositions());

            auto bestFitPars = m.GetBestFitParameters();

            elapsed = std::chrono::duration_cast<std::chrono::hours>(std::chrono::steady_clock::now() - start).count();

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
    m.fitFractionData().clear();

    LOG(INFO) << "Print Summary";
    m.PrintSummary();
    LOG(INFO) << "Generate plots";
    m.PrintAllMarginalized("output/" + m.GetSafeName() + "_plots.pdf", 2, 2);
    if (mcmc) {
        try { m.PrintParameterPlot("output/" + m.GetSafeName() + "_parameterPlots.pdf", 2, 2); }
        catch (...) { LOG(ERROR) << "Failed PrintParameterPlot"; }
        try { m.PrintCorrelationMatrix("output/" + m.GetSafeName() + "_matrix.pdf"); }
        catch (...) { LOG(ERROR) << "Failed PrintCorrelationMatrix"; }
        /*try { m.PrintCorrelationPlot("output/" + m.GetSafeName() + "_correlation_observables.pdf"); }
        catch (...) { LOG(ERROR) << "Failed PrintCorrelationPlot for obsevables"; }
        try { m.PrintCorrelationPlot("output/" + m.GetSafeName() + "_correlation.pdf", false); }
        catch (...) { LOG(ERROR) << "Failed PrintCorrelationPlot"; }*/
    }


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
                << "  \t (real, imag) = (" << real(fa->value()) << ", " << imag(fa->value()) << ")"
                << "  \t std::polar(" <<abs(fa->value()) << ", yap::rad(" << deg(arg(fa->value())) << "));";

        const double D0_mass = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/d4pi.pdl")["D0"].mass();
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

        m_gen.SetNIterationsRun(static_cast<int>(5e6 / m_gen.GetNChains()));

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
