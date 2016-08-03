// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "fit_fitFraction.h"
#include "models/d4pi.h"

#include <DecayTree.h>
#include <Filters.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MassRange.h>
#include <Model.h>
#include <PHSP.h>
#include <ZemachFormalism.h>

#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <random>

const double quad(std::vector<double> S)
{ return sqrt(std::accumulate(S.begin(), S.end(), 0., [](double a, double s) {return a + s * s;})); }

template <typename ... Types>
constexpr double quad(double s0, Types ... additional)
{ return quad({s0, additional...}); }

int main()
{
    yap::plainLogs(el::Level::Info);

    auto m = d4pi_fit_fitFraction();

    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(m.axes(), m.isp(), m.model()->finalStateParticles()));

    // generate integration data
    std::mt19937 g(0);
    std::generate_n(std::back_inserter(m.integralData()), 10000,
                    std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()),
                              m.isp()->mass()->value(), m.axes(), m2r, g,
                              std::numeric_limits<unsigned>::max()));
    LOG(INFO) << "Created " << m.integralData().size() << " data points (" << (m.integralData().bytes() * 1.e-6) << " MB)";
    // partition integration data
    m.integralPartitions() = yap::DataPartitionBlock::create(m.integralData(), 2);

    // open log file
    BCLog::OpenLog("output/" + m.GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNIterationsPreRunMax(1e6);
    m.SetNChains(4);
    // m.SetMinimumEfficiency(0.85);
    // m.SetMaximumEfficiency(0.99);

    m.SetNIterationsRun(static_cast<int>(2e4 / m.GetNChains()));

    m.WriteMarkovChain("output/" + m.GetSafeName() + "_mcmc.root", "RECREATE");

    // start timing:
    auto start = std::chrono::steady_clock::now();

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // end timing
    auto end = std::chrono::steady_clock::now();

    for (size_t i = 0; i < D->decayTrees().at(0).size(); ++i)
        std::cout << i << " = " << to_string(*D->decayTrees().at(0)[i]) << std::endl;

    m.PrintSummary();
    m.PrintAllMarginalized("output/" + m.GetSafeName() + "_plots.pdf", 2, 2);
    m.SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDetailedPosterior);
    m.PrintKnowledgeUpdatePlots("output/" + m.GetSafeName() + "_update.pdf", 2, 2, true);

    // timing:
    auto diff = end - start;
    auto ms = std::chrono::duration<double, std::micro>(diff).count();
    auto nevents = (m.GetNIterationsPreRun() + m.GetNIterationsRun()) * m.GetNChains();
    BCLog::OutSummary(std::string("Seconds = ") + std::to_string(ms / 1.e6) + " for " + std::to_string(nevents) + " iterations, " + std::to_string(m.likelihoodCalls()) + " calls");
    BCLog::OutSummary(std::to_string(ms / nevents) + " microsec / iteration");
    BCLog::OutSummary(std::to_string(ms / m.likelihoodCalls()) + " microsec / call");

    // // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
