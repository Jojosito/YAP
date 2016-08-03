// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "fit_fitFraction.h"
#include "models/d4pi.h"

#include <DecayTree.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MassRange.h>
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

    //
    // CREATE yap::Model
    //
    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    auto M = std::make_unique<yap::Model>(std::make_unique<yap::HelicityFormalism>());

    M->setFinalState(piPlus, piMinus, piPlus, piMinus);

    // use common radial size for all resonances
    double radialSize = 1.2; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D0"), radialSize);

    //
    // resonant particles
    //

    // rho
    auto rho = F.resonance(F.pdgCode("rho0"), radialSize, std::make_shared<yap::RelativisticBreitWigner>());
    rho->addChannel(piPlus, piMinus);

    // omega
    //auto omega = F.resonance(F.pdgCode("omega"), radialSize, std::make_shared<yap::RelativisticBreitWigner>());
    //omega->addChannel({piPlus, piMinus});

    // sigma / f_0(500)
    auto sigma = F.resonance(F.pdgCode("f_0(500)"), radialSize, std::make_shared<yap::RelativisticBreitWigner>());
    sigma->addChannel(piPlus, piMinus);

    // a_1
    auto a_1 = F.resonance(F.pdgCode("a_1+"), radialSize, std::make_shared<yap::BreitWigner>());
    a_1->addChannel(rho,   piPlus);
    a_1->addChannel(sigma, piPlus);

    for (auto& freeAmp : a_1->freeAmplitudes(rho, piPlus)) {
        LOG(INFO) << to_string(*freeAmp);
        if (freeAmp->spinAmplitude()->L() == 0) {
            *freeAmp = 1.; // S-wave, 1
            freeAmp->variableStatus() = yap::VariableStatus::fixed; // S-wave, fixed
        }
        else if (freeAmp->spinAmplitude()->L() == 1)
            *freeAmp = 0.; // P-wave, 0
        else if (freeAmp->spinAmplitude()->L() == 2)
            *freeAmp = std::polar(0.241, yap::rad(82.)); // D-wave
        else
            LOG(ERROR) << "wrong spin";
    }

    *a_1->freeAmplitudes(sigma, piPlus)[0] = std::polar(0.439, yap::rad(193.));

    // f_0(980) (as Flatte)
    auto piZero = F.fsp(111);
    auto Kshort = F.fsp(310);
    auto f_0_980_flatte = std::make_shared<yap::Flatte>();
    f_0_980_flatte->addChannel(0.20, piZero->mass()->value());
    f_0_980_flatte->addChannel(0.50, Kshort->mass()->value());
    auto f_0_980 = F.resonance(F.pdgCode("f_0"), radialSize, f_0_980_flatte);
    f_0_980->addChannel(piPlus, piMinus);

    // f_2(1270)
    auto f_2 = F.resonance(F.pdgCode("f_2"), radialSize, std::make_shared<yap::RelativisticBreitWigner>());
    f_2->addChannel(piPlus, piMinus);

    // pi+ pi- flat
    auto pipiFlat = F.decayingParticle(F.pdgCode("f_0"), radialSize); // just need any spin0 particle
    pipiFlat->addChannel(piPlus, piMinus);

    //
    // D0 channels
    //

    D->addChannel(rho, rho);
    D->addChannel(a_1, piMinus);
    D->addChannel(f_0_980, piPlus, piMinus);
    D->addChannel(f_2, pipiFlat);
    D->addChannel(sigma, piPlus, piMinus);

    // R pi pi
    *D->freeAmplitudes(f_0_980, piPlus, piMinus)[0] = std::polar(0.233, yap::rad(261.));
    *D->freeAmplitudes(f_2,     pipiFlat       )[0] = std::polar(0.338, yap::rad(317.));
    *D->freeAmplitudes(sigma,   piPlus, piMinus)[0] = std::polar(0.432, yap::rad(254.));

    // rho rho
    // transform into angular momentum basis
    yap::basis::canonical<double> c(yap::basis::transversity<double>(
            std::polar(0.624, yap::rad(357.)),    // A_longitudinal
            std::polar(0.157, yap::rad(120.)),    // A_parallel
            std::polar(0.384, yap::rad(163.)) )); // A_perpendicular

    for (auto& freeAmp : D->freeAmplitudes(rho, rho)) {
        LOG(INFO) << to_string(*freeAmp);
        *freeAmp = c[freeAmp->spinAmplitude()->L()];
    }

    M->addInitialStateParticle(D);



    // create bat_fit object
    fit_fitFraction m("D3PI_frac_fit", std::move(M));

    m.GetParameter("N_1").Fix(1);

    // set fit fractions to fit
    m.setFitFraction(D->decayTrees(a_1,      piMinus)[0], 43.3e-2,   quad(2.5e-2, 1.9e-2));
    m.setFitFraction(D->decayTrees(a_1,      piMinus)[1], 2.5e-2,    quad(0.5e-2, 0.4e-2));
    m.setFitFraction(D->decayTrees(a_1,      piMinus)[2], 8.3e-2,    quad(0.7e-2, 0.6e-2));
    m.setFitFraction(D->decayTrees(rho,      rho)[0],    1.1e-2,    quad(0.3e-2, 0.3e-2));
    m.setFitFraction(D->decayTrees(rho,      rho)[1],    6.4e-2,    quad(0.6e-2, 0.5e-2));
    m.setFitFraction(D->decayTrees(rho,      rho)[2],    16.88e-2,  quad(1.0e-2, 0.8e-2));

    m.setFitFraction(D->decayTrees(f_0_980, piPlus, piMinus)[0], 2.4e-2,  quad(2.4e-2, 0.4e-2));
    m.setFitFraction(D->decayTrees(f_2,     pipiFlat)[0], 4.9e-2,  quad(4.9e-2, 0.5e-2));
    m.setFitFraction(D->decayTrees(sigma,   piPlus, piMinus)[0], 8.2e-2,  quad(8.2e-2, 0.7e-2));

    // set free amplitude parameters of fit
    m.fix(D->freeAmplitudes(a_1,      piMinus)[0], 1., 0.);
    m.setPrior(D->freeAmplitudes(a_1,      piMinus)[1], new BCGaussianPrior(0.241, quad(0.033, 0.024)), new BCGaussianPrior(  82., quad(5., 4.)));
    m.setPrior(D->freeAmplitudes(f_0_980,  piPlus)[0], new BCGaussianPrior(0.439, quad(0.026, 0.021)), new BCGaussianPrior( 193., quad(4., 4.)));
    m.setPrior(D->freeAmplitudes(rho,      rho)[0],    new BCGaussianPrior(0.157, quad(0.027, 0.020)), new BCGaussianPrior( 120., quad(7., 8.)));
    m.setPrior(D->freeAmplitudes(rho,      rho)[1],    new BCGaussianPrior(0.384, quad(0.020, 0.015)), new BCGaussianPrior( 163., quad(3., 3.)));
    m.setPrior(D->freeAmplitudes(rho,      rho)[2],    new BCGaussianPrior(0.624, quad(0.023, 0.015)), new BCGaussianPrior( 357., quad(3., 3.)));

    m.setPrior(D->freeAmplitudes(f_0_980, piPlus, piMinus)[0],    new BCGaussianPrior(0.233, quad(0.019, 0.015)), new BCGaussianPrior( 261., quad(7., 4.)));
    m.setPrior(D->freeAmplitudes(f_2,     piPlus, piMinus)[0],    new BCGaussianPrior(0.338, quad(0.021, 0.016)), new BCGaussianPrior( 317., quad(4., 4.)));
    m.setPrior(D->freeAmplitudes(sigma,   piPlus, piMinus)[0],    new BCGaussianPrior(0.432, quad(0.027, 0.022)), new BCGaussianPrior( 254., quad(4., 5.)));




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
