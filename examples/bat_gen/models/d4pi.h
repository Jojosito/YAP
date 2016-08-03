// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D4PI__H
#define __BAT__D4PI__H

#include "../fit_fitFraction.h"

#include <AmplitudeBasis.h>
#include <Constants.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <Filters.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <Resonance.h>
#include <SpinAmplitudeCache.h>

#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <complex>
#include <memory>

using namespace yap;


const double quad(std::vector<double> S)
{ return sqrt(std::accumulate(S.begin(), S.end(), 0., [](double a, double s) {return a + s * s;})); }

template <typename ... Types>
constexpr double quad(double s0, Types ... additional)
{ return quad({s0, additional...}); }


inline std::unique_ptr<Model> d4pi_model()
{
    auto F = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    auto M = std::make_unique<Model>(std::make_unique<HelicityFormalism>());
    M->setFinalState(piPlus, piMinus, piPlus, piMinus);

    // use common radial size for all resonances
    double radialSize = 1.2; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D0"), radialSize);
    

    //
    // resonant particles
    //
    
    // rho
    auto rho = F.resonance(F.pdgCode("rho0"), radialSize, std::make_shared<RelativisticBreitWigner>());
    rho->addChannel(piPlus, piMinus);

    // omega
    //auto omega = F.resonance(F.pdgCode("omega"), radialSize, std::make_shared<BreitWigner>());
    //omega->addChannel({piPlus, piMinus});
    
    // sigma / f_0(500)
    auto sigma = F.resonance(F.pdgCode("f_0(500)"), radialSize, std::make_shared<RelativisticBreitWigner>());
    sigma->addChannel(piPlus, piMinus);

    // a_1
    auto a_1 = F.resonance(F.pdgCode("a_1+"), radialSize, std::make_shared<BreitWigner>());
    a_1->addChannel(rho,   piPlus);
    a_1->addChannel(sigma, piPlus);
    
    // f_0(980) (as Flatte)
    auto piZero = F.fsp(111);
    auto Kshort = F.fsp(310);
    auto f_0_980_flatte = std::make_shared<Flatte>();
    f_0_980_flatte->addChannel(0.20, piZero->mass()->value());
    f_0_980_flatte->addChannel(0.50, Kshort->mass()->value());
    auto f_0_980 = F.resonance(F.pdgCode("f_0"), radialSize, f_0_980_flatte);
    f_0_980->addChannel(piPlus, piMinus);
       
    // f_2(1270)
    auto f_2 = F.resonance(F.pdgCode("f_2"), radialSize, std::make_shared<RelativisticBreitWigner>());
    f_2->addChannel(piPlus, piMinus); 
    
    // pi+ pi- flat
    auto pipiFlat = F.decayingParticle(F.pdgCode("pi0"), radialSize); // just need any spin-0 particle
    pipiFlat->addChannel(piPlus, piMinus);   
    

    // D0 channels
    D->addChannel(rho, rho);
    D->addChannel(a_1, piMinus);
    D->addChannel(f_0_980, piPlus, piMinus);
    D->addChannel(f_2, pipiFlat);
    D->addChannel(sigma, piPlus, piMinus);
    
    M->addInitialStateParticle(D);


    ///
    /// amplitudes
    ///

    // a_1 -> rho pi
    *free_amplitude(*M, from(a_1), to(rho), l_equals(0)) = 1; // S wave
     free_amplitude(*M, from(a_1), to(rho), l_equals(0))->variableStatus() = VariableStatus::fixed;
    *free_amplitude(*M, from(a_1), to(rho), l_equals(1)) = 0.; // P wave
    *free_amplitude(*M, from(a_1), to(rho), l_equals(2)) = std::polar(0.241, rad(82.)); // D wave

    // a_1 -> sigma pi
    *free_amplitude(*M, from(a_1), to(sigma)) = std::polar(0.439, rad(193.));

    // R pi pi
    *free_amplitude(*M, to(f_0_980)) = std::polar(0.233, rad(261.));
    *free_amplitude(*M, to(f_2    )) = std::polar(0.338, rad(317.));
    *free_amplitude(*M, to(sigma  )) = std::polar(0.432, rad(254.));
    
    // rho rho
    // transform into angular momentum basis
    basis::canonical<double> c(basis::transversity<double>(
                                   std::polar(0.624, rad(357.)),    // A_longitudinal
                                   std::polar(0.157, rad(120.)),    // A_parallel
                                   std::polar(0.384, rad(163.)) )); // A_perpendicular

    for (unsigned l = 0; l < 3; ++l)
        *free_amplitude(*M, to(rho, rho), l_equals(l)) = c[l];

    
    return M;
}

inline fit_fitFraction d4pi_fit_fitFraction()
{
    fit_fitFraction m("D4pi_frac_fit", d4pi_model());

    auto& M = m.model();

    auto D       = lone_elt(M->initialStateParticles()).first;
    auto piPlus  = lone_elt(particles(*M, is_named("pi+")));
    auto piMinus = lone_elt(particles(*M, is_named("pi-")));
    auto a_1     = lone_elt(particles(*M, is_named("a_1+")));
    auto rho     = lone_elt(particles(*M, is_named("rho0")));
    auto f_0_980 = lone_elt(particles(*M, is_named("f_0")));
    auto f_2     = lone_elt(particles(*M, is_named("f_2")));
    auto sigma   = lone_elt(particles(*M, is_named("f_0(500)")));

    // set free amplitude parameters of fit
    m.fix     (free_amplitude(*M, from(a_1), to(rho), l_equals(0)), 1., 0.);
    m.setPrior(free_amplitude(*M, from(a_1), to(rho), l_equals(2)), new BCGaussianPrior(0.241, quad(0.033, 0.024)), new BCGaussianPrior( 82., quad(5., 4.)));
    m.setPrior(free_amplitude(*M, from(a_1), to(sigma)),            new BCGaussianPrior(0.439, quad(0.026, 0.021)), new BCGaussianPrior(193., quad(4., 4.)));

    // todo transform amplitudes and errors
    m.setPrior(free_amplitude(*M, to(rho), l_equals(0)), new BCGaussianPrior(0.157, quad(0.027, 0.020)), new BCGaussianPrior(120., quad(7., 8.)));
    m.setPrior(free_amplitude(*M, to(rho), l_equals(1)), new BCGaussianPrior(0.384, quad(0.020, 0.015)), new BCGaussianPrior(163., quad(3., 3.)));
    m.setPrior(free_amplitude(*M, to(rho), l_equals(2)), new BCGaussianPrior(0.624, quad(0.023, 0.015)), new BCGaussianPrior(357., quad(3., 3.)));

    m.setPrior(free_amplitude(*M, to(f_0_980)), new BCGaussianPrior(0.233, quad(0.019, 0.015)), new BCGaussianPrior(261., quad(7., 4.)));
    m.setPrior(free_amplitude(*M, to(f_2)),     new BCGaussianPrior(0.338, quad(0.021, 0.016)), new BCGaussianPrior(317., quad(4., 4.)));
    m.setPrior(free_amplitude(*M, to(sigma)),   new BCGaussianPrior(0.432, quad(0.027, 0.022)), new BCGaussianPrior(254., quad(4., 5.)));

    // set fit fractions to fit
    m.setFitFraction(lone_elt(D->decayTreeSet(to(a_1), l_equals(0))), 43.3e-2, quad(2.5e-2, 1.9e-2)); // S
    m.setFitFraction(lone_elt(D->decayTreeSet(to(a_1), l_equals(2))),  2.5e-2, quad(0.5e-2, 0.4e-2)); // D
    m.setFitFraction(lone_elt(D->decayTreeSet(to(a_1), l_equals(1))),  8.3e-2, quad(0.7e-2, 0.6e-2)); // sigma

    m.setFitFraction(lone_elt(D->decayTreeSet(to(rho), l_equals(0))),  1.1e-2, quad(0.3e-2, 0.3e-2));
    m.setFitFraction(lone_elt(D->decayTreeSet(to(rho), l_equals(1))),  6.4e-2, quad(0.6e-2, 0.5e-2));
    m.setFitFraction(lone_elt(D->decayTreeSet(to(rho), l_equals(2))), 16.8e-2, quad(1.0e-2, 0.8e-2));

    m.setFitFraction(lone_elt(D->decayTreeSet(to(f_0_980))), 2.4e-2, quad(0.5e-2, 0.4e-2));
    m.setFitFraction(lone_elt(D->decayTreeSet(to(f_2))),     4.9e-2, quad(0.6e-2, 0.5e-2));
    m.setFitFraction(lone_elt(D->decayTreeSet(to(sigma))),   8.2e-2, quad(0.9e-2, 0.7e-2));

    return m;
}

#endif
