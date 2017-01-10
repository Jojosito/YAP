// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__A1__H
#define __BAT__A1__H

#include "d4pi.h"

#include "../bat_fit.h"
#include "../tools.h"

#include <a1MassShape.h>
#include <Attributes.h>
#include <AmplitudeBasis.h>
#include <BreitWigner.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MathUtilities.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <SpinAmplitudeCache.h>


#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCConstantPrior.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>


#include <complex>
#include <memory>

using namespace yap;

//-------------------------
inline std::unique_ptr<Model> a1()
{
    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus  = FinalStateParticle::create(T[211]);
    auto piMinus = FinalStateParticle::create(T[-211]);

    auto M = std::make_unique<yap::Model>(std::make_unique<yap::HelicityFormalism>());
    M->setFinalState(piPlus, piMinus, piPlus);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto a_1 = DecayingParticle::create(T["a_1+"], radialSize);
    
    //
    // resonant particles
    //
    
    // rho
    auto rho = DecayingParticle::create(T["rho0"], radialSize, std::make_shared<RelativisticBreitWigner>(T["rho0"]));
    rho->addStrongDecay(piPlus, piMinus);

    // omega
    //auto omega = DecayingParticle::create(T["omega"], radialSize, std::make_shared<BreitWigner>(T["omega"]));
    //omega->addStrongDecay({piPlus, piMinus});
    
    // sigma / f_0(500)
    auto sigma = DecayingParticle::create(T["f_0(500)"], radialSize, std::make_shared<BreitWigner>(T["f_0(500)"]));
    sigma->addStrongDecay(piPlus, piMinus);
    

    // a_1 -> sigma pi 
    a_1->addStrongDecay(sigma, piPlus);
    *free_amplitude(*a_1, to(sigma)) = std::polar(scale_a_rho_sigma * 0.439, rad(193.));

    // a_1 -> rho pi
    a_1->addStrongDecay(rho,   piPlus);
    // S wave
    *free_amplitude(*a_1, to(rho), l_equals(0)) = 1;
    // D wave
    *free_amplitude(*a_1, to(rho), l_equals(2)) = std::polar(scale_a_rho_pi_D * 0.241, rad(82.));


    /*LOG(INFO) << "a_1 decay trees:";
    LOG(INFO) << to_string(a_1->decayTrees());

    LOG(INFO) << std::endl << "Free amplitudes: ";
    for (const auto& fa : free_amplitudes(*M, yap::is_not_fixed()))
        LOG(INFO) << yap::to_string(*fa) << "  \t (mag, phase) = (" << abs(fa->value()) << ", " << deg(arg(fa->value())) << "°)"
            << "  \t (real, imag) = (" << real(fa->value()) << ", " << imag(fa->value()) << ")";
     */

    return M;
}

//-------------------------
inline bat_fit a1_fit()
{
    bat_fit m("a1_fit", a1());

    return m;
}

#endif
