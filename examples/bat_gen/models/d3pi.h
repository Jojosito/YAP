// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D3PI__H
#define __BAT__D3PI__H

#include "../bat_fit.h"
#include "../fit_fitFraction.h"
#include "../tools.h"

#include <BreitWigner.h>
#include <Constants.h>
#include <container_utils.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <Filters.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <make_unique.h>
#include <Model.h>
#include <Parameter.h>
#include <Particle.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <PoleMass.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <Resonance.h>
#include <SpinAmplitudeCache.h>

#include <BAT/BCGaussianPrior.h>

#include <complex>
#include <memory>

using namespace std;
using namespace yap;

inline unique_ptr<Model> d3pi(unique_ptr<Model> M)
{
    auto F = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    M->setFinalState(piPlus, piMinus, piPlus);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);

    // rho
    auto rho = F.resonance(113, radialSize, make_shared<RelativisticBreitWigner>(775.49e-3, 149.4e-3));
    rho->addChannel(piPlus, piMinus);
    D->addChannel(rho, piPlus);

    // f_2(1270)
    auto f_2 = F.resonance(225, radialSize, make_shared<RelativisticBreitWigner>());
    f_2->addChannel(piPlus, piMinus);
    D->addChannel(f_2, piPlus);

    // f_0(980)
    auto f_0_980_flatte = make_shared<Flatte>(0.965);
    f_0_980_flatte->add(FlatteChannel(0.406, *piPlus, *piMinus));
    f_0_980_flatte->add(FlatteChannel(0.406 * 2, *F.fsp(321), *F.fsp(-321))); // K+K-
    auto f_0_980 = Resonance::create("f_0(980)", QuantumNumbers(0, 0), radialSize, f_0_980_flatte);
    f_0_980->addChannel(piPlus, piMinus);
    D->addChannel(f_0_980, piPlus);

    // f_0(1370)
    auto f_0_1370 = Resonance::create("f_0(1370)", F["f_0"], radialSize, make_unique<RelativisticBreitWigner>(1.350, 0.265));
    f_0_1370->addChannel(piPlus, piMinus);
    D->addChannel(f_0_1370, piPlus);

    // f_0(1500)
    auto f_0_1500 = F.resonance(F.pdgCode("f_0(1500)"), radialSize, make_unique<RelativisticBreitWigner>());
    f_0_1500->addChannel(piPlus, piMinus);
    D->addChannel(f_0_1500, piPlus);

    // sigma a.k.a. f_0(500)
    auto sigma = F.resonance(F.pdgCode("f_0(500)"), radialSize, make_unique<PoleMass>(complex<double>(0.470, -0.220)));
    sigma->addChannel(piPlus, piMinus);
    D->addChannel(sigma, piPlus);

    M->addInitialStateParticle(D);

    // Add channels to D
    *free_amplitude(*M, to(rho))      = polar(1., 0.);
    *free_amplitude(*M, to(f_0_980))  = polar(1.4, rad(12.));
    *free_amplitude(*M, to(f_2))      = polar(2.1, rad(-123.));
    *free_amplitude(*M, to(f_0_1370)) = polar(1.3, rad(-21.));
    *free_amplitude(*M, to(f_0_1500)) = polar(1.1, rad(-44.));
    *free_amplitude(*M, to(sigma))    = polar(3.7, rad(-3.));

    return M;
}

inline bat_fit d3pi_fit(string name, unique_ptr<Model> M, vector<vector<unsigned> > pcs = {})
{
    bat_fit m(name, d3pi(move(M)), pcs);

    auto is_rho = is_named("rho0");

    for (const auto& p : particles(*m.model(), is_resonance)) {
        auto fa = free_amplitude(*m.model(), to(p));
        if (is_rho(*p))
            m.fix(fa, real(fa->value()), imag(fa->value()));
        /* else */
        /*     m.setPrior(fa, 0.5 * abs(fa->value()), 1.5 * abs(fa->value()), */
        /*                deg(arg(fa->value())) - 15., deg(arg(fa->value())) + 15.); */
    }

    /* auto rho = dynamic_pointer_cast<Resonance>(particle(*m.model(), is_named("rho0"))); */
    /* m.addParameter("rho_mass", static_pointer_cast<BreitWigner>(rho->massShape())->mass(), 0.5, 1.2); */
    /* m.GetParameters().Back().SetPriorConstant(); */
    /* m.addParameter("rho_width", static_pointer_cast<BreitWigner>(rho->massShape())->width(), 0.1, 0.2); */
    /* m.GetParameters().Back().SetPriorConstant(); */
    // m.addParameter("rho_radialSize", rho->radialSize(), 1, 5);
    // m.GetParameters().Back().SetPriorConstant();
    
    return m;
}

inline fit_fitFraction d3pi_fit_fitFractions(string name, unique_ptr<Model> M, vector<vector<unsigned> > pcs = {})
{
    fit_fitFraction m(name, d3pi(move(M)), pcs);
    
    // find particles
    auto D = static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("D+")));
    auto rho      = static_pointer_cast<Resonance>(particle(*m.model(), is_named("rho0")));
    auto f_2      = static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_2")));
    auto f_0_980  = static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_0(980)")));
    auto f_0_1370 = static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_0(1370)")));
    auto f_0_1500 = static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_0(1500)")));
    auto sigma    = static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_0(500)")));

    // set fit fractions to fit
    m.setFitFraction(decay_tree(*D, to(rho)),      20e-2,   quad(2.3e-2, 0.9e-2));
    m.setFitFraction(decay_tree(*D, to(f_2)),      18.2e-2, quad(2.6e-2, 0.7e-2));
    m.setFitFraction(decay_tree(*D, to(f_0_980)),  4.1e-2,  quad(0.9e-2, 0.3e-2));
    m.setFitFraction(decay_tree(*D, to(f_0_1370)), 2.6e-2,  quad(1.8e-2, 0.6e-2));
    m.setFitFraction(decay_tree(*D, to(f_0_1500)), 3.4e-2,  quad(1.0e-2, 0.8e-2));
    m.setFitFraction(decay_tree(*D, to(sigma)),    41.8e-2, quad(1.4e-2, 2.5e-2));

    // set free amplitude parameters of fit
    m.fix(free_amplitude(*m.model(), to(rho)),      1.,  0.);
    m.fix(free_amplitude(*m.model(), to(f_2)),      2.1, -123.);
    m.fix(free_amplitude(*m.model(), to(f_0_980)),  1.4, 12.);
    m.fix(free_amplitude(*m.model(), to(f_0_1370)), 1.3, -21.);
    m.fix(free_amplitude(*m.model(), to(f_0_1500)), 1.1, -44.);
    m.fix(free_amplitude(*m.model(), to(sigma)),    3.7, -3);
    /* m.setPrior(free_amplitude(*m.model(), to(f_2)),      new BCGaussianPrior(2.1, quad(0.2, 0.1)), new BCGaussianPrior(-123., quad(6.,   3.))); */
    /* m.setPrior(free_amplitude(*m.model(), to(f_0_980)),  new BCGaussianPrior(1.4, quad(0.2, 0.2)), new BCGaussianPrior(  12., quad(12., 10.))); */
    /* m.setPrior(free_amplitude(*m.model(), to(f_0_1370)), new BCGaussianPrior(1.3, quad(0.4, 0.2)), new BCGaussianPrior( -21., quad(15., 14.))); */
    /* m.setPrior(free_amplitude(*m.model(), to(f_0_1500)), new BCGaussianPrior(1.1, quad(0.3, 0.2)), new BCGaussianPrior( -44., quad(13., 16.))); */
    /* m.setPrior(free_amplitude(*m.model(), to(sigma)),    new BCGaussianPrior(3.7, quad(0.3, 0.2)), new BCGaussianPrior(  -3., quad(4.,   2.))); */

    /* // add shape parameters */
    /* m.addParameter("f_0_980_mass", static_pointer_cast<Flatte>(f_0_980->massShape())->mass(), 0.953 - 3 * 0.02, 0.953 + 3 * 0.02); */
    /* m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.953, 0.02)); */
    /* /\* m.addParameter("f_0_980_coupling", f_0_980_flatte->channels()[0].Coupling, 0.329 - 3 * 0.096, 0.329 + 3 * 0.096); *\/ */
    /* /\* m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.329, 0.096)); *\/ */
    /* m.addParameter("f_0_1370_mass", static_pointer_cast<BreitWigner>(f_0_1370->massShape())->mass(), 1.259 - 3 * 0.055, 1.259 + 3 * 0.055); */
    /* m.GetParameters().Back().SetPrior(new BCGaussianPrior(1.259, 0.055)); */
    /* m.addParameter("f_0_1370_width", dynamic_pointer_cast<BreitWigner>(f_0_1370->massShape())->width(), */
    /*                0.298 - 3 * 0.021, 0.298 + 3 * 0.021); */
    /* m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.298, 0.021)); */
    /* m.addParameter("sigma_mass", dynamic_pointer_cast<PoleMass>(sigma->massShape())->mass(), */
    /*                complex<double>(0.466, -0.223) - 3. * complex<double>(0.018, 0.028), */
    /*                complex<double>(0.466, -0.223) + 3. * complex<double>(0.018, 0.028)); */
    /* m.GetParameters()[m.GetParameters().Size() - 2].SetPrior(new BCGaussianPrior(0.466, 0.018)); */
    /* m.GetParameters().Back().SetPrior(new BCGaussianPrior(-0.223, 0.028)); */

    return m;
}

    
#endif
