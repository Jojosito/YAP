#ifndef __test_helper_functions__H
#define __test_helper_functions__H

#include <Attributes.h>
#include <BreitWigner.h>
#include <DataSet.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <PHSP.h>

#include <Group.h>
#include <Sort.h>

/// generate a model with 4 final state particles
//-------------------------
inline std::shared_ptr<yap::Model> d4pi(std::vector<int> fsps = {211, -211, 211, -211})
{
    auto M = std::make_shared<yap::Model>(std::make_unique<yap::HelicityFormalism>());

    auto T = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");
    double radialSize = 1.;

    // initial state particle
    auto D = yap::DecayingParticle::create(T[421], radialSize);

    // final state particles
    yap::FinalStateParticleVector FSP;
    std::transform(fsps.begin(), fsps.end(), std::back_inserter(FSP), [&T](int pdg){return yap::FinalStateParticle::create(T[pdg]);});


    auto piPlus  = yap::FinalStateParticle::create(T[211]);
    auto piMinus = yap::FinalStateParticle::create(T[-211]);

    // Set final-state particles
    std::vector<std::shared_ptr<yap::FinalStateParticle>> FSPs;
    for (auto i : fsps) {
        DEBUG(i <<  "  ");
        FSPs.push_back(i > 0 ? piPlus : piMinus);
    }

    M->setFinalState(FSPs);

    // rho
    auto rho = yap::DecayingParticle::create(T[113], radialSize, std::make_shared<yap::BreitWigner>(T[113]));
    rho->addStrongDecay(piPlus, piMinus);

    // sigma / f_0(500)
    auto sigma = yap::DecayingParticle::create(T[9000221], radialSize, std::make_shared<yap::BreitWigner>(T[9000221]));
    sigma->addStrongDecay(piPlus, piMinus);

    // a_1
    auto a_1 = yap::DecayingParticle::create(T[20213], radialSize, std::make_shared<yap::BreitWigner>(T[20213]));
    a_1->addStrongDecay(sigma, piPlus);
    a_1->addStrongDecay(rho,   piPlus);

    // D's channels
    D->addWeakDecay(rho, rho);
    D->addWeakDecay(a_1, piMinus);

    M->lock();
    return M;
}

/// generate a model with 4 final state particles
//-------------------------
inline std::shared_ptr<yap::Model> d4pi_a1(bool rho_pi_S, bool rho_pi_D, bool sigma_pi, bool plusMinus)
{
    auto M = std::make_shared<yap::Model>(std::make_unique<yap::HelicityFormalism>());

    auto T = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");
    double radialSize = 1.;

    // initial state particle
    auto D = yap::DecayingParticle::create(T[plusMinus ? 421 : -421], radialSize);

    // final state particles
    auto piPlus  = yap::FinalStateParticle::create(T[211]);
    auto piMinus = yap::FinalStateParticle::create(T[-211]);

    // Set final-state particles
    M->setFinalState(piPlus, piMinus, piPlus, piMinus);

    // rho
    auto rho = yap::DecayingParticle::create(T[113], radialSize, std::make_shared<yap::BreitWigner>(T[113]));
    rho->addStrongDecay(piPlus, piMinus);

    // sigma / f_0(500)
    auto sigma = yap::DecayingParticle::create(T[9000221], radialSize, std::make_shared<yap::BreitWigner>(T[9000221]));
    sigma->addStrongDecay(piPlus, piMinus);

    // a_1
    auto a_1_plus = yap::DecayingParticle::create(T["a_1+"], radialSize, std::make_shared<yap::BreitWigner>(T["a_1+"]));
    auto a_1_minus = yap::DecayingParticle::create(T["a_1-"], radialSize, std::make_shared<yap::BreitWigner>(T["a_1-"]));

    if (rho_pi_S or rho_pi_D) {
        a_1_plus->addStrongDecay(rho,   piPlus);
        a_1_minus->addStrongDecay(rho,   piMinus);

        if (not rho_pi_S) {
            *free_amplitude(*a_1_plus, yap::to(rho), yap::l_equals(0)) = 0.;
            *free_amplitude(*a_1_minus, yap::to(rho), yap::l_equals(0)) = 0.;
        }

        if (not rho_pi_D) {
            *free_amplitude(*a_1_plus, yap::to(rho), yap::l_equals(2)) = 0.;
            *free_amplitude(*a_1_minus, yap::to(rho), yap::l_equals(2)) = 0.;
        }

    }

    if (sigma_pi) {
        a_1_plus->addStrongDecay(sigma, piPlus);
        a_1_minus->addStrongDecay(sigma, piMinus);
    }

    // D's channels
    if (plusMinus)
        D->addWeakDecay(a_1_plus, piMinus);
    else
        D->addWeakDecay(a_1_minus, piPlus);
        //D->addWeakDecay(piPlus, a_1_minus);  // does not make a difference

    M->lock();
    return M;
}

//-------------------------
template <typename Formalism>
inline std::shared_ptr<yap::Model> dkkp(int pdg_D, std::vector<int> fsps)
{
    auto T = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    yap::FinalStateParticleVector FSP;
    std::transform(fsps.begin(), fsps.end(), std::back_inserter(FSP), [&T](int pdg){return yap::FinalStateParticle::create(T[pdg]);});

    double radial_size = 3;

    auto M = std::make_shared<yap::Model>(std::make_unique<Formalism>());
    M->setFinalState(FSP);
    
    auto D = yap::DecayingParticle::create(T[pdg_D], radial_size);

    auto piPlus = lone_elt(filter(FSP, yap::is_named("pi+")));
    auto kPlus  = lone_elt(filter(FSP, yap::is_named("K+")));
    auto kMinus = lone_elt(filter(FSP, yap::is_named("K-")));
        
    for (unsigned j = 0; j < 3; ++j) {
        auto res = yap::DecayingParticle::create("res_" + std::to_string(j), yap::QuantumNumbers(0, j * 2), radial_size,
                                                 std::make_shared<yap::BreitWigner>(0.750 + 0.250 * j, 0.025));
        res->addStrongDecay(piPlus, kMinus);
        D->addWeakDecay(res, kPlus);
    }

    *free_amplitude(*M, yap::from(*D), yap::l_equals(0)) = 0.5;
    *free_amplitude(*M, yap::from(*D), yap::l_equals(1)) = 1.;
    *free_amplitude(*M, yap::from(*D), yap::l_equals(2)) = 30.;

    M->lock();
    return M;
}

//-------------------------
template <typename Formalism>
std::shared_ptr<yap::Model> d3pi()
{
    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    auto M = std::make_shared<yap::Model>(std::make_unique<Formalism>());

    auto T = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    // initial state particle
    auto D = yap::DecayingParticle::create(T["D+"], radialSize);

    // final state particles
    auto piPlus  = yap::FinalStateParticle::create(T[211]);
    auto piMinus = yap::FinalStateParticle::create(T[-211]);

    // set final state
    M->setFinalState(piPlus, piMinus, piPlus);

    // rho
    auto rho = yap::DecayingParticle::create(T["rho0"], radialSize, std::make_shared<yap::BreitWigner>(T["rho0"]));
    rho->addStrongDecay(piPlus, piMinus);

    // Add channels to D
    D->addWeakDecay(rho, piPlus);

    M->lock();
    return M;
}

//-------------------------
inline yap::DataSet generate_data(yap::Model& M, unsigned nPoints)
{
    auto T = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");
    
    auto isp_mass = T[M.initialStates()[0]->name()].mass();

    auto A = M.massAxes();
    auto m2r = yap::squared(yap::mass_range(isp_mass, A, M.finalStateParticles()));

    yap::DataSet data(M.createDataSet());

    std::mt19937 g(0);
    // fill data set with nPoints points
    std::generate_n(std::back_inserter(data), nPoints,
                    std::bind(yap::phsp<std::mt19937>, std::cref(M), isp_mass, A, m2r, g, std::numeric_limits<unsigned>::max()));

    return data;
}



#endif

