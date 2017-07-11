#include <catch.hpp>
#include <catch_capprox.hpp>

#include "helperFunctions.h"

#include <Attributes.h>
#include <BreitWigner.h>
#include <container_utils.h>
#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayTree.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <MassRange.h>
#include <Model.h>
#include <Parameter.h>
#include <PDL.h>
#include <ZemachFormalism.h>

#include <cmath>

/**
 * Test that the amplitude remains the same after swapping the order of the final state particles
 */
inline std::complex<double> calculate_model(double isp_mass, yap::Model& M, const yap::MassAxes& A, std::vector<yap::FourVector<double>> P)
{
    auto data = M.createDataSet();
    // add point
    data.push_back(P);

    // return amplitude
    M.calculate(data);
    return amplitude(M.components()[0].decayTrees(), data[0]);
}


TEST_CASE( "d4pia1" )
{

    // disable debug logs in test
    //yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    auto D_mass = F["D0"].mass();

    // create models
    for (bool rho_pi_S : {false, true}) {
        for (bool rho_pi_D : {false, true}) {
            for (bool sigma_pi : {false, true}) {
                if ((rho_pi_S || rho_pi_D || sigma_pi) == false)
                    continue;

                DEBUG((rho_pi_S ? "rho_pi_S; " : "")
                        << (rho_pi_D ? "rho_pi_D; " : "")
                        << (sigma_pi ? "sigma_pi; " : "") );


                auto M_plus = d4pi_a1(rho_pi_S, rho_pi_D, sigma_pi, true);
                auto M_minus = d4pi_a1(rho_pi_S, rho_pi_D, sigma_pi, false);

                // get pi pi mass ranges
                auto m2_pipi_range = yap::squared(yap::mass_range(D_mass, M_plus->massAxes(), M_plus->finalStateParticles()))[0];
                m2_pipi_range[0] *= 1.01;
                m2_pipi_range[1] *= 0.99;

                // create random number engine for generation of points
                std::mt19937 g(0);


                for (unsigned i = 0; i < 100; ++i) {

                    auto P = yap::phsp(*M_plus, D_mass, M_plus->massAxes(), yap::squared(yap::mass_range(D_mass, M_plus->massAxes(), M_plus->finalStateParticles())), g, 1000);

                    if (P.empty()) {
                        continue;
                    }

                    DEBUG("\namp_plus");
                    auto amp_plus =  calculate_model(D_mass, *M_plus,  M_plus->massAxes(), P);
                    DEBUG("\namp_minus");
                    auto amp_minus = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[1], P[0], P[3], P[2]});


                    DEBUG(amp_plus << "; " << amp_minus << " \t"
                            << (amp_plus == Catch::Detail::CApprox(amp_minus) ? "S" :
                                    (-amp_plus == Catch::Detail::CApprox(amp_minus) ? "O" : "X")));

                    if (amp_plus != amp_plus)
                        REQUIRE(amp_minus != amp_minus);
                    else
                        REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_minus) );
                }
            }
        }
    }

}

