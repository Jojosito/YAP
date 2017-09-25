#include <catch.hpp>
#include <catch_capprox.hpp>

#include "../examples/bat_gen/models/d4pi.h"

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


TEST_CASE( "d4pia1_2" )
{

    // disable debug logs in test
    //yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    auto D_mass = F["D0"].mass();

    // create models
    for (bool rho_pi_S : {false, true})
        for (bool rho_pi_D : {false, true})
            for (bool pipiS_pi : {false, true}) {
                if ((rho_pi_S || rho_pi_D || pipiS_pi) == false)
                    continue;

                LOG(INFO) << "Waves used: " << (rho_pi_S ? "rho_pi_S; " : "")
                                << (rho_pi_D ? "rho_pi_D; " : "")
                                << (pipiS_pi ? "pipiS_pi; " : "");

                // which waves to include in the model
                const_cast<bool&>(d4pi_a_rho_pi_S)  = rho_pi_S;
                const_cast<bool&>(d4pi_a_rho_pi_D)  = rho_pi_D;
                const_cast<bool&>(d4pi_a_pipiS_pi)  = pipiS_pi;

                const_cast<bool&>(d4pi_pm_shared)   = false;

                const_cast<bool&>(d4pi_a1_plus )   = true;
                const_cast<bool&>(d4pi_a1_minus)   = false;

                auto M_plus = d4pi();

                // which waves to include in the model
                const_cast<bool&>(d4pi_a1_plus )   = false;
                const_cast<bool&>(d4pi_a1_minus)   = true;

                auto M_minus = d4pi();

                // which waves to include in the model
                const_cast<bool&>(d4pi_a1_plus )   = true;
                const_cast<bool&>(d4pi_a1_minus)   = true;

                const_cast<bool&>(d4pi_pm_shared)   = true;

                auto M_plusMinus = d4pi();

                // get pi pi mass ranges
                auto m2_pipi_range = yap::squared(yap::mass_range(D_mass, M_plus->massAxes(), M_plus->finalStateParticles()))[0];
                m2_pipi_range[0] *= 1.01;
                m2_pipi_range[1] *= 0.99;

                // create random number engine for generation of points
                std::mt19937 g(143425);


                for (unsigned i = 0; i < 100; ++i) {

                    auto P = yap::phsp(*M_plus, D_mass, M_plus->massAxes(), yap::squared(yap::mass_range(D_mass, M_plus->massAxes(), M_plus->finalStateParticles())), g, 1000);

                    if (P.empty()) {
                        continue;
                    }

                    const auto amp_plus =  calculate_model(D_mass, *M_plus,  M_plus->massAxes(), P);
                    const auto amp_minus = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[1], P[0], P[3], P[2]});
                    const auto amp_plusMinus =  calculate_model(D_mass, *M_plusMinus,  M_plus->massAxes(), P);

                    LOG(INFO) << amp_plus << "; " << amp_minus << "; " << amp_plusMinus << " \t"
                            << (amp_plus == Catch::Detail::CApprox(amp_minus) ? "S" :
                                    (-amp_plus == Catch::Detail::CApprox(amp_minus) ? "O" : "X"));

                    // can have slight differences due to bowler mass shape
                    if (d4pi_a1_bowler) {
                        REQUIRE ( fabs((real(amp_plus) - real(amp_minus)) / (real(amp_plus) + real(amp_minus))) < 0.01);
                        REQUIRE ( fabs((imag(amp_plus) - imag(amp_minus)) / (imag(amp_plus) + imag(amp_minus))) < 0.01);
                    }
                    else
                        REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_minus) );

                    auto amp_plus_swap = calculate_model(D_mass, *M_plus, M_plus->massAxes(), {P[2], P[3], P[0], P[1]});
                    REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_plus_swap) );
                    amp_plus_swap = calculate_model(D_mass, *M_plus, M_plus->massAxes(), {P[0], P[3], P[2], P[1]});
                    REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_plus_swap) );
                    amp_plus_swap = calculate_model(D_mass, *M_plus, M_plus->massAxes(), {P[0], P[1], P[2], P[3]});
                    REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_plus_swap) );

                    auto amp_minus_swap = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[3], P[2], P[1], P[0]});
                    REQUIRE ( amp_minus == Catch::Detail::CApprox(amp_minus_swap) );
                    amp_minus_swap = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[1], P[2], P[3], P[0]});
                    REQUIRE ( amp_minus == Catch::Detail::CApprox(amp_minus_swap) );
                    amp_minus_swap = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[3], P[0], P[1], P[2]});
                    REQUIRE ( amp_minus == Catch::Detail::CApprox(amp_minus_swap) );

                    auto amp_plusMinus_swap =  calculate_model(D_mass, *M_plusMinus,  M_plus->massAxes(), {P[2], P[3], P[0], P[1]});
                    REQUIRE ( amp_plusMinus == Catch::Detail::CApprox(amp_plusMinus_swap) );
                    amp_plusMinus_swap =  calculate_model(D_mass, *M_plusMinus,  M_plus->massAxes(), {P[1], P[0], P[3], P[2]});
                    REQUIRE ( amp_plusMinus == Catch::Detail::CApprox(amp_plusMinus_swap) );
                    amp_plusMinus_swap =  calculate_model(D_mass, *M_plusMinus,  M_plus->massAxes(), {P[0], P[3], P[2], P[1]});
                    REQUIRE ( amp_plusMinus == Catch::Detail::CApprox(amp_plusMinus_swap) );
                }
            }
}

TEST_CASE( "d4pia1_3" )
{

    // disable debug logs in test
    //yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    auto D_mass = F["D0"].mass();

    // create models
    auto M_plusMinus = d4pi();

    // which waves to include in the model
    const_cast<bool&>(d4pi_pm_shared)   = false;

    const_cast<bool&>(d4pi_a1_plus )   = true;
    const_cast<bool&>(d4pi_a1_minus)   = false;

    auto M_plus = d4pi();

    // which waves to include in the model
    const_cast<bool&>(d4pi_a1_plus )   = false;
    const_cast<bool&>(d4pi_a1_minus)   = true;

    auto M_minus = d4pi();


    // get pi pi mass ranges
    auto m2_pipi_range = yap::squared(yap::mass_range(D_mass, M_plus->massAxes(), M_plus->finalStateParticles()))[0];
    m2_pipi_range[0] *= 1.01;
    m2_pipi_range[1] *= 0.99;

    // create random number engine for generation of points
    std::mt19937 g(143425);


    for (unsigned i = 0; i < 100; ++i) {

        auto P = yap::phsp(*M_plus, D_mass, M_plus->massAxes(), yap::squared(yap::mass_range(D_mass, M_plus->massAxes(), M_plus->finalStateParticles())), g, 1000);

        if (P.empty()) {
            continue;
        }

        const auto amp_plus =  calculate_model(D_mass, *M_plus,  M_plus->massAxes(), P);
        const auto amp_minus = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[1], P[0], P[3], P[2]});
        const auto amp_plusMinus =  calculate_model(D_mass, *M_plusMinus,  M_plus->massAxes(), P);

        LOG(INFO) << amp_plus << "; " << amp_minus << "; " << amp_plusMinus << " \t"
        << (amp_plus == Catch::Detail::CApprox(amp_minus) ? "S" :
                (-amp_plus == Catch::Detail::CApprox(amp_minus) ? "O" : "X"));

        // can have slight differences due to bowler mass shape
        if (d4pi_a1_bowler) {
            REQUIRE ( fabs((real(amp_plus) - real(amp_minus)) / (real(amp_plus) + real(amp_minus))) < 0.01);
            REQUIRE ( fabs((imag(amp_plus) - imag(amp_minus)) / (imag(amp_plus) + imag(amp_minus))) < 0.01);
        }
        else
            REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_minus) );

        auto amp_plus_swap = calculate_model(D_mass, *M_plus, M_plus->massAxes(), {P[2], P[3], P[0], P[1]});
        REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_plus_swap) );
        amp_plus_swap = calculate_model(D_mass, *M_plus, M_plus->massAxes(), {P[0], P[3], P[2], P[1]});
        REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_plus_swap) );
        amp_plus_swap = calculate_model(D_mass, *M_plus, M_plus->massAxes(), {P[0], P[1], P[2], P[3]});
        REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_plus_swap) );

        auto amp_minus_swap = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[3], P[2], P[1], P[0]});
        REQUIRE ( amp_minus == Catch::Detail::CApprox(amp_minus_swap) );
        amp_minus_swap = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[1], P[2], P[3], P[0]});
        REQUIRE ( amp_minus == Catch::Detail::CApprox(amp_minus_swap) );
        amp_minus_swap = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[3], P[0], P[1], P[2]});
        REQUIRE ( amp_minus == Catch::Detail::CApprox(amp_minus_swap) );

        auto amp_plusMinus_swap =  calculate_model(D_mass, *M_plusMinus,  M_plus->massAxes(), {P[2], P[3], P[0], P[1]});
        REQUIRE ( amp_plusMinus == Catch::Detail::CApprox(amp_plusMinus_swap) );
        amp_plusMinus_swap =  calculate_model(D_mass, *M_plusMinus,  M_plus->massAxes(), {P[1], P[0], P[3], P[2]});
        REQUIRE ( amp_plusMinus == Catch::Detail::CApprox(amp_plusMinus_swap) );
        amp_plusMinus_swap =  calculate_model(D_mass, *M_plusMinus,  M_plus->massAxes(), {P[0], P[3], P[2], P[1]});
        REQUIRE ( amp_plusMinus == Catch::Detail::CApprox(amp_plusMinus_swap) );
    }
}

/*
TEST_CASE( "d4pi_components" )
{

    // disable debug logs in test
    //yap::disableLogs(el::Level::Debug);
    yap::plainLogs(el::Level::Debug);

    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    auto D_mass = F["D0"].mass();

    // create models
    unsigned i_component(0);
    while (true) {

        unsigned i_check(0);

        // which waves to include in the model
        const_cast<bool&>(a_rho_pi_S)  = false;
        const_cast<bool&>(a_rho_pi_D)  = false;
        const_cast<bool&>(a_sigma_pi)  = false;
        const_cast<bool&>(rho_rho  )  = i_component == i_check++;
        const_cast<bool&>(f_0_pipi  )  = i_component == i_check++;
        const_cast<bool&>(f_2_pipi  )  = i_component == i_check++;
        const_cast<bool&>(sigma_pipi)  = i_component == i_check++;

        const_cast<bool&>(omega_omega) = i_component == i_check++;
        const_cast<bool&>(a1_omega   ) = false;

        // CLEO waves
        const_cast<bool&>(sigma_f_0_1370) = i_component == i_check++;
        const_cast<bool&>(sigma_rho) = i_component == i_check++; // problem
        const_cast<bool&>(f_2_f_2  ) = i_component == i_check++;

        const_cast<bool&>(pi1300_pi_pi_pi) = i_component == i_check++;


        const_cast<bool&>(flat_4pi)    = i_component == i_check++;


        const_cast<bool&>(bg_flat_4pi) = i_component == i_check++;
        const_cast<bool&>(bg_rho     ) = i_component == i_check++; // problem
        const_cast<bool&>(bg_a1      ) = i_component == i_check++;

        const_cast<bool&>(a1_bowler)  = false;
        const_cast<bool&>(a1_shared)  = true;  // share a1+ and a1- free amplitudes

        const_cast<bool&>(a1_plus )    = false;
        const_cast<bool&>(a1_minus)    = false;

        const_cast<bool&>(bg_only)     = false; // fix D admixture to 0

        if (++i_component > i_check)
            break;

        std::unique_ptr<Model> M_plus;
        try {
            M_plus = d4pi();
        }
        catch (yap::exceptions::Exception& e) {
            continue;
        }

        if (M_plus->components().empty())
            continue;


        auto M_minus = d4pi();

        // get pi pi mass ranges
        auto m2_pipi_range = yap::squared(yap::mass_range(D_mass, M_plus->massAxes(), M_plus->finalStateParticles()))[0];
        m2_pipi_range[0] *= 1.01;
        m2_pipi_range[1] *= 0.99;

        // create random number engine for generation of points
        std::mt19937 g(143425);


        for (unsigned i = 0; i < 10; ++i) {

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

            REQUIRE ( (amp_plus == Catch::Detail::CApprox(amp_minus) or amp_plus == Catch::Detail::CApprox(-amp_minus)) );

            auto amp_plus_swap = calculate_model(D_mass, *M_plus, M_plus->massAxes(), {P[2], P[3], P[0], P[1]});
            REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_plus_swap) );
            amp_plus_swap = calculate_model(D_mass, *M_plus, M_plus->massAxes(), {P[0], P[3], P[2], P[1]});
            REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_plus_swap) );
            amp_plus_swap = calculate_model(D_mass, *M_plus, M_plus->massAxes(), {P[0], P[1], P[2], P[3]});
            REQUIRE ( amp_plus == Catch::Detail::CApprox(amp_plus_swap) );

            auto amp_minus_swap = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[3], P[2], P[1], P[0]});
            REQUIRE ( amp_minus == Catch::Detail::CApprox(amp_minus_swap) );
            amp_minus_swap = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[1], P[2], P[3], P[0]});
            REQUIRE ( amp_minus == Catch::Detail::CApprox(amp_minus_swap) );
            amp_minus_swap = calculate_model(D_mass, *M_minus, M_plus->massAxes(), {P[3], P[0], P[1], P[2]});
            REQUIRE ( amp_minus == Catch::Detail::CApprox(amp_minus_swap) );
        }
    }

}
*/
