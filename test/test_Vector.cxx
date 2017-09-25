#include <catch.hpp>

#include <FourVector.h>
#include <ThreeVector.h>
#include <logging.h>
#include <MathUtilities.h>

#include <cmath>

TEST_CASE( "Vector" )
{

    yap::CoordinateSystem<double, 3> three_axes({yap::ThreeVector<double>({1., 0., 0.}),
                                                 yap::ThreeVector<double>({0., 1., 0.}),
                                                 yap::ThreeVector<double>({0., 0., 1.})});
    
    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    SECTION( "Initialization" ) {
        yap::Vector<double, 3> v;
        yap::Vector<double, 3> v0({0, 0, 0});
        REQUIRE( v == v0);
    }

    SECTION( "ThreeVector" ) {

        const auto zero = yap::ThreeVector<double>({0, 0, 0});
        const auto v1 = yap::ThreeVector<double>({1, 2, 3});
        const auto v2 = yap::ThreeVector<double>({4, 5, 6});

        SECTION( "addition-assignment" ) {
            yap::ThreeVector<double> v1_copy(v1);
            v1_copy += v2;
            REQUIRE( v1_copy == yap::ThreeVector<double>({5, 7, 9}) );
            REQUIRE( not (v1_copy == zero) );
        }


        SECTION( "subtraction-assignment" ) {
            yap::ThreeVector<double> v1_copy(v1);
            v1_copy -= v2;
            REQUIRE( v1_copy == yap::ThreeVector<double>({ -3, -3, -3}) );
            REQUIRE( not (v1_copy == zero) );
        }

        SECTION( "multiplication-assignment" ) {
            yap::ThreeVector<double> v1_copy(v1);
            v1_copy *= 3.;
            REQUIRE( v1_copy == yap::ThreeVector<double>({3, 6, 9}) );
            REQUIRE( not (v1_copy == zero) );
        }


        SECTION( "arithmetic operations" ) {

            // +
            REQUIRE( v1 + v2 == yap::ThreeVector<double>({5, 7, 9}) );

            // -
            REQUIRE( v1 - v2 == yap::ThreeVector<double>({ -3, -3, -3}) );

            // inner product
            REQUIRE( v1 * v2 == 32 );

            // norm
            REQUIRE( norm(v1) == 14 );
            REQUIRE( norm(v2) == 77 );

            // abs
            REQUIRE( abs(v1) == sqrt(14) );
            REQUIRE( abs(v2) == sqrt(77) );

            // cross product
            auto x = yap::ThreeVector<double>({1, 0, 0});
            auto y = yap::ThreeVector<double>({0, 1, 0});
            auto z = yap::ThreeVector<double>({0, 0, 1});
            REQUIRE( cross(x, y) == z );

            // (unary) minus
            REQUIRE( -v1 == -1. * v1 );
            REQUIRE( -v1 + v1 == zero );


        }

        SECTION( "angles" ) {

            auto a = yap::ThreeVector<double>({2, 0, 0});
            auto b = yap::ThreeVector<double>({0, 2, 0});
            auto c = yap::ThreeVector<double>({0, 0, 2});
            auto z = yap::ThreeVector<double>({0, 0, 0});

            REQUIRE(angle(a, a) == 0.);
            REQUIRE(angle(a, -a) == yap::pi());
            REQUIRE(angle(a, b) == 0.5 * yap::pi());
            REQUIRE(angle(a, c) == 0.5 * yap::pi());
            REQUIRE(angle(a, z) == 0.);

            REQUIRE( yap::theta(c, three_axes) == Approx(0.) );
            REQUIRE( yap::theta(a, three_axes) == Approx(yap::pi() / 2.) );
            REQUIRE( yap::theta(b, three_axes) == Approx(yap::pi() / 2.) );

            REQUIRE( yap::phi(a, three_axes) == Approx(0.) );
            REQUIRE( yap::phi(b, three_axes) == Approx(yap::pi() / 2.) );
        }

        SECTION( "cross product" ) {
            REQUIRE( cross(three_axes[0], three_axes[1]) == three_axes[2]);
            REQUIRE( cross(three_axes[1], three_axes[0]) == -three_axes[2]);

            REQUIRE( cross(three_axes[1], three_axes[2]) == three_axes[0]);
            REQUIRE( cross(three_axes[2], three_axes[1]) == -three_axes[0]);

            REQUIRE( cross(three_axes[2], three_axes[0]) == three_axes[1]);
            REQUIRE( cross(three_axes[0], three_axes[2]) == -three_axes[1]);
        }

        SECTION( "constants" ) {

            // axes:
            REQUIRE( cross(three_axes[0], three_axes[1]) == three_axes[2] );
            REQUIRE( is_right_handed(three_axes) );
        }

    }

    SECTION( "FourVector" ) {

        const auto v1 = yap::FourVector<double>({4, 3, 2, 1});
        const auto v2 = yap::FourVector<double>({8, 7, 6, 5});

        SECTION( "addition-assignment" ) {
            yap::FourVector<double> v1_copy(v1);
            v1_copy += v2;
            REQUIRE( v1_copy == yap::FourVector<double>({12, 10, 8, 6}) );
            REQUIRE( not (v1_copy == v1) );
        }

        SECTION( "subtraction-assignment" ) {
            yap::FourVector<double> v1_copy(v1);
            v1_copy -= v2;
            REQUIRE( v1_copy == yap::FourVector<double>({ -4, -4, -4, -4}) );
            REQUIRE( not (v1_copy == v1) );
        }

        SECTION( "multiplication-assignment" ) {
            yap::FourVector<double> v1_copy(v1);
            v1_copy *= 2.;
            REQUIRE( v1_copy == yap::FourVector<double>({8, 6, 4, 2}) );
            REQUIRE( not (v1_copy == v1) );
        }

        SECTION( "arithmetic operations" ) {

            // +
            REQUIRE( v1 + v2 == yap::FourVector<double>({12, 10, 8, 6}) );

            // -
            REQUIRE( v1 - v2 == yap::FourVector<double>({ -4, -4, -4, -4}) );

            // * (four-vector inner product)
            REQUIRE( v1 * v2 == -6 );

            // norm
            REQUIRE( norm(v1) == 2 );
            REQUIRE( norm(v2) == -46 );

            // abs
            REQUIRE( abs(v1) == sqrt(2) );
            REQUIRE( abs(yap::FourVector<double>({1, 0, 0, 0})) == 1 );
            REQUIRE( abs(yap::FourVector<double>({1, 1, 0, 0})) == 0 );

        }

        // unit
        //REQUIRE( unit(yap::FourVector<double>({1,0,0,0})) == yap::FourVector<double>({1,0,0,0}) );

        // *

        // -
        SECTION("minus") {
            REQUIRE( -v1 == yap::FourVector<double>({4, -3, -2, -1}) );

            std::vector<yap::FourVector<double>> VV = {v1, v2};
            std::vector<yap::FourVector<double>> mVV = { -v1, -v2};
            REQUIRE( -VV == mVV );
        }

        // cross

        // angle

    }

}
