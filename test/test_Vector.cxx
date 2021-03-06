#include <catch.hpp>

#include <Constants.h>
#include <FourVector.h>
#include <ThreeVector.h>
#include <logging.h>

#include <cmath>

TEST_CASE( "Vector" )
{

    // disable logs in text
    yap::disableLogs(el::Level::Global);

    SECTION( "ThreeVector" ) {

        auto v1 = yap::ThreeVector<double>({1, 2, 3});
        auto v2 = yap::ThreeVector<double>({4, 5, 6});

        SECTION( "addition-assignment" ) {
            v1 += v2;
            REQUIRE( v1 == yap::ThreeVector<double>({5, 7, 9}) );
        }

        SECTION( "subtraction-assignment" ) {
            v1 -= v2;
            REQUIRE( v1 == yap::ThreeVector<double>({ -3, -3, -3}) );
        }

        SECTION( "multiplication-assignment" ) {
            v1 *= 3.;
            REQUIRE( v1 == yap::ThreeVector<double>({3, 6, 9}) );
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

        }

        SECTION( "constants" ) {

            // axes:
            REQUIRE( cross(yap::ThreeAxis_X, yap::ThreeAxis_Y) == yap::ThreeAxis_Z );
            REQUIRE( isRightHanded(yap::ThreeAxes) );
        }

    }

    SECTION( "FourVector" ) {

        auto v1 = yap::FourVector<double>({4, 3, 2, 1});
        auto v2 = yap::FourVector<double>({8, 7, 6, 5});

        SECTION( "addition-assignment" ) {
            v1 += v2;
            REQUIRE( v1 == yap::FourVector<double>({12, 10, 8, 6}) );
        }

        SECTION( "subtraction-assignment" ) {
            v1 -= v2;
            REQUIRE( v1 == yap::FourVector<double>({ -4, -4, -4, -4}) );
        }

        SECTION( "multiplication-assignment" ) {
            v1 *= 2.;
            REQUIRE( v1 == yap::FourVector<double>({8, 6, 4, 2}) );
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

        // cross

        // angle

    }

}
