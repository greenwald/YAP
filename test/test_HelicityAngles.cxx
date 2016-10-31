#include <catch.hpp>
#include <catch_capprox.hpp>

#include "helperFunctions.h"

#include <logging.h>
#include <BreitWigner.h>
#include <DataSet.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <LorentzTransformation.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <MathUtilities.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <PHSP.h>
#include <Resonance.h>

#include <assert.h>
#include <cmath>
#include <random>

/*
 * Test the calculation of helicity angles (with passive rotation in YAP)
 * against rootPWA style (with active rotations of four momenta)
 */

yap::CoordinateSystem<double, 3> three_axes({yap::ThreeVector<double>({1., 0., 0.}),
                                             yap::ThreeVector<double>({0., 1., 0.}),
                                             yap::ThreeVector<double>({0., 0., 1.})});

yap::FourMatrix<double> transformation_to_helicityFrame(const yap::FourVector<double>& daughter)
{
    // rotate to put x-y component of daughter in y direction
    // rotate to put daughter in z direction
    auto R = rotation(three_axes[1], -yap::theta(vect(daughter), three_axes))
             * rotation(three_axes[2], -yap::phi(vect(daughter), three_axes));

    // return boost * rotations
    return lorentzTransformation(-(R * daughter)) * lorentzTransformation(R);
}

void calculate_helicity_angles(const yap::Model& M,
                               yap::ParticleCombinationMap<yap::spherical_angles<double> >& hel_angles,
                               const std::shared_ptr<const yap::ParticleCombination>& pc,
                               std::vector<yap::FourVector<double> > momenta)
{
    // loop over daughters
    for (const auto& d : pc->daughters()) {

        yap::FourVector<double> p;
        // construct 4-vector of daughter
        for (const auto& i : d->indices())
            p += momenta[i];

        // if not yet calculated
        if (hel_angles.find(pc) == hel_angles.end())
            hel_angles[pc] = angles(vect(p), three_axes);

        if (d->daughters().empty())
            continue;

        // next helicity frame
        const auto L = transformation_to_helicityFrame(p);

        for (unsigned i : d->indices())
            momenta[i] = L * momenta[i];

        calculate_helicity_angles(M, hel_angles, d, momenta);
    }
}


TEST_CASE( "HelicityAngles" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    yap::ParticleFactory factory = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    double D_mass = factory["D+"].mass();

    auto M = d3pi<yap::HelicityFormalism>();

    // choose default Dalitz axes
    auto A = M->massAxes();
    // get mass^2 ranges
    auto m2r = yap::squared(yap::mass_range(D_mass, A, M->finalStateParticles()));

    REQUIRE( M->consistent() );

    // create DataSet
    auto data = M->createDataSet();

    // create random number engine for generation of points
    std::mt19937 g(0);

    for (unsigned int iEvt = 0; iEvt < 100; ++iEvt) {

        // generate random phase space point (with 100 attempts before failing)
        auto momenta = yap::phsp(*M, D_mass, A, m2r, g, 100);
        if (momenta.empty())
            continue;

        // boost into total rest frame
        momenta = lorentzTransformation(-momenta) * momenta;

        data.push_back(momenta);
        const auto dp = data.back();

        yap::ParticleCombinationMap<yap::spherical_angles<double> > hel_angles;

        for (auto& pc : M->initialStates().begin()->first->particleCombinations()) {
            REQUIRE( M->fourMomenta()->m(dp, pc) == Approx(D_mass) );

            calculate_helicity_angles(*M, hel_angles, pc, momenta); // YAP
        }

        // compare results
        for (auto& kv : hel_angles) {
            REQUIRE( cos(M->helicityAngles()(dp, data, kv.first).phi)   == Approx(cos(kv.second.phi)) );
            REQUIRE( M->helicityAngles()(dp, data, kv.first).theta == Approx(kv.second.theta) );
        }
    }
}
