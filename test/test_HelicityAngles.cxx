#include <catch.hpp>
#include <catch_capprox.hpp>

#include <logging.h>
#include <BreitWigner.h>
#include <Constants.h>
#include <DataSet.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <LorentzTransformation.h>
#include <make_unique.h>
#include <MassAxes.h>
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

yap::FourMatrix<double> transformation_to_helicityFrame(const yap::FourVector<double>& daughter)
{
    // rotate to put x-y component of daughter in y direction
    // rotate to put daughter in z direction
    auto R = rotation(yap::ThreeAxis_Y, -yap::theta(vect(daughter), yap::ThreeAxes))
             * rotation(yap::ThreeAxis_Z, -yap::phi(vect(daughter), yap::ThreeAxes));

    // return boost * rotations
    return lorentzTransformation(-(R * daughter)) * lorentzTransformation(R);
}

void calculate_helicity_angles(const yap::Model& M,
                               std::map<const std::shared_ptr<yap::ParticleCombination>, std::array<double, 2> >& phi_theta,
                               const std::shared_ptr<yap::ParticleCombination>& pc,
                               std::vector<yap::FourVector<double> > momenta)
{
    // loop over daughters
    for (const auto& d : pc->daughters()) {

        auto p = yap::FourVector_0;
        // construct 4-vector of daughter
        for (const auto& i : d->indices())
            p += momenta[i];

        // if not yet calculated
        if (phi_theta.find(pc) == phi_theta.end())
            phi_theta[pc] = angles(vect(p), yap::ThreeAxes);

        if (d->daughters().empty())
            continue;

        // next helicity frame
        const auto L = transformation_to_helicityFrame(p);

        for (unsigned i : d->indices())
            momenta[i] = L * momenta[i];

        calculate_helicity_angles(M, phi_theta, d, momenta);
    }
}


TEST_CASE( "HelicityAngles" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    yap::Model M(std::make_unique<yap::HelicityFormalism>());

    yap::ParticleFactory factory = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    double D_mass = factory["D+"].Mass;

    // initial state particle
    auto D = factory.decayingParticle(factory.pdgCode("D+"), radialSize);

    // final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    // set final state
    M.setFinalState({piPlus, piMinus, piPlus});

    // rho
    auto rho = factory.resonance(113, radialSize, std::make_shared<yap::BreitWigner>());
    rho->addChannel({piPlus, piMinus});

    // Add channels to D
    D->addChannel({rho, piPlus});

    M.addInitialStateParticle(D);

    // choose default Dalitz axes
    auto A = M.massAxes();
    // get mass^2 ranges
    auto m2r = yap::squared(yap::mass_range(D_mass, A, M.finalStateParticles()));

    REQUIRE( M.consistent() );

    // create DataSet
    auto data = M.createDataSet();

    // create random number engine for generation of points
    std::mt19937 g(0);

    for (unsigned int iEvt = 0; iEvt < 100; ++iEvt) {

        // generate random phase space point (with 100 attempts before failing)
        auto momenta = yap::phsp(M, D_mass, A, m2r, g, 100);
        if (momenta.empty())
            continue;

        // boost into total rest frame
        momenta = lorentzTransformation(-momenta) * momenta;

        data.push_back(momenta);
        const auto dp = data.back();

        std::map<const std::shared_ptr<yap::ParticleCombination>, std::array<double, 2> > phi_theta;

        for (auto& pc : D->particleCombinations()) {
            REQUIRE( M.fourMomenta()->m(dp, pc) == Approx(D_mass) );

            calculate_helicity_angles(M, phi_theta, pc, momenta); // YAP
        }

        // compare results
        for (auto& kv : phi_theta) {
            REQUIRE( cos(M.helicityAngles()->phi(dp, kv.first))   == Approx(cos(kv.second[0])) );
            REQUIRE( M.helicityAngles()->theta(dp, kv.first) == Approx(kv.second[1]) );
        }
    }
}
