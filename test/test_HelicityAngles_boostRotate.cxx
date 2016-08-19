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
#include <Rotation.h>

#include <assert.h>
#include <cmath>
#include <random>

/*
 * Test the calculation of helicity angles
 * Downstream theta angles must be the same regardless of rotations/boosts of the final state particles.
 * Phi can change, since it only affects the phase of the amplitude, and it will change in the same way for all amplitudes
 */

TEST_CASE( "HelicityAngles_boostRotate" )
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

    // choose default Dalitz coordinates
    auto A = M.massAxes();
    // get mass^2 ranges
    auto m2r = yap::squared(yap::mass_range(D_mass, A, M.finalStateParticles()));

    // create DataSet
    auto data = M.createDataSet();

    // create random number engine for generation of points
    std::mt19937 g(0);

    // create random number generators
    std::uniform_real_distribution<double> uniform_0_pi(0, yap::pi<double>());
    std::uniform_real_distribution<double> uniform_m99_p99(-0.99, 0.99);

    for (unsigned int iEvt = 0; iEvt < 1000; ++iEvt) {
        std::map<std::shared_ptr<yap::ParticleCombination>, std::vector<double>> resultingThetas;

        // generate random phase space point (with 100 attempts before failing)
        auto momenta = yap::phsp(M, D_mass, A, m2r, g, 100);
        if (momenta.empty())
            continue;

        for (unsigned iTrans = 0; iTrans < 7; ++iTrans) {

            double angle = uniform_0_pi(g);
            double boost = uniform_m99_p99(g);

            for (auto& p : momenta) {

                // testing. Theta of downstream helicity angles must stay the same
                switch (iTrans) {
                    case 1:
                    case 2:
                    case 3:
                        // rotate around axis: case 1, x; case 2, y; case 3, z
                        p = yap::rotation(yap::ThreeAxes[iTrans - 1], angle) * p;
                        break;
                    case 4:
                    case 5:
                    case 6:
                        // boost in direction of axis: case 4, x; case 5, y; case 6, z
                        p = yap::lorentzTransformation(yap::ThreeAxes[iTrans - 3] * boost) * p;
                        break;
                    default:
                        break;
                }
            }

            data.push_back(momenta);
            const auto dp = data.back();

            // compare results
            for (auto& pc_i : M.helicityAngles()->symmetrizationIndices())
                if (pc_i.first->indices().size() < M.finalStateParticles().size())
                    resultingThetas[pc_i.first].push_back(M.helicityAngles()->theta(dp, pc_i.first));
        }


        // check if thetas are equal
        // Phi can change, since it only affects the phase of the amplitude, and it will change in the same way for all amplitudes
        for (auto& kv : resultingThetas) {
            double theta_0 = kv.second[0];
            for (auto theta : kv.second)
                REQUIRE(theta == Approx(theta_0));
        }
    }
}
