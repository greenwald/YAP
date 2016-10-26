#include "BreitWigner.h"
#include "DataPartition.h"
#include "DataPoint.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "make_unique.h"
#include "Particle.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "Resonance.h"
#include "WignerD.h"

#include <TGenPhaseSpace.h>
#include <TLorentzVector.h>

#include <assert.h>
#include <iostream>
#include <string>

//#include <callgrind.h>

#include "logging.h"

int main( int argc, char** argv)
{

    //yap::disableLogs(el::Level::Debug);
    yap::plainLogs(el::Level::Debug);

    unsigned max2L(2 * 4);

    yap::ParticleFactory factory((std::string)::getenv("YAPDIR") + "/evt.pdl");

    // initial state particle
    double radialSize = 1.;
    auto D = factory.createInitialStateParticle(421, radialSize);

    // final state particles
    auto piPlus = factory.createFinalStateParticle(211, {0, 2});
    auto piMinus = factory.createFinalStateParticle(-211, {1, 3});

    // rho rho
    auto rho = factory.createResonance(113, radialSize, std::make_unique<yap::BreitWigner>());
    rho->addChannels(piPlus, piMinus, max2L);

    D->addChannels(rho, rho, max2L);

    // omega omega
    auto omega = factory.createResonance(223, radialSize, std::make_unique<yap::BreitWigner>());
    omega->addChannels(piPlus, piMinus, max2L);

    D->addChannels(omega, omega, max2L);

    // rho omega
    D->addChannels(rho, omega, max2L);

    // a_1 channels
    auto sigma = factory.createResonance(9000221, radialSize, std::make_unique<yap::BreitWigner>());
    sigma->addChannels(piPlus, piMinus, max2L);

    auto a_1 = factory.createResonance(20213, radialSize, std::make_unique<yap::BreitWigner>());
    a_1->addChannels(sigma, piPlus, max2L);

    a_1->addChannels(rho, piPlus, max2L);

    D->addChannels(a_1, piMinus, max2L);


    // R pi pi channels
    //yap::Resonance* f_0_980 = factory.createResonanceBreitWigner(9000221, radialSize);
    //factory.createChannel(f_0_980, piPlus, piMinus, 0);


    // consistency and optimizations
    assert(D->prepare());
    std::cout << "consistent! \n";

    // print stuff

    std::cout << "\n" << D->particleCombinations().size() << " D symmetrizations \n";
    /*for (auto& pc : D->particleCombinations())
        std::cout << std::string(*pc) << "\n";
    std::cout << "\n";*/

    std::cout << "\nFour momenta symmetrizations with " << D->fourMomenta().maxSymmetrizationIndex() + 1 << " indices \n";
    /*for (auto& pc : D->fourMomenta().particleCombinations())
        std::cout << std::string(*pc) << ": " << D->fourMomenta().symmetrizationIndex(pc) << "\n";*/

    std::cout << "\nHelicity angles symmetrizations with " << D->helicityAngles().maxSymmetrizationIndex() + 1 << " indices \n";
    /*for (auto& pc : D->helicityAngles().particleCombinations())
        std::cout << std::string(*pc) << ": " << D->helicityAngles().symmetrizationIndex(pc) << "\n";*/

    yap::multiline(LOG(INFO), to_decay_string(*D));

    D->printSpinAmplitudes();
    D->printDataAccessors(false);
    //D->printDataAccessors();



    // create random dataPoint
    TLorentzVector P(0., 0., 0., D->mass()->value());
    Double_t masses[4] = { piPlus->mass()->value(), piMinus->mass()->value(),
                           piPlus->mass()->value(), piMinus->mass()->value()
                         };

    LOG(INFO) << "create dataPoint";

    TGenPhaseSpace event;
    event.SetDecay(P, 4, masses);
    event.Generate();

    std::vector<TLorentzVector> momenta;
    for (unsigned i = 0; i < 4; ++i)
        momenta.push_back(*event.GetDecay(i));

    D->addDataPoint(momenta);
    yap::DataPoint& d = D->dataSet()[0];


    //
    // preparation
    //

    // calculate once
    D->logOfSquaredAmplitude(d, 0);

    // change masses and update calculation Statuses
    std::map<std::shared_ptr<const yap::ParticleCombination>, double> pairMassSquares = D->fourMomenta().pairMassSquares(d);
    pairMassSquares.erase(pairMassSquares.begin());

    for (auto& kv : pairMassSquares) {
        std::cout << std::string(*kv.first) << ": " << kv.second << " GeV\n";
    }

    assert(D->fourMomenta().setMassSquares(d, pairMassSquares));
    D->updateGlobalCalculationStatuses();


    D->fourMomenta().printMasses(d);



    // Dalitz coordinates: m^2_{12}, m^2_{13}, m^2_{14}, m^2_{23}, m^2_{24}

    //
    // generate
    //

    for (unsigned i = 0; i < 10; ++i) {

        pairMassSquares.begin()->second += 0.001;

        for (auto& kv : pairMassSquares) {
            std::cout << std::string(*kv.first) << ": \t" << kv.second << " GeV^2\n";
        }

        if (! D->fourMomenta().setMassSquares(d, pairMassSquares)) {
            continue;
        }

        std::cout << "meh nr " << i << "\n";
        D->fourMomenta().printMasses(d);

        D->calculate(d);

        std::cout << "meh nr " << i << " after calculate\n";
        D->fourMomenta().printMasses(d);

        double logA = D->logOfSquaredAmplitude(d, 0u);

        std::cout << logA << " \n";
    }



    std::cout << "alright! \n";
}
