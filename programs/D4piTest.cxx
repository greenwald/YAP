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
#include "SpinUtilities.h"
#include "WignerD.h"

#include "TGenPhaseSpace.h"
#include "TLorentzVector.h"

#include <assert.h>
#include <iostream>
#include <string>

//#include <callgrind.h>

#include "logging.h"
//INITIALIZE_EASYLOGGINGPP

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
    //yap::ParticleCombination::printParticleCombinationSet();

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

    D->printDecayChain();
    std::cout << "\n";

    D->printSpinAmplitudes();
    D->printDataAccessors(false);
    //D->printDataAccessors();



    // create pseudo data
    TLorentzVector P(0., 0., 0., D->mass()->value());
    Double_t masses[4] = { piPlus->mass()->value(), piMinus->mass()->value(),
                           piPlus->mass()->value(), piMinus->mass()->value()
                         };

    LOG(INFO) << "create dataPoints";

    for (unsigned int iEvt = 0; iEvt < 4; ++iEvt) {
        TGenPhaseSpace event;
        event.SetDecay(P, 4, masses);
        event.Generate();

        std::vector<TLorentzVector> momenta;
        for (unsigned i = 0; i < 4; ++i)
            momenta.push_back(*event.GetDecay(i));

        assert(D->addDataPoint(momenta));
    }

    LOG(INFO) << "done creating dataPoints";

    D->dataSet()[0].printDataSize();

    // create data partitions
    D->setDataPartitions(yap::createDataPartitionsBlock(D->dataSet(), 2));

    // to test amplitude calculation, set all free amps to 1
    auto freeAmps = D->freeAmplitudes();

    LOG(INFO) << freeAmps.size() << " free amplitudes";

    for (auto& a : freeAmps)
        a->setValue(yap::Complex_1);

    //CALLGRIND_START_INSTRUMENTATION

    // do several loops over all dataPartitions
    for (unsigned i = 0; i < 2; ++i) {

        // change amplitudes
        for (auto& a : freeAmps)
            a->setValue(0.9 * a->value());
        DEBUG("===================================================================================================================== ");

        double logA(0);

        if (false) {
            // multi threaded
            logA = D->sumOfLogsOfSquaredAmplitudes();
        } else {
            // update global calculationStatuses before looping over partitions
            D->updateGlobalCalculationStatuses();

            // loop over partitions
            for (yap::DataPartitionBase* partition : D->dataPartitions()) {
                DEBUG("calculate logA for partition " << partition->index() << " ---------------------------------------------------------------------------");
                logA += D->partialSumOfLogsOfSquaredAmplitudes(partition);
            }

            // set parameter flags to unchanged after looping over all partitions
            D->setParameterFlagsToUnchanged();
        }

        LOG(INFO) << "logA = " << logA;
    }

    //CALLGRIND_STOP_INSTRUMENTATION

    /*
        for (auto& a : freeAmps)
            a->setValue(0.5 * a->value());

        std::cout << "try second calculation after changing free amps! ===================================================================================================================== \n";

        D->logLikelihood(d);

        // only change some amps
        unsigned i(0);
        for (auto& a : freeAmps) {
            if (i++ % 2 == 0)
                a->setValue(1.);
        }

        std::cout << "try second calculation after changing free amps! ===================================================================================================================== \n";

        D->logLikelihood(d);

        // set to zero
        for (auto& a : freeAmps) {
            a->setValue(0.);
        }

        std::cout << "try second calculation after changing free amps! ===================================================================================================================== \n";

        D->logLikelihood(d);
    */


    std::cout << "alright! \n";
}
