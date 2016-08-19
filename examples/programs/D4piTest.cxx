#include "BreitWigner.h"
#include "DataPartition.h"
#include "DataPoint.h"
#include "DataSet.h"
#include "DecayChannel.h"
#include "DecayTree.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
#include "FourVector.h"
#include "FreeAmplitude.h"
#include "HelicityAngles.h"
#include "HelicityFormalism.h"
#include "logging.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "MassShapeWithNominalMass.h"
#include "Model.h"
#include "Parameter.h"
#include "Particle.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "PDL.h"
#include "PHSP.h"
#include "Resonance.h"
#include "SpinAmplitudeCache.h"
#include "WignerD.h"

#include <assert.h>
#include <iostream>
#include <memory>
#include <random>
#include <string>

//#include <callgrind.h>

int main( int argc, char** argv)
{

    //yap::disableLogs(el::Level::Debug);
    yap::plainLogs(el::Level::Debug);

    yap::Model M(std::make_unique<yap::HelicityFormalism>());

    yap::ParticleFactory factory = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    double radialSize = 1.;

    // initial state particle
    auto D = factory.decayingParticle(421, radialSize);

    auto D_mass = factory[421].Mass;

    // final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    // Set final-state particles
    M.setFinalState(piPlus, piMinus, piPlus, piMinus);

    // sigma / f_0(500)
    auto sigma = factory.resonance(9000221, radialSize, std::make_shared<yap::BreitWigner>());
    sigma->addChannel(piPlus, piMinus);

    // rho
    auto rho = factory.resonance(113, radialSize, std::make_shared<yap::BreitWigner>());
    rho->addChannel(piPlus, piMinus);

    // omega
    auto omega = factory.resonance(223, radialSize, std::make_shared<yap::BreitWigner>());
    omega->addChannel(piPlus, piMinus);

    // a_1
    auto a_1 = factory.resonance(20213, radialSize, std::make_shared<yap::BreitWigner>());
    a_1->addChannel(sigma, piPlus);
    a_1->addChannel(rho,   piPlus);

    // D's channels
    D->addChannel(rho, rho);
    D->addChannel(omega, omega);
    D->addChannel(rho, omega);
    D->addChannel(a_1, piMinus);
    D->addChannel(sigma, piPlus, piMinus);
    D->addChannel(piPlus, piMinus, piPlus, piMinus);

    // R pi pi channels
    //yap::Resonance* f_0_980 = factory.resonanceBreitWigner(9000221, radialSize);
    //factory.createChannel(f_0_980, piPlus, piMinus, 0);

    // InitialStateParticles
    M.addInitialStateParticle(D);
    // add other background particles
    M.addInitialStateParticle(a_1);
    M.addInitialStateParticle(rho);

    // check consistency
    if (M.consistent())
        LOG(INFO) << "consistent!";
    else
        LOG(INFO) << "inconsistent!";

    M.lock();

    // print stuff
    //yap::ParticleCombination::printParticleCombinationSet();

    for (auto& isp : M.initialStateParticles()) {
        std::cout << "\n" << isp.first->particleCombinations().size() << " " << *isp.first << " symmetrizations \n";
        for (auto& pc : isp.first->particleCombinations())
            std::cout << *pc << "\n";
        std::cout << "\n";
    }

    std::cout << "\nFour momenta symmetrizations with " << M.fourMomenta()->nSymmetrizationIndices() << " indices \n";
    for (auto& pc_i : M.fourMomenta()->symmetrizationIndices())
        std::cout << *pc_i.first << ": " << pc_i.second << "\n";

    std::cout << "\nHelicity angles symmetrizations with " << M.helicityAngles()->nSymmetrizationIndices() << " indices \n";
    for (auto& pc_i : M.helicityAngles()->symmetrizationIndices())
        std::cout << *pc_i.first << ": " << pc_i.second << "\n";

    D->printDecayChain();
    std::cout << "\n";

    std::cout << *M.spinAmplitudeCache() << std::endl;
    M.printDataAccessors(false);

    LOG(INFO) << "create dataPoints";

    // create default mass axes
    
    auto A = M.massAxes();
    auto m2r = yap::squared(yap::mass_range(D_mass, A, M.finalStateParticles()));

    // create data set
    yap::DataSet data = M.createDataSet();
    // create random number engine for generation of points
    std::mt19937 g(0);
    // fill data set with 1 point
    std::generate_n(std::back_inserter(data), 1,
                    std::bind(yap::phsp<std::mt19937>, std::cref(M), D_mass, A, m2r, g, std::numeric_limits<unsigned>::max()));

    LOG(INFO) << "Size of DataPoint: " + std::to_string(data[0].bytes()) + " byte (for " + std::to_string(data[0].nDataAccessors()) + " data accessors";
    LOG(INFO) << "Printing data:";
    for (unsigned d = 0; d < data.size(); ++d) {
        LOG(INFO) << "  DataPoint " << d;
        for (auto& v : M.fourMomenta()->finalStateMomenta(data[d]))
            LOG(INFO) << yap::to_string(v);
    }

    // create data partitions
    unsigned nChains = 1;
    auto parts = yap::DataPartitionWeave::create(data, nChains);

    //CALLGRIND_START_INSTRUMENTATION

    // create uniform random distributions
    std::uniform_real_distribution<double> uniform;
    std::uniform_real_distribution<double> uniform2(0.95, 1.052631579);

    // do several loops over all dataPartitions
    for (unsigned i = 0; i < 100; ++i) {

        // change amplitudes
        if (uniform(g) > 0.5)
            for (auto& isp_b : M.initialStateParticles())
                for (auto& m_dtv : isp_b.first->decayTrees())
                    for (auto& dt : m_dtv.second)
                        if (dt->freeAmplitude()->variableStatus() != yap::VariableStatus::fixed and uniform(g) > 0.5)
                            *dt->freeAmplitude() = uniform2(g) * dt->freeAmplitude()->value();

        // change masses
        if (uniform(g) > 0.5)
            for (auto& d : particles(M, yap::has_mass))
                if (uniform(g) > 0.5 and mass_parameter(*d).variableStatus() != yap::VariableStatus::fixed) {
                    DEBUG("change mass for " << to_string(*d));
                    mass_parameter(*d) = uniform2(g) * mass_parameter(*d).value();
                }

        DEBUG("===================================================================================================================== ");

        double logA = sum_of_log_intensity(M, parts);
        M.setParameterFlagsToUnchanged();

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


    LOG(INFO) << to_string(D->decayTrees());

    std::cout << "alright! \n";
}
