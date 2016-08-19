#include "logging.h"
#include "BreitWigner.h"
#include "DataSet.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
#include "HelicityAngles.h"
#include "HelicityFormalism.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "PDL.h"
#include "Resonance.h"
#include "SpinAmplitudeCache.h"

#include <memory>
#include <vector>

int main( int argc, char** argv)
{

    yap::plainLogs(el::Level::Debug);

    yap::Model M(std::make_unique<yap::HelicityFormalism>());

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    yap::ParticleFactory factory = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    auto D_mass = factory["D+"].Mass;

    // create final state particles
    auto kPlus  = factory.fsp(+321);
    auto kMinus = factory.fsp(-321);
    auto piPlus = factory.fsp(+211);

    M.setFinalState(kPlus, kMinus, piPlus);

    // create initial state particle and set final state
    auto D = factory.decayingParticle(factory.pdgCode("D+"), radialSize);

    // create a phi
    auto phi = factory.resonance(factory["phi"].PDG, radialSize, std::make_shared<yap::BreitWigner>());
    phi->addChannel(kPlus, kMinus);

    // Add channels to D
    D->addChannel(phi, piPlus);

    M.addInitialStateParticle(D);

    // check consistency
    if (M.consistent())
        LOG(INFO) << "consistent!";
    else
        LOG(INFO) << "inconsistent!";

    // print stuff
    //yap::ParticleCombination::printParticleCombinationSet();

    std::cout << "\n" << D->particleCombinations().size() << " D symmetrizations \n";

    std::cout << "\nFour momenta symmetrizations with " << M.fourMomenta()->nSymmetrizationIndices() << " indices \n";

    std::cout << "\nHelicity angle symmetrizations with " << M.helicityAngles()->nSymmetrizationIndices() << " indices \n";

    D->printDecayChain();
    std::cout << "\n";

    std::cout << *M.spinAmplitudeCache() << std::endl;
    M.printDataAccessors(false);

    // choose default Dalitz coordinates
    const yap::MassAxes massAxes = M.massAxes();

    std::vector<double> m2(massAxes.size(), 1);

    // create data set with one empty data point
    auto data = M.createDataSet(1);

    DEBUG("BEFORE");
    M.fourMomenta()->printMasses(data[0]);

    LOG(INFO) << "setting squared mass ...";
    auto P = calculate_four_momenta(D_mass, M, massAxes, m2);
    if (P.empty())
        LOG(INFO) << "... outside phase space";
    else {
        LOG(INFO) << "... inside phase space";
        M.setFinalStateMomenta(data[0], P, data);
    }

    DEBUG("AFTER");
    M.fourMomenta()->printMasses(data[0]);

    std::cout << "alright! \n";
}
