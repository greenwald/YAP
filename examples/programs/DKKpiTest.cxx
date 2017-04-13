#include <logging.h>
#include <BreitWigner.h>
#include <DataSet.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <MassAxes.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <SpinAmplitudeCache.h>

#include <memory>
#include <vector>

int main( int argc, char** argv)
{

    yap::plainLogs(el::Level::Debug);

    yap::Model M(std::make_unique<yap::HelicityFormalism>());

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    auto T = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    auto D_mass = T["D+"].mass();

    // create final state particles
    auto kPlus  = yap::FinalStateParticle::create(T[+321]);
    auto kMinus = yap::FinalStateParticle::create(T[-321]);
    auto piPlus = yap::FinalStateParticle::create(T[+211]);

    M.setFinalState(kPlus, kMinus, piPlus);

    // create initial state particle and set final state
    auto D = yap::DecayingParticle::create(T["D+"], radialSize);

    // create a phi
    auto phi = yap::DecayingParticle::create(T["phi"], radialSize, std::make_shared<yap::BreitWigner>(T["phi"]));
    phi->addStrongDecay(kPlus, kMinus);

    // Add channels to D
    D->addWeakDecay(phi, piPlus);

    // check consistency
    if (M.consistent())
        LOG(INFO) << "consistent!";
    else
        LOG(INFO) << "inconsistent!";

    // print stuff

    FLOG(INFO) << "";
    FLOG(INFO) << D->particleCombinations().size() << " D symmetrizations";

    FLOG(INFO) << "";
    FLOG(INFO) << "Four momenta symmetrizations with " << M.fourMomenta()->nSymmetrizationIndices() << " indices";

    MULTILINE(FLOG(INFO),to_decay_string(*D));
    FLOG(INFO) << "";

    MULTILINE(FLOG(INFO),to_string(*M.spinAmplitudeCache()));

    // choose default Dalitz coordinates
    const yap::MassAxes massAxes = M.massAxes();

    std::vector<double> m2(massAxes.size(), 1);

    // create data set with one empty data point
    auto data = M.createDataSet(1);

    DEBUG("BEFORE");
    MULTILINE(FLOG(INFO),M.fourMomenta()->massesString(data[0]));

    LOG(INFO) << "setting squared mass ...";
    auto P = calculate_four_momenta(D_mass, M, massAxes, m2);
    if (P.empty())
        LOG(INFO) << "... outside phase space";
    else {
        LOG(INFO) << "... inside phase space";
        M.setFinalStateMomenta(data[0], P, data);
    }

    DEBUG("AFTER");
    MULTILINE(FLOG(INFO),M.fourMomenta()->massesString(data[0]));

    FLOG(INFO) << "alright! \n";
}
