#ifndef __test_helper_functions__H
#define __test_helper_functions__H

#include <Attributes.h>
#include <BreitWigner.h>
#include <DataSet.h>
#include <FinalStateParticle.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <PHSP.h>
#include <Resonance.h>

/// generate a model with 4 final state particles
//-------------------------
inline std::shared_ptr<yap::Model> d4pi()
{
    auto M = std::make_shared<yap::Model>(std::make_unique<yap::HelicityFormalism>());

    yap::ParticleFactory factory = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");
    double radialSize = 1.;

    // initial state particle
    auto D = factory.decayingParticle(421, radialSize);

    // final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    // Set final-state particles
    M->setFinalState(piPlus, piMinus, piPlus, piMinus);

    // rho
    auto rho = factory.resonance(113, radialSize, std::make_shared<yap::BreitWigner>());
    rho->addStrongDecay(piPlus, piMinus);

    // sigma / f_0(500)
    auto sigma = factory.resonance(9000221, radialSize, std::make_shared<yap::BreitWigner>());
    sigma->addStrongDecay(piPlus, piMinus);

    // a_1
    auto a_1 = factory.resonance(20213, radialSize, std::make_shared<yap::BreitWigner>());
    a_1->addStrongDecay(sigma, piPlus);
    a_1->addStrongDecay(rho,   piPlus);

    // D's channels
    D->addWeakDecay(rho, rho);
    D->addWeakDecay(a_1, piMinus);

    M->addInitialState(D);

    return M;
}

//-------------------------
template <typename Formalism>
inline std::shared_ptr<yap::Model> dkkp(int pdg_D, std::vector<int> fsps)
{
    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    yap::FinalStateParticleVector FSP;
    std::transform(fsps.begin(), fsps.end(), std::back_inserter(FSP), [&F](int pdg){return F.fsp(pdg);});

    double radial_size = 3;

    auto M = std::make_shared<yap::Model>(std::make_unique<Formalism>());
    M->setFinalState(FSP);
    
    auto D = F.decayingParticle(pdg_D, radial_size);

    auto piPlus = lone_elt(filter(FSP, yap::is_named("pi+")));
    auto kPlus  = lone_elt(filter(FSP, yap::is_named("K+")));
    auto kMinus = lone_elt(filter(FSP, yap::is_named("K-")));
    
    for (unsigned j = 0; j < 3; ++j) {
        auto res = yap::Resonance::create("res_" + std::to_string(j), yap::QuantumNumbers(0, j * 2), radial_size,
                                          std::make_shared<yap::BreitWigner>(0.750 + 0.250 * j, 0.025));
        res->addStrongDecay(piPlus, kMinus);
        D->addWeakDecay(res, kPlus);
        switch (j) {
        case 0:
            *free_amplitude(*D, yap::to(res)) = 0.5;
            break;
        case 1:
            *free_amplitude(*D, yap::to(res)) = 1;
            break;
        case 2:
            *free_amplitude(*D, yap::to(res)) = 30;
            break;
        }
    }

    M->addInitialState(D);

    return M;
}

//-------------------------
template <typename Formalism>
std::shared_ptr<yap::Model> d3pi()
{
    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    auto M = std::make_shared<yap::Model>(std::make_unique<Formalism>());

    yap::ParticleFactory factory = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    // initial state particle
    auto D = factory.decayingParticle(factory.pdgCode("D+"), radialSize);

    // final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    // set final state
    M->setFinalState(piPlus, piMinus, piPlus);

    // rho
    auto rho = factory.resonance(113, radialSize, std::make_shared<yap::BreitWigner>());
    rho->addStrongDecay(piPlus, piMinus);

    // Add channels to D
    D->addWeakDecay(rho, piPlus);

    M->addInitialState(D);

    return M;
}

//-------------------------
inline yap::DataSet generate_data(yap::Model& M, unsigned nPoints)
{
    yap::ParticleFactory factory = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");
    
    auto isp_mass = factory[M.initialStates().begin()->first->name()].mass();

    auto A = M.massAxes();
    auto m2r = yap::squared(yap::mass_range(isp_mass, A, M.finalStateParticles()));

    yap::DataSet data(M.createDataSet());

    std::mt19937 g(0);
    // fill data set with nPoints points
    std::generate_n(std::back_inserter(data), nPoints,
                    std::bind(yap::phsp<std::mt19937>, std::cref(M), isp_mass, A, m2r, g, std::numeric_limits<unsigned>::max()));

    return data;
}



#endif

