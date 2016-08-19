#ifndef __test_helper_functions__H
#define __test_helper_functions__H

#include "BreitWigner.h"
#include "DataSet.h"
#include "FinalStateParticle.h"
#include "HelicityFormalism.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleFactory.h"
#include "PDL.h"
#include "PHSP.h"
#include "Resonance.h"

/// generate a model with 4 final state particles
//-------------------------
inline std::shared_ptr<yap::Model> create_model()
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
    rho->addChannel(piPlus, piMinus);

    // sigma / f_0(500)
    auto sigma = factory.resonance(9000221, radialSize, std::make_shared<yap::BreitWigner>());
    sigma->addChannel(piPlus, piMinus);

    // a_1
    auto a_1 = factory.resonance(20213, radialSize, std::make_shared<yap::BreitWigner>());
    a_1->addChannel(sigma, piPlus);
    a_1->addChannel(rho,   piPlus);

    // D's channels
    D->addChannel(rho, rho);
    D->addChannel(a_1, piMinus);

    M->addInitialStateParticle(D);

    return M;
}

//-------------------------
inline yap::DataSet generate_data(yap::Model& M, unsigned nPoints)
{
    yap::ParticleFactory factory = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");
    
    auto isp_mass = factory[M.initialStateParticles().begin()->first->name()].Mass;

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

