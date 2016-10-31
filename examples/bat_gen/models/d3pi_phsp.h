// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D3PI_PHSP__H
#define __BAT__D3PI_PHSP__H

#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <Model.h>
#include <ParticleFactory.h>
#include <PDL.h>

#include <memory>

inline std::unique_ptr<yap::Model> d3pi_phsp(std::unique_ptr<yap::Model> M)
{
    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    M->setFinalState(piPlus, piMinus, piPlus);

    auto D = F.decayingParticle(F.pdgCode("D+"), 3 /*Gev^-1*/);
    D->addWeakDecay(piPlus, piMinus, piPlus);

    M->addInitialState(D);

    return M;
}

#endif
