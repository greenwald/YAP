// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__DKKPI__H
#define __BAT__DKKPI__H

#include "BreitWigner.h"
#include "Constants.h"
#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "make_unique.h"
#include <Model.h>
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "PDL.h"
#include "QuantumNumbers.h"
#include "RelativisticBreitWigner.h"
#include "Resonance.h"
#include <SpinAmplitudeCache.h>

#include <complex>
#include <memory>

using namespace std;
using namespace yap;

inline unique_ptr<Model> dkkpi(unique_ptr<SpinAmplitudeCache> SAC)
{
    auto F = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto kPlus  = F.fsp(+321);
    auto kMinus = F.fsp(-321);
    auto piPlus = F.fsp(+211);

    auto M = make_unique<Model>(move(SAC));

    M->setFinalState(kPlus, kMinus, piPlus);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);

    auto KK0 = Resonance::create(QuantumNumbers(0, 0), "KK0", radialSize, make_shared<BreitWigner>(1.1, 0.075));
    KK0->addChannel(kPlus, kMinus);
    D->addChannel(KK0, piPlus);
    *free_amplitude(*D, to(KK0)) = polar(0.15, rad(39.));
    
    auto KK1 = Resonance::create(QuantumNumbers(2, 0), "KK1", radialSize, make_shared<BreitWigner>(1.35, 0.125));
    KK1->addChannel(kPlus, kMinus);
    D->addChannel(KK1, piPlus);
    *free_amplitude(*D, to(KK1)) = 1.;
    
    auto KK2 = Resonance::create(QuantumNumbers(4, 0), "KK2", radialSize, make_shared<BreitWigner>(1.6, 0.100));
    KK2->addChannel(kPlus, kMinus);
    D->addChannel(KK2, piPlus);
    *free_amplitude(*D, to(KK2)) = polar(10., rad(-12.));
    
    auto piK0 = Resonance::create(QuantumNumbers(0, 0), "piK0", radialSize, make_shared<BreitWigner>(0.75, 0.085));
    piK0->addChannel(piPlus, kMinus);
    D->addChannel(piK0, kPlus);
    *free_amplitude(*D, to(piK0)) = polar(0.23, rad(112.));

    auto piK1 = Resonance::create(QuantumNumbers(2, 0), "piK1", radialSize, make_shared<BreitWigner>(1.0, 0.125));
    piK1->addChannel(piPlus, kMinus);
    D->addChannel(piK1, kPlus);
    *free_amplitude(*D, to(piK1)) = polar(1.2, rad(-76.));

    auto piK2 = Resonance::create(QuantumNumbers(4, 0), "piK2", radialSize, make_shared<BreitWigner>(1.25, 0.065));
    piK2->addChannel(piPlus, kMinus);
    D->addChannel(piK2, kPlus);
    *free_amplitude(*D, to(piK2)) = polar(7.8, rad(56.));
    
    M->addInitialStateParticle(D);

    return M;
}

inline bat_fit dkkpi_fit(string name, unique_ptr<SpinAmplitudeCache> SAC, vector<vector<unsigned> > pcs = {})
{
    bat_fit m(name, dkkpi(move(SAC)), pcs);
    
    auto KK0 = dynamic_pointer_cast<Resonance>(particle(*m.model(), is_named("KK0")));
    auto KK1 = dynamic_pointer_cast<Resonance>(particle(*m.model(), is_named("KK1")));
    auto KK2 = dynamic_pointer_cast<Resonance>(particle(*m.model(), is_named("KK2")));

    auto piK0 = dynamic_pointer_cast<Resonance>(particle(*m.model(), is_named("piK0")));
    auto piK1 = dynamic_pointer_cast<Resonance>(particle(*m.model(), is_named("piK1")));
    auto piK2 = dynamic_pointer_cast<Resonance>(particle(*m.model(), is_named("piK2")));

    m.setPrior(free_amplitude(*m.model(), to(KK0)), 0.05, 0.25, -160., 160.);
    m.fix(free_amplitude(*m.model(), to(KK1)), 1., 0.);
    m.setPrior(free_amplitude(*m.model(), to(KK2)), 1., 25., -160., 160.);

    m.setPrior(free_amplitude(*m.model(), to(piK0)), 0.05, 0.45, -160., 160.);
    m.setPrior(free_amplitude(*m.model(), to(piK1)), 0.2, 3., -160., 160.);
    m.setPrior(free_amplitude(*m.model(), to(piK2)), 1., 25., -160., 160.);

    // fix some phases
    m.GetParameter(m.findFreeAmplitude(free_amplitude(*m.model(), to(piK0))) + 1).Fix(112.);
    m.GetParameter(m.findFreeAmplitude(free_amplitude(*m.model(), to(piK1))) + 1).Fix(-76.);
    m.GetParameter(m.findFreeAmplitude(free_amplitude(*m.model(), to(piK2))) + 1).Fix(56.);

    /* m.addParameter("KK0_mass", KK0->mass(), 1.1 - 2e-3, 1.1 + 2e-3); */
    /* m.GetParameters().Back().SetPriorConstant(); */
    
    /* m.addParameter("KK0_width", static_pointer_cast<RelativisticBreitWigner>(KK0->massShape())->width(), 25e-3 - 5e-3, 25e-3 + 5e-3); */
    /* m.GetParameters().Back().SetPriorConstant(); */

    /* m.addParameter("KK1_mass", KK1->mass(), 1.35 - 10e-3, 1.35 + 10e-3); */
    /* m.GetParameters().Back().SetPriorConstant(); */
    
    /* m.addParameter("KK1_width", static_pointer_cast<BreitWigner>(KK1->massShape())->width(), 50e-3 - 5e-3, 50e-3 + 5e-3); */
    /* m.GetParameters().Back().SetPriorConstant(); */

    return m;
}

#endif
