// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D3PI__H
#define __BAT__D3PI__H

#include "../bat_fit.h"

#include <BreitWigner.h>
#include <Constants.h>
#include <container_utils.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <Filters.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <make_unique.h>
#include <Model.h>
#include <Parameter.h>
#include <Particle.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <PoleMass.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <Resonance.h>
#include <SpinAmplitudeCache.h>

#include <complex>
#include <memory>

using namespace std;
using namespace yap;

inline unique_ptr<Model> d3pi(unique_ptr<SpinAmplitudeCache> SAC)
{
    auto F = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    auto M = make_unique<Model>(move(SAC));
    M->setFinalState(piPlus, piMinus, piPlus);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D+"), radialSize);

    // rho
    /* auto rho = F.resonance(113, radialSize, make_shared<RelativisticBreitWigner>()); */
    auto rho = Resonance::create(F.quantumNumbers(113), "rho", radialSize, make_shared<RelativisticBreitWigner>(775.49e-3, 149.4e-3));
    rho->addChannel(piPlus, piMinus);
    D->addChannel(rho, piPlus);

    // f_2(1270)
    auto f_2 = F.resonance(225, radialSize, make_shared<RelativisticBreitWigner>());
    f_2->addChannel(piPlus, piMinus);
    D->addChannel(f_2, piPlus);

    // f_0(980)
    auto f_0_980_flatte = make_shared<Flatte>(0.965);
    f_0_980_flatte->add(FlatteChannel(0.406, *piPlus, *piMinus));
    f_0_980_flatte->add(FlatteChannel(0.406 * 2, *F.fsp(321), *F.fsp(-321))); // K+K-
    auto f_0_980 = Resonance::create(QuantumNumbers(0, 0), "f_0_980", radialSize, f_0_980_flatte);
    f_0_980->addChannel(piPlus, piMinus);
    D->addChannel(f_0_980, piPlus);

    // f_0(1370)
    auto f_0_1370 = Resonance::create(F["f_0"], "f_0_1370", radialSize, make_unique<RelativisticBreitWigner>(1.350, 0.265));
    f_0_1370->addChannel(piPlus, piMinus);
    D->addChannel(f_0_1370, piPlus);

    // f_0(1500)
    auto f_0_1500 = F.resonance(F.pdgCode("f_0(1500)"), radialSize, make_unique<RelativisticBreitWigner>());
    f_0_1500->addChannel(piPlus, piMinus);
    D->addChannel(f_0_1500, piPlus);

    // sigma a.k.a. f_0(500)
    auto sigma = F.resonance(F.pdgCode("f_0(500)"), radialSize, make_unique<PoleMass>(complex<double>(0.470, -0.220)));
    sigma->addChannel(piPlus, piMinus);
    D->addChannel(sigma, piPlus);

    M->addInitialStateParticle(D);

    // Add channels to D
    *free_amplitude(*M, to(rho))      = polar(1., 0.);
    *free_amplitude(*M, to(f_0_980))  = polar(1.4, rad(12.));
    *free_amplitude(*M, to(f_2))      = polar(2.1, rad(-123.));
    *free_amplitude(*M, to(f_0_1370)) = polar(1.3, rad(-21.));
    *free_amplitude(*M, to(f_0_1500)) = polar(1.1, rad(-44.));
    *free_amplitude(*M, to(sigma))    = polar(3.7, rad(-3.));

    return M;
}

inline bat_fit d3pi_fit(string name, unique_ptr<SpinAmplitudeCache> SAC, vector<vector<unsigned> > pcs = {})
{
    bat_fit m(name, d3pi(move(SAC)), pcs);

    auto rho = dynamic_pointer_cast<Resonance>(particle(*m.model(), is_named("rho")));

    // m.fix(D->freeAmplitudes(rho, piPlus)[0], 1, 0);
    // m.setPrior(D->freeAmplitudes(f_0_980,  piPlus)[0], 0.,  100., -180, 180);
    // m.GetParameter(m.findFreeAmplitude(D->freeAmplitudes(f_0_980, piPlus)[0]) + 1).Fix(12.);
    // m.setPrior(D->freeAmplitudes(f_2,      piPlus)[0], 1.,  3., -130, -110);
    // m.setPrior(D->freeAmplitudes(f_0_1370, piPlus)[0], 1.,  2., -30, -10);
    // m.setPrior(D->freeAmplitudes(f_0_1500, piPlus)[0], 0.5, 2., -50, -30);
    // m.setPrior(D->freeAmplitudes(sigma,    piPlus)[0], 3.,  5., -10, 10);

    m.addParameter("rho_mass", static_pointer_cast<BreitWigner>(rho->massShape())->mass(), 0.5, 1.2);
    m.GetParameters().Back().SetPriorConstant();
    m.addParameter("rho_width", static_pointer_cast<BreitWigner>(rho->massShape())->width(), 0.1, 0.2);
    m.GetParameters().Back().SetPriorConstant();
    // m.addParameter("rho_radialSize", rho->radialSize(), 1, 5);
    // m.GetParameters().Back().SetPriorConstant();
    
    return m;
}

#endif
