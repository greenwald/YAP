// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D3PI__H
#define __BAT__D3PI__H

#include "B_particle_list.h"

#include <Attributes.h>
#include <BreitWigner.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <FreeAmplitude.h>
#include <make_unique.h>
#include <MathUtilities.h>
#include <Model.h>
#include <Parameter.h>
#include <Particle.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <QuantumNumbers.h>
#include <SpinAmplitudeCache.h>

#include <complex>
#include <memory>

using namespace std;
using namespace yap;

inline unique_ptr<Model> B0_Dpipi0(unique_ptr<Model> M)
{
    auto T = B_particle_list();

    // use common radial size for all resonances
    double r = 3.; // [GeV^-1]

    // final state
    auto Dminus = FinalStateParticle::create(T["D-"]);
    auto piplus = FinalStateParticle::create(T["pi+"]);
    auto pizero = FinalStateParticle::create(T["pi0"]);

    M->setFinalState(Dminus, piplus, pizero);

    // initial state particle
    auto B0 = DecayingParticle::create(T["B0"], r);

    auto rho770  = DecayingParticle::create(T["rho(770)+"],  r, make_shared<BreitWigner>(T["rho(770)+"]));
    auto rho1450 = DecayingParticle::create(T["rho(1450)+"], r, make_shared<BreitWigner>(T["rho(1450)+"]));
    auto rho1700 = DecayingParticle::create(T["rho(1700)+"], r, make_shared<BreitWigner>(T["rho(1700)+"]));
    
    auto D2400z = DecayingParticle::create(T["D_0*(2400)0"], r, make_shared<BreitWigner>(T["D_0*(2400)0"]));
    auto D2420z = DecayingParticle::create(T["D_1(2420)0"],  r, make_shared<BreitWigner>(T["D_1(2420)0"]));
    auto D2460z = DecayingParticle::create(T["D_2*(2460)0"], r, make_shared<BreitWigner>(T["D_2*(2460)0"]));

    auto D2400m = DecayingParticle::create(T["D_0*(2400)-"], r, make_shared<BreitWigner>(T["D_0*(2400)-"]));
    auto D2420m = DecayingParticle::create(T["D_1(2420)-"],  r, make_shared<BreitWigner>(T["D_1(2420)-"]));
    auto D2460m = DecayingParticle::create(T["D_2*(2460)-"], r, make_shared<BreitWigner>(T["D_2*(2460)-"]));

    for (auto& rho : {rho770/*, rho1450, rho1700*/}) {
        rho->addStrongDecay(piplus, pizero);
        B0->addWeakDecay(Dminus, rho);
    }

    for (auto& Dstar0 : {D2400z, /*D2420z,*/ D2460z}) {
        Dstar0->addStrongDecay(Dminus, piplus);
        B0->addWeakDecay(Dstar0, pizero);
    }

    for (auto& Dstarminus : {D2400m, /*D2420m,*/ D2460m}) {
        Dstarminus->addStrongDecay(Dminus, pizero);
        B0->addWeakDecay(Dstarminus, piplus);
    }

    // add nonresonant decay
    B0->addWeakDecay(Dminus, piplus, pizero);
    
    *free_amplitude(*B0, to(rho770))  = polar(1., 0.);
    // *free_amplitude(*B0, yap::to(rho1450)) = polar(1.4, yap::rad(12.));
    // *free_amplitude(*B0, yap::to(rho1700)) = polar(1.4, yap::rad(12.));
    *free_amplitude(*B0, to(D2400z)) = complex<double>(-0.194, -0.152);
    *free_amplitude(*B0, to(D2460z)) = complex<double>(0.048, -0.004);
    *free_amplitude(*B0, to(D2400m)) = complex<double>(0.042, 0.018);
    *free_amplitude(*B0, to(D2460m)) = complex<double>(0.123, 0.028);
    *free_amplitude(*B0, to(Dminus, piplus, pizero)) = complex<double>(-0.045, 0);
    
    return M;
}

#endif
