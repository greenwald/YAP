// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D4PI__H
#define __BAT__D4PI__H

#include "../bat_fit.h"
#include "../fit_fitFraction.h"
#include "../tools.h"

#include <Attributes.h>
#include <AmplitudeBasis.h>
#include <BreitWigner.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MathUtilities.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <Resonance.h>
#include <SpinAmplitudeCache.h>

#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <complex>
#include <memory>

using namespace yap;

bool a_rho_pi_S  = true;
bool a_rho_pi_D  = true;
bool a_rho_sigma = true;
bool rho_rho     = true;
bool f_0_pipi    = true;
bool f_2_pipi    = true;
bool sigma_pipi  = true;

// scaling to reproduce (approximately) the fit fractions of the FOCUS model
double scale_rho_rho     = 0.7158375483;
double scale_a_rho_pi_D  = 4.2395231085;
double scale_a_rho_sigma = 1.0739853519;
double scale_f_0_pipi    = 0.3762733548;
double scale_f_2_pipi    = 12.687283302;
double scale_sigma_pipi  = 0.204794164 ;


inline std::unique_ptr<Model> d4pi()
{
    auto F = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus = F.fsp(211);
    auto piMinus = F.fsp(-211);

    auto M = std::make_unique<yap::Model>(std::make_unique<yap::HelicityFormalism>());
    M->setFinalState(piPlus, piMinus, piPlus, piMinus);

    // use common radial size for all resonances
    double radialSize = 1.2; // [GeV^-1]

    // initial state particle
    auto D = F.decayingParticle(F.pdgCode("D0"), radialSize);
    
    //
    // resonant particles
    //
    
    // rho
    auto rho = F.resonance(F.pdgCode("rho0"), radialSize, std::make_shared<RelativisticBreitWigner>());
    rho->addStrongDecay(piPlus, piMinus);

    // omega
    //auto omega = F.resonance(F.pdgCode("omega"), radialSize, std::make_shared<BreitWigner>());
    //omega->addStrongDecay({piPlus, piMinus});

    // sigma / f_0(500)
    auto sigma = F.resonance(F.pdgCode("f_0(500)"), radialSize, std::make_shared<BreitWigner>());
    sigma->addStrongDecay(piPlus, piMinus);

    // a_1
    auto a_1 = F.resonance(F.pdgCode("a_1+"), radialSize, std::make_shared<BreitWigner>());
    if (a_rho_pi_S or a_rho_pi_D)
        a_1->addStrongDecay(rho,   piPlus);
    if (a_rho_sigma)
        a_1->addStrongDecay(sigma, piPlus);

    // a_1 -> sigma pi 
    if (a_rho_sigma)
        *free_amplitude(*a_1, to(sigma)) = std::polar(scale_a_rho_sigma * 0.439, rad(193.));

    // f_0(980) (as Flatte)
    auto f_0_980_flatte = std::make_shared<Flatte>();
    f_0_980_flatte->add(FlatteChannel(0.20, *piPlus, *piMinus));
    f_0_980_flatte->add(FlatteChannel(0.50, *F.fsp(321), *F.fsp(-321))); // K+K-
    auto f_0_980 = F.resonance(F.pdgCode("f_0"), radialSize, f_0_980_flatte);
    f_0_980->addStrongDecay(piPlus, piMinus);

    // f_2(1270)
    auto f_2 = F.resonance(F.pdgCode("f_2"), radialSize, std::make_shared<BreitWigner>());
    f_2->addStrongDecay(piPlus, piMinus); 

    // pi+ pi- flat
    auto pipiFlat = DecayingParticle::create("pipiFlat", QuantumNumbers(0, 0), radialSize);
    pipiFlat->addStrongDecay(piPlus, piMinus);
    
    //
    // D0 channels
    //
    if (rho_rho) {
        D->addWeakDecay(rho, rho);

        amplitude_basis::canonical<double> c(amplitude_basis::transversity<double>(
                                                 std::polar(scale_rho_rho * 0.624, rad(357.)),    // A_longitudinal
                                                 std::polar(scale_rho_rho * 0.157, rad(120.)),    // A_parallel
                                                 std::polar(scale_rho_rho * 0.384, rad(163.)) )); // A_perpendicular
        
        for (auto& fa : free_amplitudes(*D, to(rho, rho)))
            *fa = static_cast<std::complex<double> >(c[fa->spinAmplitude()->L()]);
    }
    if (a_rho_pi_S or a_rho_pi_D) {
        D->addWeakDecay(a_1, piMinus);
        free_amplitude(*D, to(a_1))->variableStatus() = VariableStatus::fixed;

        auto a_rho_S = free_amplitude(*a_1, to(rho), l_equals(0));
        auto a_rho_P = free_amplitude(*a_1, to(rho), l_equals(1));
        auto a_rho_D = free_amplitude(*a_1, to(rho), l_equals(2));

        // S wave
        if (a_rho_pi_S) {
            *a_rho_S = 1;// will be fixed in d4pi_fit
            LOG(INFO) << "set a_rho_pi_S to 1";
        }
        else{
            *a_rho_S = 0.;
            a_rho_S->variableStatus() = VariableStatus::fixed;
            LOG(INFO) << "fixed a_rho_S to 0";
        }

        // P wave
        *a_rho_P = 0.;
        a_rho_P->variableStatus() = VariableStatus::fixed;
        LOG(INFO) << "fixed a_rho_P to 0";

        // D wave
        if (a_rho_pi_D)
            *a_rho_D = std::polar(scale_a_rho_pi_D * 0.241, rad(82.));
        else {
            *a_rho_D = 0.;
            a_rho_D->variableStatus() = VariableStatus::fixed;
            LOG(INFO) << "fixed a_rho_D to 0";
        }
    }
    if (f_0_pipi) {
        D->addWeakDecay(f_0_980, piPlus, piMinus);
        *free_amplitude(*D, to(f_0_980, piPlus, piMinus)) = std::polar(scale_f_0_pipi * 0.233, rad(261.));
    }
    if (f_2_pipi) {
        D->addWeakDecay(f_2, pipiFlat);
        *free_amplitude(*D, to(f_2,     pipiFlat       )) = std::polar(scale_f_2_pipi * 0.338, rad(317.));
    }
    if (sigma_pipi) {
        D->addWeakDecay(sigma, piPlus, piMinus);
        *free_amplitude(*D, to(sigma,   piPlus, piMinus)) = std::polar(scale_sigma_pipi * 0.432, rad(254.));
    }
    
    M->addInitialState(D);

    LOG(INFO) << "D Decay trees:";
    LOG(INFO) << to_string(D->decayTrees());

    LOG(INFO) << std::endl << "Free amplitudes: ";
    for (const auto& fa : free_amplitudes(*M, yap::is_not_fixed()))
        LOG(INFO) << yap::to_string(*fa) << "  \t (mag, phase) = (" << abs(fa->value()) << ", " << deg(arg(fa->value())) << "°)";

    return M;
}

inline bat_fit d4pi_fit(std::string name, std::vector<std::vector<unsigned> > pcs = {})
{
    bat_fit m(name, d4pi(), pcs);

    // find particles
    auto D     = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("D0")));
    auto rho   = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("rho0")));
    auto a_1   = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("a_1+")));

    auto fixed_amp = free_amplitude(*a_1, to(rho), l_equals(0));
    if (fixed_amp) {
        m.fix(fixed_amp, abs(fixed_amp->value()), deg(arg(fixed_amp->value())));
        LOG(INFO) << "fixed amplitude " << to_string(*fixed_amp);
    }

    LOG(INFO) << "setting priors";
    //unsigned i = 0;
    for (const auto& fa : m.freeAmplitudes()) {
        double re = real(fa->value());
        double im = imag(fa->value());
        //double ab = abs(fa->value());
        //double ar = deg(arg(fa->value()));
        double rangeLo = 0.5;
        double rangeHi = 1.5;
        //m.setPriors(fa, new ConstantPrior(rangeLo*ab, rangeHi*ab),
        //        new ConstantPrior(ar - (1.-rangeLo) * 360, ar + (rangeHi-1.) * 360));
        m.setRealImagRanges(fa, std::min(rangeLo*re, rangeHi*re), std::max(rangeLo*re, rangeHi*re),
                std::min(rangeLo*im, rangeHi*im), std::max(rangeLo*im, rangeHi*im));
        //m.setAbsArgRanges(fa, rangeLo*ab, rangeHi*ab,
         //       ar - (1.-rangeLo) * 360, ar + (rangeHi-1.) * 360);
/*
        if (++i%3 != 0) {
            m.fix(fa, abs(fa->value()), deg(arg(fa->value())));
            LOG(INFO) << "fixed amplitude " << to_string(*fa);
        }
*/
    }

    return m;
}

inline fit_fitFraction d4pi_fit_fitFraction()
{
    // create bat_fit object
    fit_fitFraction m("D4PI_frac_fit", d4pi());

    //double D_mass = 1.8648400;

    m.GetParameter("N_1").Fix(1);

    // find particles
    auto D     = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("D0")));
    auto rho   = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("rho0")));
    auto sigma = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_0(500)")));
    auto a_1   = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("a_1+")));
    auto f_0   = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_0")));
    auto f_2   = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_2")));

    LOG(INFO) << m.model()->initialStates().at(D).begin()->first;

    // set fit fractions to fit
    /// \todo does not yet work with more than one decayTree
    //m.setFitFraction(decay_trees(*D, yap::from(D), yap::to(a_1), yap::l_equals(0)), 43.3e-2,   quad(2.5e-2, 1.9e-2));
    //m.setFitFraction(decay_trees(*D, yap::from(D), yap::to(a_1), yap::l_equals(1)), 2.5e-2,    quad(0.5e-2, 0.4e-2));
    m.setFitFraction(decay_tree(*D, yap::from(D), yap::to(f_0), yap::l_equals(0)), 8.3e-2,    quad(0.7e-2, 0.6e-2));

    amplitude_basis::canonical<double> c(amplitude_basis::transversity<double>(
            complex_basis::cartesian<double>(std::complex<double>( 1.1e-2),  quad(0.3e-2, 0.3e-2)),
            complex_basis::cartesian<double>(std::complex<double>( 6.4e-2),  quad(0.6e-2, 0.5e-2)),
            complex_basis::cartesian<double>(std::complex<double>(16.88e-2), quad(1.0e-2, 0.8e-2))));

    /// \todo does not yet work with more than one decayTree
    //for (unsigned l = 0; l<3; ++l)
    //    m.setFitFraction(decay_trees(*D, yap::from(D), yap::to(rho), yap::l_equals(l)), real(c.amplitudes()[l]), c.covariance()[l][l][0][0]);

    m.setFitFraction(decay_tree(*D, yap::from(D), yap::to(f_0)),   2.4e-2,  quad(2.4e-2, 0.4e-2));
    m.setFitFraction(decay_tree(*D, yap::from(D), yap::to(f_2)),   4.9e-2,  quad(4.9e-2, 0.5e-2));
    m.setFitFraction(decay_tree(*D, yap::from(D), yap::to(sigma)), 8.2e-2,  quad(8.2e-2, 0.7e-2));

    // set free amplitude parameters of fit
    m.fix(free_amplitude(*D, yap::from(D), yap::to(a_1)), 1., 0.);
    //m.fix(free_amplitude(*D, yap::from(a_1), yap::to(rho), yap::l_equals(1)), 0., 0.);
    m.setPriors(free_amplitude(*D, yap::from(a_1), yap::to(rho), yap::l_equals(2)), new BCGaussianPrior(0.241, quad(0.033, 0.024)), new BCGaussianPrior( 82., quad(5.,   4.)));

    m.setPriors(free_amplitude(*D, yap::from(D), yap::to(f_0), yap::l_equals(0)), new BCGaussianPrior(0.493, quad(0.026, 0.021)), new BCGaussianPrior(193., quad(4.,   4.)));

    // polar -> cartesian; transversity -> canonical
    amplitude_basis::canonical<double> can(amplitude_basis::transversity<double>(
            complex_basis::cartesian<double>(complex_basis::polar<double>(0.624, rad(357.), {quad(0.023, 0.015), quad(3., 3.)})), // A_longitudinal
            complex_basis::cartesian<double>(complex_basis::polar<double>(0.157, rad(120.), {quad(0.027, 0.020), quad(7., 8.)})), // A_parallel
            complex_basis::cartesian<double>(complex_basis::polar<double>(0.384, rad(163.), {quad(0.020, 0.015), quad(3., 3.)})))); // A_perpendicular

    for (unsigned l = 0; l<3; ++l) {
        // cartesian -> polar
        complex_basis::polar<double> polar(can[l]);

        m.setPriors(free_amplitude(*D, yap::from(D), yap::to(rho), yap::l_equals(l)),
                new BCGaussianPrior(polar.value()[0], polar.covariance()[0][0]),
                new BCGaussianPrior(polar.value()[1], polar.covariance()[1][1]));
    }


    m.setPriors(free_amplitude(*D, yap::from(D), yap::to(f_0)),   new BCGaussianPrior(0.233, quad(0.019, 0.015)), new BCGaussianPrior(261., quad(7., 3.)));
    m.setPriors(free_amplitude(*D, yap::from(D), yap::to(f_2)),   new BCGaussianPrior(0.338, quad(0.021, 0.016)), new BCGaussianPrior(317., quad(4., 4.)));
    m.setPriors(free_amplitude(*D, yap::from(D), yap::to(sigma)), new BCGaussianPrior(0.432, quad(0.027, 0.022)), new BCGaussianPrior(254., quad(4., 5.)));

    return m;
}

#endif
