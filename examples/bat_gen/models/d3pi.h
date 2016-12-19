// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D3PI__H
#define __BAT__D3PI__H

#include "../bat_fit.h"
#include "../ConstantPrior.h"
#include "../fit_fitFraction.h"
#include "../tools.h"

#include <Attributes.h>
#include <BreitWigner.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FreeAmplitude.h>
#include <make_unique.h>
#include <MathUtilities.h>
#include <Model.h>
#include <Parameter.h>
#include <Particle.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <PoleMass.h>
#include <QuantumNumbers.h>
#include <RelativisticBreitWigner.h>
#include <SpinAmplitudeCache.h>
#include <StepFunction.h>

#include <BAT/BCGaussianPrior.h>
#include <BAT/BCSplitGaussianPrior.h>

#include <complex>
#include <memory>

using namespace std;
using namespace yap;

inline unique_ptr<Model> d3pi(unique_ptr<Model> M)
{
    auto T = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus  = FinalStateParticle::create(T[211]);
    auto piMinus = FinalStateParticle::create(T[-211]);

    M->setFinalState(piPlus, piMinus, piPlus);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = DecayingParticle::create(T["D+"], radialSize);

    // rho
    auto rho = DecayingParticle::create(T[113], radialSize, make_shared<RelativisticBreitWigner>(775.49e-3, 149.4e-3));
    rho->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(rho, piPlus);

    // f_2(1270)
    auto f_2 = DecayingParticle::create(T[225], radialSize, make_shared<RelativisticBreitWigner>(T[225]));
    f_2->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(f_2, piPlus);

    // f_0(980)
    auto f_0_980_flatte = make_shared<Flatte>(0.965);
    f_0_980_flatte->add(FlatteChannel(0.406 * 0.406, *piPlus, *piMinus));
    f_0_980_flatte->add(FlatteChannel(0.406 * 0.406 * 4, T[321], T[-321])); // K+K-
    auto f_0_980 = DecayingParticle::create("f_0(980)", T["f_0"].quantumNumbers(), radialSize, f_0_980_flatte);
    f_0_980->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(f_0_980, piPlus);

    // f_0(1370)
    auto f_0_1370 = DecayingParticle::create("f_0(1370)", T["f_0"].quantumNumbers(), radialSize, make_unique<RelativisticBreitWigner>(1.350, 0.265));
    f_0_1370->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(f_0_1370, piPlus);

    // f_0(1500)
    auto f_0_1500 = DecayingParticle::create(T["f_0(1500)"], radialSize, make_unique<RelativisticBreitWigner>(T["f_0(1500)"]));
    f_0_1500->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(f_0_1500, piPlus);

    // sigma a.k.a. f_0(500)
    auto sigma = DecayingParticle::create(T["f_0(500)"], radialSize, make_unique<PoleMass>(complex<double>(0.470, -0.220)));
    sigma->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(sigma, piPlus);

    // Add channels to D
    *free_amplitude(*M, to(rho))      = polar(1., 0.);
    *free_amplitude(*M, to(f_0_980))  = polar(1.4, rad(12.));
    *free_amplitude(*M, to(f_2))      = polar(2.1, rad(-123.));
    *free_amplitude(*M, to(f_0_1370)) = polar(1.3, rad(-21.));
    *free_amplitude(*M, to(f_0_1500)) = polar(1.1, rad(-44.));
    *free_amplitude(*M, to(sigma))    = polar(3.7, rad(-3.));

    /* *free_amplitude(*M, to(f_0_fake)) = polar(1., 60.); */

    return M;
}

inline unique_ptr<Model> d3pi_mi(unique_ptr<Model> M)
{
    auto T = read_pdl_file((string)::getenv("YAPDIR") + "/data/evt.pdl");

    // final state particles
    auto piPlus  = FinalStateParticle::create(T[211]);
    auto piMinus = FinalStateParticle::create(T[-211]);

    M->setFinalState(piPlus, piMinus, piPlus);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    // initial state particle
    auto D = DecayingParticle::create(T["D+"], radialSize);

    // rho
    auto rho = DecayingParticle::create(T[113], radialSize, make_shared<RelativisticBreitWigner>(775.49e-3, 149.4e-3));
    rho->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(rho, piPlus);

    // f_2(1270)
    auto f_2 = DecayingParticle::create(T[225], radialSize, make_shared<RelativisticBreitWigner>(T[225]));
    f_2->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(f_2, piPlus);
    
    // model independent pipi S wave
    std::vector<double> bins;
    double m2_low = pow(2 * T[211].mass(), 2);
    double m2_high = pow(T["D+"].mass() - T[211].mass(), 2);
    double m2_f0 = pow(T["f_0(1500)"].mass(), 2);
    double mw_f0 = T["f_0(1500)"].mass() * T["f_0(1500)"].massShapeParameters()[0];
    double phi_min = atan((mw_f0 / (m2_f0 - m2_low) - (m2_f0 - m2_low) / mw_f0) / 2);
    double phi_max = atan((mw_f0 / (m2_f0 - m2_high) - (m2_f0 - m2_high) / mw_f0) / 2) + yap::pi();
    unsigned Nbins = 20;
    for (unsigned b = 0; b <= Nbins; ++b) {
        auto amp = std::polar(0.5 / mw_f0, phi_min + (phi_max - phi_min) * b / 20) + std::complex<double>(0, 0.5 / mw_f0);
        bins.push_back(m2_f0 - mw_f0 / tan(arg(amp)));
    }
    auto pipi_Swave = DecayingParticle::create("pipi_Swave", T["f_0(500)"].quantumNumbers(), radialSize, make_shared<StepFunction>(bins));
    pipi_Swave->addStrongDecay(piPlus, piMinus);
    D->addWeakDecay(pipi_Swave, piPlus);

	M->lock();

    // Add channels to D
    *free_amplitude(*M, to(rho))        = polar(1., 0.);
    *free_amplitude(*M, to(f_2))        = polar(2.1, rad(-123.));
    *free_amplitude(*M, to(pipi_Swave)) = polar(3., rad(-44.));

    auto step_function = static_pointer_cast<StepFunction>(pipi_Swave->massShape());
    double m_f0 = T["f_0(1500)"].mass();
    double w_f0 = T["f_0(1500)"].massShapeParameters()[0];
    for (size_t i = 0; i < step_function->steps().size(); ++i) {
        double m2 = 0.5 * (step_function->lowEdges()[i]->value() + step_function->lowEdges()[i + 1]->value());
        *step_function->freeAmplitudes()[i] = 1.1 / (m_f0 * m_f0 - m2 - 1_i * m_f0 * w_f0);
    }
    
    return M;
}

inline bat_fit d3pi_fit(string name, unique_ptr<Model> M, vector<vector<unsigned> > pcs = {})
{
    bat_fit m(name, d3pi_mi(move(M)), pcs);

    auto is_to_rho = to(std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("rho0"))));
    auto is_to_pipiS = to(std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("pipi_Swave"))));
    auto is_from_pipiS = from(std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("pipi_Swave"))));
    
    for (const auto& fa : free_amplitudes(*m.model(), is_not_fixed())) {
        m.setPriors(fa, new ConstantPrior(0, 5), new ConstantPrior(-180, 180));
        m.setRealImagRanges(fa, -2, 2, -2, 2);
        m.setAbsArgRanges(fa, 0, 2, -180, 180);
        /* if (!is_from_pipiS(fa)) */
        if (is_to_rho(fa) or is_to_pipiS(fa))
            m.fix(fa, real(fa->value()), imag(fa->value()));
    }

    /* auto f_0_1370 = dynamic_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(1370)"))); */
    /* m.addParameter("mass(f_0(1370))", static_pointer_cast<BreitWigner>(f_0_1370->massShape())->mass(), 1.2, 1.5); */
    /* m.GetParameters().Back().SetPriorConstant(); */
    /* m.addParameter("width(f_0(1370))", static_pointer_cast<BreitWigner>(f_0_1370->massShape())->width(), .2, .5); */
    /* m.GetParameters().Back().SetPriorConstant(); */

    /* auto f_0_1500 = dynamic_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(1500)"))); */
    /* m.addParameter("mass(f_0(1500))", static_pointer_cast<BreitWigner>(f_0_1500->massShape())->mass(), 1.475, 1.535); */
    /* m.GetParameters().Back().SetPrior(new BCGaussianPrior(1.505, 6e-3)); */
    /* m.addParameter("width(f_0(1500))", static_pointer_cast<BreitWigner>(f_0_1500->massShape())->width(), .075, .143); */
    /* m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.109, 7e-3)); */

    /* auto f_2 = dynamic_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_2"))); */
    /* m.addParameter("mass(f_2)", static_pointer_cast<BreitWigner>(f_2->massShape())->mass(), 1.270, 1.280); */
    /* m.GetParameters().Back().SetPrior(new BCGaussianPrior(1.2751, 1.2e-3)); */
    /* m.addParameter("width(f_2)", static_pointer_cast<BreitWigner>(f_2->massShape())->width(), .170, .200); */
    /* m.GetParameters().Back().SetPrior(new BCSplitGaussianPrior(0.185, 2.4e-3, 2.9e-3)); */

    /* auto rho = dynamic_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("rho0"))); */
    /* m.addParameter("rho_mass", static_pointer_cast<BreitWigner>(rho->massShape())->mass(), 0.5, 1.2); */
    /* m.GetParameters().Back().SetPriorConstant(); */
    /* m.addParameter("rho_width", static_pointer_cast<BreitWigner>(rho->massShape())->width(), 0.1, 0.2); */
    /* m.GetParameters().Back().SetPriorConstant(); */
    // m.addParameter("rho_radialSize", rho->radialSize(), 1, 5);
    // m.GetParameters().Back().SetPriorConstant();
    
    return m;
}

inline fit_fitFraction d3pi_fit_fitFractions(string name, unique_ptr<Model> M, vector<vector<unsigned> > pcs = {})
{
    fit_fitFraction m(name, d3pi(move(M)), pcs);
    
    // find particles
    auto D        = static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("D+")));
    auto rho      = static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("rho0")));
    auto f_2      = static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_2")));
    auto f_0_980  = static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(980)")));
    auto f_0_1370 = static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(1370)")));
    auto f_0_1500 = static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(1500)")));
    auto sigma    = static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(500)")));

    // set fit fractions to fit
    m.setFitFraction(decay_tree(*D, to(rho)),      20e-2,   quad(2.3e-2, 0.9e-2));
    m.setFitFraction(decay_tree(*D, to(f_2)),      18.2e-2, quad(2.6e-2, 0.7e-2));
    m.setFitFraction(decay_tree(*D, to(f_0_980)),  4.1e-2,  quad(0.9e-2, 0.3e-2));
    m.setFitFraction(decay_tree(*D, to(f_0_1370)), 2.6e-2,  quad(1.8e-2, 0.6e-2));
    m.setFitFraction(decay_tree(*D, to(f_0_1500)), 3.4e-2,  quad(1.0e-2, 0.8e-2));
    m.setFitFraction(decay_tree(*D, to(sigma)),    41.8e-2, quad(1.4e-2, 2.5e-2));

    // set free amplitude parameters of fit
    m.fix(free_amplitude(*m.model(), to(rho)),      1.,  0.);
    m.fix(free_amplitude(*m.model(), to(f_2)),      2.1, -123.);
    m.fix(free_amplitude(*m.model(), to(f_0_980)),  1.4, 12.);
    m.fix(free_amplitude(*m.model(), to(f_0_1370)), 1.3, -21.);
    m.fix(free_amplitude(*m.model(), to(f_0_1500)), 1.1, -44.);
    m.fix(free_amplitude(*m.model(), to(sigma)),    3.7, -3);
    /* m.setPrior(free_amplitude(*m.model(), to(f_2)),      new BCGaussianPrior(2.1, quad(0.2, 0.1)), new BCGaussianPrior(-123., quad(6.,   3.))); */
    /* m.setPrior(free_amplitude(*m.model(), to(f_0_980)),  new BCGaussianPrior(1.4, quad(0.2, 0.2)), new BCGaussianPrior(  12., quad(12., 10.))); */
    /* m.setPrior(free_amplitude(*m.model(), to(f_0_1370)), new BCGaussianPrior(1.3, quad(0.4, 0.2)), new BCGaussianPrior( -21., quad(15., 14.))); */
    /* m.setPrior(free_amplitude(*m.model(), to(f_0_1500)), new BCGaussianPrior(1.1, quad(0.3, 0.2)), new BCGaussianPrior( -44., quad(13., 16.))); */
    /* m.setPrior(free_amplitude(*m.model(), to(sigma)),    new BCGaussianPrior(3.7, quad(0.3, 0.2)), new BCGaussianPrior(  -3., quad(4.,   2.))); */

    /* // add shape parameters */
    /* m.addParameter("f_0_980_mass", static_pointer_cast<Flatte>(f_0_980->massShape())->mass(), 0.953 - 3 * 0.02, 0.953 + 3 * 0.02); */
    /* m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.953, 0.02)); */
    /* /\* m.addParameter("f_0_980_coupling", f_0_980_flatte->channels()[0].Coupling, 0.329 - 3 * 0.096, 0.329 + 3 * 0.096); *\/ */
    /* /\* m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.329, 0.096)); *\/ */
    /* m.addParameter("f_0_1370_mass", static_pointer_cast<BreitWigner>(f_0_1370->massShape())->mass(), 1.259 - 3 * 0.055, 1.259 + 3 * 0.055); */
    /* m.GetParameters().Back().SetPrior(new BCGaussianPrior(1.259, 0.055)); */
    /* m.addParameter("f_0_1370_width", dynamic_pointer_cast<BreitWigner>(f_0_1370->massShape())->width(), */
    /*                0.298 - 3 * 0.021, 0.298 + 3 * 0.021); */
    /* m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.298, 0.021)); */
    /* m.addParameter("sigma_mass", dynamic_pointer_cast<PoleMass>(sigma->massShape())->mass(), */
    /*                complex<double>(0.466, -0.223) - 3. * complex<double>(0.018, 0.028), */
    /*                complex<double>(0.466, -0.223) + 3. * complex<double>(0.018, 0.028)); */
    /* m.GetParameters()[m.GetParameters().Size() - 2].SetPrior(new BCGaussianPrior(0.466, 0.018)); */
    /* m.GetParameters().Back().SetPrior(new BCGaussianPrior(-0.223, 0.028)); */

    return m;
}

    
#endif
