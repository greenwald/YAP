// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "fit_fitFraction.h"
#include "models/d3pi.h"

#include <DecayTree.h>
#include <Filters.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MassRange.h>
#include <Model.h>
#include <PHSP.h>
#include <ZemachFormalism.h>

#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <algorithm>
#include <array>
#include <chrono>
#include <random>

const double quad(std::vector<double> S)
{ return sqrt(std::accumulate(S.begin(), S.end(), 0., [](double a, double s) {return a + s * s;})); }

template <typename ... Types>
constexpr double quad(double s0, Types ... additional)
{ return quad({s0, additional...}); }

int main()
{
    yap::plainLogs(el::Level::Info);

    // create bat_fit object
    fit_fitFraction m("D3PI_frac_fit", d3pi(std::make_unique<yap::ZemachFormalism>()));

    double D_mass = 1.86961;

    m.GetParameter("N_1").Fix(1);

    // find particles
    auto D        = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("D+")));
    auto rho      = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("rho")));
    auto f_2      = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_2")));
    auto f_0_980  = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_0_980")));
    auto f_0_1370 = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_0_1370")));
    auto f_0_1500 = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("f_0_1500")));
    auto sigma    = std::static_pointer_cast<Resonance>(particle(*m.model(), is_named("sigma")));

    // set fit fractions to fit
    m.setFitFraction(decay_tree(*D, yap::to(rho)),      20e-2,   quad(2.3e-2, 0.9e-2));
    m.setFitFraction(decay_tree(*D, yap::to(f_2)),      18.2e-2, quad(2.6e-2, 0.7e-2));
    m.setFitFraction(decay_tree(*D, yap::to(f_0_980)),  4.1e-2,  quad(0.9e-2, 0.3e-2));
    m.setFitFraction(decay_tree(*D, yap::to(f_0_1370)), 2.6e-2,  quad(1.8e-2, 0.6e-2));
    m.setFitFraction(decay_tree(*D, yap::to(f_0_1500)), 3.4e-2,  quad(1.0e-2, 0.8e-2));
    m.setFitFraction(decay_tree(*D, yap::to(sigma)),    41.8e-2, quad(1.4e-2, 2.5e-2));

    // set free amplitude parameters of fit
    m.fix(free_amplitude(*m.model(), yap::to(rho)), 1., 0.);
    m.setPrior(free_amplitude(*m.model(), yap::to(f_2)),      new BCGaussianPrior(2.1, quad(0.2, 0.1)), new BCGaussianPrior(-123., quad(6.,   3.)));
    m.setPrior(free_amplitude(*m.model(), yap::to(f_0_980)),  new BCGaussianPrior(1.4, quad(0.2, 0.2)), new BCGaussianPrior(  12., quad(12., 10.)));
    m.setPrior(free_amplitude(*m.model(), yap::to(f_0_1370)), new BCGaussianPrior(1.3, quad(0.4, 0.2)), new BCGaussianPrior( -21., quad(15., 14.)));
    m.setPrior(free_amplitude(*m.model(), yap::to(f_0_1500)), new BCGaussianPrior(1.1, quad(0.3, 0.2)), new BCGaussianPrior( -44., quad(13., 16.)));
    m.setPrior(free_amplitude(*m.model(), yap::to(sigma)),    new BCGaussianPrior(3.7, quad(0.3, 0.2)), new BCGaussianPrior(  -3., quad(4.,   2.)));

    // add shape parameters
    m.addParameter("f_0_980_mass", std::static_pointer_cast<yap::Flatte>(f_0_980->massShape())->mass(), 0.953 - 3 * 0.02, 0.953 + 3 * 0.02);
    m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.953, 0.02));
    // m.addParameter("f_0_980_coupling", f_0_980_flatte->channels()[0].Coupling, 0.329 - 3 * 0.096, 0.329 + 3 * 0.096);
    // m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.329, 0.096));
    m.addParameter("f_0_1370_mass", std::static_pointer_cast<yap::BreitWigner>(f_0_1370->massShape())->mass(), 1.259 - 3 * 0.055, 1.259 + 3 * 0.055);
    m.GetParameters().Back().SetPrior(new BCGaussianPrior(1.259, 0.055));
    m.addParameter("f_0_1370_width", std::dynamic_pointer_cast<yap::BreitWigner>(f_0_1370->massShape())->width(),
                   0.298 - 3 * 0.021, 0.298 + 3 * 0.021);
    m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.298, 0.021));
    m.addParameter("sigma_mass", std::dynamic_pointer_cast<yap::PoleMass>(sigma->massShape())->mass(),
                   std::complex<double>(0.466, -0.223) - 3. * std::complex<double>(0.018, 0.028),
                   std::complex<double>(0.466, -0.223) + 3. * std::complex<double>(0.018, 0.028));
    m.GetParameters()[m.GetParameters().Size() - 2].SetPrior(new BCGaussianPrior(0.466, 0.018));
    m.GetParameters().Back().SetPrior(new BCGaussianPrior(-0.223, 0.028));

    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(D_mass, m.axes(), m.model()->finalStateParticles()));

    // generate integration data
    std::mt19937 g(0);
    std::generate_n(std::back_inserter(m.integralData()), 10000,
                    std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()),
                              D_mass, m.axes(), m2r, g,
                              std::numeric_limits<unsigned>::max()));
    LOG(INFO) << "Created " << m.integralData().size() << " data points (" << (m.integralData().bytes() * 1.e-6) << " MB)";
    // partition integration data
    m.integralPartitions() = yap::DataPartitionBlock::create(m.integralData(), 2);

    // open log file
    BCLog::OpenLog("output/" + m.GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m.SetPrecision(BCEngineMCMC::kMedium);
    m.SetNIterationsPreRunMax(1e6);
    m.SetNChains(4);
    // m.SetMinimumEfficiency(0.85);
    // m.SetMaximumEfficiency(0.99);

    m.SetNIterationsRun(static_cast<int>(2e4 / m.GetNChains()));

    m.WriteMarkovChain("output/" + m.GetSafeName() + "_mcmc.root", "RECREATE");

    // start timing:
    auto start = std::chrono::steady_clock::now();

    // run MCMC, marginalizing posterior
    m.MarginalizeAll(BCIntegrate::kMargMetropolis);

    // end timing
    auto end = std::chrono::steady_clock::now();

    for (size_t i = 0; i < D->decayTrees().at(0).size(); ++i)
        std::cout << i << " = " << to_string(*D->decayTrees().at(0)[i]) << std::endl;

    m.PrintSummary();
    m.PrintAllMarginalized("output/" + m.GetSafeName() + "_plots.pdf", 2, 2);
    m.SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDetailedPosterior);
    m.PrintKnowledgeUpdatePlots("output/" + m.GetSafeName() + "_update.pdf", 2, 2, true);

    // timing:
    auto diff = end - start;
    auto ms = std::chrono::duration<double, std::micro>(diff).count();
    auto nevents = (m.GetNIterationsPreRun() + m.GetNIterationsRun()) * m.GetNChains();
    BCLog::OutSummary(std::string("Seconds = ") + std::to_string(ms / 1.e6) + " for " + std::to_string(nevents) + " iterations, " + std::to_string(m.likelihoodCalls()) + " calls");
    BCLog::OutSummary(std::to_string(ms / nevents) + " microsec / iteration");
    BCLog::OutSummary(std::to_string(ms / m.likelihoodCalls()) + " microsec / call");

    // // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    return 0;
}
