// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include <DecayingParticle.h>
#include <FinalStateParticle.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <ZemachFormalism.h>

#include "bat_gen.h"
#include "models/d3pi.h"
#include "models/d3pi_phsp.h"
#include "models/d4pi.h"
#include "models/dkkpi.h"
#include "models/D_K0pi0pi0.h"

#include <chrono>
#include <ratio>

int main()
{
    yap::plainLogs(el::Level::Info);

    std::vector<bat_gen*> test_models = {
        // new bat_gen("D3PI_PHSP", d3pi_phsp(std::make_unique<yap::ZemachFormalism>()), 1.86961)
        // new bat_gen("D3PI", d3pi(std::make_unique<yap::ZemachFormalism>()), 1.86961)
        // new bat_gen("DKSPIPI_Zemach", D_K0pi0pi0(std::make_unique<yap::ZemachFormalism>()), 1.86961)
        // new bat_gen("DKSPIPI_Helicity", D_K0pi0pi0(std::make_unique<yap::HelicityFormalism>()), 1.86961)
        new bat_gen("DKKPI", dkkpi(std::make_unique<yap::ZemachFormalism>()), 1.86961)
        // new bat_gen("DKKPI", dkkpi(std::make_unique<yap::HelicityFormalism>()), 1.86961)
        // new bat_gen("D4PI", d4pi(), 1.86961)
    };

    for (bat_gen* m : test_models) {

        // open log file
        BCLog::OpenLog("output/" + m->GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

        // set precision
        m->SetPrecision(BCEngineMCMC::kMedium);
        m->SetNChains(4);
        m->SetMinimumEfficiency(0.15);
        m->SetMaximumEfficiency(0.35);

        m->SetNIterationsRun(static_cast<int>(4e6 / m->GetNChains()));

        m->WriteMarkovChain("output/" + m->GetSafeName() + "_mcmc.root", "RECREATE");

        // start timing:
        auto start = std::chrono::steady_clock::now();

        // run MCMC, marginalizing posterior
        m->MarginalizeAll(BCIntegrate::kMargMetropolis);

        // end timing
        auto end = std::chrono::steady_clock::now();

        // m->PrintAllMarginalized(m->GetSafeName() + "_plots.pdf");

        // timing:
        auto diff = end - start;
        auto ms = std::chrono::duration<double, std::micro>(diff).count();
        auto nevents = (m->GetNIterationsPreRun() + m->GetNIterationsRun()) * m->GetNChains();
        BCLog::OutSummary(std::string("Seconds = ") + std::to_string(ms / 1.e6) + " for " + std::to_string(nevents) + " iterations, " + std::to_string(m->likelihoodCalls()) + " calls");
        BCLog::OutSummary(std::to_string(ms / nevents) + " microsec / iteration");
        BCLog::OutSummary(std::to_string(ms / m->likelihoodCalls()) + " microsec / call");

        // close log file
        BCLog::OutSummary("Exiting");
        BCLog::CloseLog();
    }

    return 0;
}
