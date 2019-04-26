// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

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
#include "tools.h"

#include <chrono>

using namespace std;
using namespace yap;

int main()
{
    plainLogs(el::Level::Info);

    auto* m = new bat_gen("D3PI", d3pi(yap_model<ZemachFormalism>()), 1.86961);
    //    auto* m = new bat_gen("DKKPI", dkkpi(yap_model<HelicityFormalism>()), 1.86961),

    // open log file
    BCLog::OpenLog("output/" + m->GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m->SetPrecision(BCEngineMCMC::kMedium);
    m->SetNChains(4);
    m->SetMinimumEfficiency(0.95);
    m->SetMaximumEfficiency(0.99);
    m->SetInitialPositionAttemptLimit(1e5);

        m->SetNIterationsRun(static_cast<int>(1e5 / m->GetNChains()));

        m->WriteMarkovChain("output/" + m->GetSafeName() + "_mcmc.root", "RECREATE", true, false);

        // start timing:
        auto start = chrono::steady_clock::now();

        // run MCMC, marginalizing posterior
        m->MarginalizeAll(BCIntegrate::kMargMetropolis);

        // end timing
        auto end = chrono::steady_clock::now();

        // timing:
        auto diff = end - start;
        auto ms = chrono::duration<double, micro>(diff).count();
        auto nevents = (m->GetNIterationsPreRun() + m->GetNIterationsRun()) * m->GetNChains();
        BCLog::OutSummary(string("Seconds = ") + to_string(ms / 1.e6) + " for " + to_string(nevents) + " iterations, " + to_string(m->likelihoodCalls()) + " calls");
        BCLog::OutSummary(to_string(ms / nevents) + " microsec / iteration");
        BCLog::OutSummary(to_string(ms / m->likelihoodCalls()) + " microsec / call");

        // close log file
        BCLog::OutSummary("Exiting");
        BCLog::CloseLog();
    }

    return 0;
}
