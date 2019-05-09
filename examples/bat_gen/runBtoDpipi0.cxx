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
#include "models/B0_Dpipi0.h"
#include "tools.h"

#include <chrono>

using namespace std;
using namespace yap;

int main()
{
    plainLogs(el::Level::Info);

    auto* m = new bat_gen("B0_Dpipi0", B0_Dpipi0(yap_model<ZemachFormalism>()), 5.27953);
    // auto* m = new bat_gen("B0_Dpipi0", B0_Dpipi0(yap_model<HelicityFormalism>()), 5.27963);

    // open log file
    BCLog::OpenLog("output/" + m->GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m->SetPrecision(BCEngineMCMC::kMedium);
    m->SetNChains(4);
    m->SetMinimumEfficiency(0.50);
    m->SetMaximumEfficiency(0.99);
    m->SetInitialPositionAttemptLimit(1e5);

    m->SetNIterationsRun(static_cast<int>(10e6 / m->GetNChains()));

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

    return 0;
}
