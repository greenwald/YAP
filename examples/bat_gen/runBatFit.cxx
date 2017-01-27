// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_fit.h"
#include "fit_base.h"
#include "mi_fit.h"
//#include "hist.h"
#include "models/d3pi.h"
#include "models/d4pi.h"
#include "models/dkkpi.h"
#include "tools.h"

#include <Exceptions.h>
#include <HelicityFormalism.h>
#include <MassRange.h>
#include <PHSP.h>
#include <RelativisticBreitWigner.h>
#include <ZemachFormalism.h>
#include <logging.h>
#include <make_unique.h>

#include <BAT/BCAux.h>
#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <random>

int main()
{
    yap::plainLogs(el::Level::Info);

    unsigned i_model = 3;

    std::string model_name;
    switch (i_model) {
        case 0:
            model_name = "D3PI";
            break;
        case 1:
            model_name = "DKKPI";
            break;
        case 2:
            model_name = "D4PI";
            break;
        case 3:
            model_name = "D3PI";
            break;
        default:
            throw yap::exceptions::Exception("No valid model requested.", "runBatFit::main");
    }

    // open file
    auto file = TFile::Open(("output/" + model_name + "_mcmc.root").data(), "READ");
    if (file->IsZombie())
        throw yap::exceptions::Exception("could not open file", "main");

    TTree* t_mcmc = nullptr;
    file->GetObject((model_name + "_mcmc").data(), t_mcmc);
    if (!t_mcmc)
        throw yap::exceptions::Exception("could not retrieve mcmc tree", "main");

    TTree* t_pars = nullptr;
    file->GetObject((model_name + "_parameters").data(), t_pars);
    if (!t_pars)
        throw yap::exceptions::Exception("could not retrieve mcmc tree", "main");

    // create model
    fit_base* m = nullptr;
    double D_mass(0);
    switch (i_model) {
        case 0:
            m = new bat_fit(d3pi_fit(model_name + "_fit", yap_model<yap::ZemachFormalism>(), find_mass_axes(*t_pars)));
            D_mass = 1.86961; // D+
            break;
        case 1:
            m = new bat_fit(dkkpi_fit(model_name + "_fit", yap_model<yap::HelicityFormalism>(), find_mass_axes(*t_pars)));
            D_mass = 1.86961; // D+
            break;
        case 2:
            m = new bat_fit(d4pi_fit(model_name + "_fit"));
            D_mass = 1.8648400; // D0
            break;
        case 3:
            m = new mi_fit(model_name + "_fit", yap_model<yap::ZemachFormalism>(), find_mass_axes(*t_pars));
            D_mass = 1.86961; // D+
            break;
        default:
            m = nullptr;
    }
    if (!m)
        throw exceptions::Exception("No model loaded", "runBatFit::main");

    // load fit data and partition it
    LOG(INFO) << "Load data";
    load_data(m->fitData(), *m->model(), m->axes(), D_mass, *t_mcmc, 5e4);
    m->fitPartitions() = yap::DataPartitionBlock::create(m->fitData(), 4);

    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(D_mass, m->axes(), m->model()->finalStateParticles()));

    // generate integration data
    std::mt19937 g(0);
    unsigned n_integrationPoints = 1e5;
    if (true) {
        m->integrationPointGenerator() = std::bind(yap::phsp<std::mt19937>, std::cref(*m->model()), D_mass, m->axes(), m2r, g, std::numeric_limits<unsigned>::max());
        // m->setNIntegrationPoints(4e4, 4e4);
        m->setNIntegrationPoints(n_integrationPoints, 1e4, 4);
        LOG(INFO) << "Generating integration points on the fly";
    } else {
        // get FSP mass ranges
        auto m2r = yap::squared(mass_range(D_mass, m->axes(), m->model()->finalStateParticles()));
        
        // generate integration data
        std::generate_n(std::back_inserter(m->integralData()), n_integrationPoints,
                        std::bind(yap::phsp<std::mt19937>, std::cref(*m->model()), D_mass, m->axes(), m2r, g, std::numeric_limits<unsigned>::max()));
        m->integralPartitions() = yap::DataPartitionBlock::create(m->integralData(), 4);
        LOG(INFO) << "Created " << m->integralData().size() << " data points (" << (m->integralData().bytes() * 1.e-6) << " MB)";
    }

    // TH2D* h2_fit_data = hist2(*m->model()->fourMomenta(), m->axes(), m2r, m->fitData());
    // TH2D* h2_int_data = hist2(*m->model()->fourMomenta(), m->axes(), m2r, m->integralData());

    // TH1D* h1_fit_data = hist1(*m->model()->fourMomenta(), m->axes()[0], m2r[0], m->fitData());
    // TH1D* h1_int_data = hist1(*m->model()->fourMomenta(), m->axes()[0], m2r[0], m->integralData());

    // TCanvas* C = new TCanvas("C","C");
    // C->Divide(1,2);
    // C->cd(1);
    // h2_fit_data->Draw("colz");
    // C->cd(2);
    // h2_int_data->Draw("colz");
    // C->Print("data.pdf");

    // return 0;

    // open log file
    BCLog::OpenLog("output/" + m->GetSafeName() + "_log.txt", BCLog::detail, BCLog::detail);

    // set precision
    m->SetPrecision(BCEngineMCMC::kMedium);
    m->SetNIterationsPreRunMax(1e6);
    m->SetNChains(4);
    // m->SetMinimumEfficiency(0.85);
    // m->SetMaximumEfficiency(0.99);

    m->SetNIterationsRun(static_cast<int>(50e3 / m->GetNChains()));

    m->SetInitialPositionAttemptLimit(100000);
    
    // m->WriteMarkovChain("output/" + m->GetSafeName() + "_mcmc.root", "RECREATE", true, false);

    // start timing:
    auto start = std::chrono::steady_clock::now();

    // run MCMC, marginalizing posterior
    m->MCMCUserInitialize();
    // m->MarginalizeAll(BCIntegrate::kMargMetropolis);

    BCLog::OutSummary("Optimizing ...");

    std::vector<double> guess;
    if (dynamic_cast<mi_fit*>(m))
        guess = static_cast<mi_fit*>(m)->truth();
    guess.resize(m->GetParameters().Size(), 1);
    for (size_t i = 0; i < guess.size(); ++i)
        if (m->GetParameter(i).Fixed())
            guess[i] = m->GetParameter(i).GetFixedValue();

    m->FindMode(guess);
    // BCLog::OutSummary("Optimizing more ...");
    // m->FindMode(m->GetBestFitParameters());
    BCLog::OutSummary("... done optimizing");

    // end timing
    auto end = std::chrono::steady_clock::now();

    m->PrintSummary();
    // m->PrintAllMarginalized("output/" + m->GetSafeName() + "_plots.pdf", 2, 2);

    // m->SetNIterationsRun(static_cast<int>(10e3 / m->GetNChains()));
    // // m->SetKnowledgeUpdateDrawingStyle(BCAux::kKnowledgeUpdateDetailedPosterior);
    // m->PrintKnowledgeUpdatePlots("output/" + m->GetSafeName() + "_update.pdf", 2, 2, true);

    // timing:
    auto diff = end - start;
    auto ms = std::chrono::duration<double, std::micro>(diff).count();
    auto nevents = (m->GetNIterationsPreRun() + m->GetNIterationsRun()) * m->GetNChains();
    BCLog::OutSummary(std::string("Seconds = ") + std::to_string(ms / 1.e6) + " for " + std::to_string(nevents) + " iterations, " + std::to_string(m->likelihoodCalls()) + " calls");
    BCLog::OutSummary(std::to_string(ms / nevents) + " microsec / iteration");
    BCLog::OutSummary(std::to_string(ms / m->likelihoodCalls()) + " microsec / call");

    // Print best fit:
    // for Model Independent
    if (dynamic_cast<mi_fit*>(m))
        dynamic_cast<mi_fit*>(m)->printSwave("output/mi_result.pdf", m->GetBestFitParameters(), m->GetBestFitParameterErrors());

    // // close log file
    BCLog::OutSummary("Exiting");
    BCLog::CloseLog();

    delete m;

    return 0;
}
