#include "mi_fit.h"

#include "models/d3pi.h"

#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TLine.h>

#include <BAT/BCParameter.h>
#include <BAT/BCCauchyPrior.h>

#include <DecayingParticle.h>
#include <FreeAmplitude.h>
#include <logging.h>
#include <Model.h>
#include <Parameter.h>
#include <StepFunction.h>

// -----------------------
mi_fit::mi_fit(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs)
    : fit_base(name, d3pi_mi(std::move(M)), pcs)
{
    
    auto pipi_Swave = std::static_pointer_cast<yap::DecayingParticle>(particle(*model(), is_named("pipi_Swave")));
    auto step_function = std::static_pointer_cast<yap::StepFunction>(pipi_Swave->massShape());

    FAtoPiPiSwave_ = free_amplitude(*model(), to(pipi_Swave));
    AddParameter("arg_pipi_Swave", -180, 180, "arg(#pi#pi S)");
    GetParameters().Back().SetPriorConstant();
    GetParameters().Back().Fix(-44.);

    double min_edge = step_function->steps().front()->lowEdge()->value();
    double max_edge = step_function->steps().back()->highEdge()->value();

    // loop over steps
    for (size_t i = 0; i < step_function->steps().size(); ++i) {
        double low_edge = step_function->lowEdges()[i]->value();
        double high_edge = step_function->lowEdges()[i + 1]->value();
        LOG(INFO) << "Adding bin " << i << " = [" << low_edge << ", " << high_edge << ") GeV^2 = ["
                  << static_cast<int>(1.e3 * (low_edge - min_edge) / (max_edge - min_edge)) * 0.1 << ", "
                  << static_cast<int>(1.e3 * (high_edge - min_edge) / (max_edge - min_edge)) * 0.1 << ") %";
        AddParameter("amp(bin_" + std::to_string(i) + ")", 0, 25,
                     "amp [" + std::to_string(low_edge) + ", " + std::to_string(high_edge) + "]");
        GetParameters().Back().SetPriorConstant();
        AddParameter("rel_phase(bin_" + std::to_string(i) + ")", 0, 18,
                     "rel. phase [" + std::to_string(low_edge) + ", " + std::to_string(high_edge) + "]");
        // GetParameters().Back().SetPriorConstant();
        GetParameters().Back().SetPrior(new BCCauchyPrior(0, 9.));
        if (i == 0)
            GetParameters().Back().Fix(0);
    }
    FAfromPiPiSwave_ = step_function->freeAmplitudes();
    
}

//-------------------------
void mi_fit::setParameters(const std::vector<double>& p)
{
    // set overall phase
    *FAtoPiPiSwave_ = std::polar(1., yap::rad(p[0]));

    // set step function amplitudes
    double phase = 0;
    for (size_t i = 0; i < FAfromPiPiSwave_.size(); ++i) {
        phase += yap::rad(p[i * 2 + 2]);
        *FAfromPiPiSwave_[i] = std::polar(p[i * 2 + 1], phase);
    }

    fit_base::setParameters(p);
}

//-------------------------
void mi_fit::printSwave(const std::string& filename, const std::vector<double>& pars, const std::vector<double>& uncs)
{
    if (pars.empty()) {
        if (GetBestFitParameters().empty())
            throw;
        printSwave(filename, GetBestFitParameters(), GetBestFitParameterErrors());
        return;
    }

    auto pipi_Swave = std::static_pointer_cast<yap::DecayingParticle>(particle(*model(), is_named("pipi_Swave")));
    auto step_function = std::static_pointer_cast<yap::StepFunction>(pipi_Swave->massShape());
    std::vector<double> bins;
    bins.reserve(step_function->lowEdges().size());
    for (const auto& l : step_function->lowEdges())
        bins.push_back(l->value());

    TH1D* h1_amp = new TH1D("h1_amp", ";M^{2} [GeV^{2}];amp(A);", bins.size() - 1, &bins[0]);
    TH1D* h1_arg = new TH1D("h1_arg", ";M^{2} [GeV^{2}];arg(A);", bins.size() - 1, &bins[0]);
    TH1D* h1_rel = new TH1D("h1_rel", ";M^{2} [GeV^{2}];relative phase;", bins.size() - 1, &bins[0]);

    TH1D* h1_amp_tr = new TH1D("h1_amp_tr", ";M^{2} [GeV^{2}];amp(A);", bins.size() - 1, &bins[0]);
    TH1D* h1_arg_tr = new TH1D("h1_arg_tr", ";M^{2} [GeV^{2}];arg(A);", bins.size() - 1, &bins[0]);
    TH1D* h1_rel_tr = new TH1D("h1_rel_tr", ";M^{2} [GeV^{2}];relative phase;", bins.size() - 1, &bins[0]);

    TGraph* g_argand = new TGraph();
    g_argand->SetTitle(";re(A);im(A)");

    TGraph* g_argand_tr = new TGraph();
    g_argand_tr->SetTitle(";re(A);im(A)");

    TH1D* h1_re = new TH1D("h1_re",   ";M [GeV];re(A);", bins.size() - 1, &bins[0]);
    TH1D* h1_im = new TH1D("h1_im",   ";M [GeV];im(A);", bins.size() - 1, &bins[0]);

    double m_f0 = 1.505;
    double w_f0 = 0.109;
    double a_f0 = 3.;

    double phase = 0;
    double phase_unc = 0;
    std::complex<double> amp_tr_old = 0;
    for (size_t i = 0; i < step_function->steps().size(); ++i) {
        phase += pars[i * 2 + 2];

        auto amp = std::polar(pars[i * 2 + 1], yap::rad(phase));
        
        h1_amp->SetBinContent(i + 1, abs(amp));
        h1_arg->SetBinContent(i + 1, yap::deg(arg(amp)));
        h1_rel->SetBinContent(i + 1, pars[i * 2 + 2]);
        h1_re->SetBinContent(i + 1, real(amp));
        h1_im->SetBinContent(i + 1, imag(amp));

        g_argand->SetPoint(i, real(amp), imag(amp));

        auto amp_tr = a_f0 / (m_f0 * m_f0 - h1_amp->GetXaxis()->GetBinCenter(i + 1) - 1_i * m_f0 * w_f0);
        if (i == 0)
            amp_tr_old = amp_tr;
        
        h1_amp_tr->SetBinContent(i + 1, abs(amp_tr));
        h1_arg_tr->SetBinContent(i + 1, yap::deg(arg(amp_tr)));
        h1_rel_tr->SetBinContent(i + 1, yap::deg(arg(amp_tr) - arg(amp_tr_old)));
        g_argand_tr->SetPoint(i, real(amp_tr), imag(amp_tr));

        amp_tr_old = amp_tr;
        
        if (!uncs.empty()) {
            phase_unc += pow(uncs[i * 2 + 2], 2);
            h1_amp->SetBinError(i + 1, uncs[i * 2 + 1]);
            h1_arg->SetBinError(i + 1, sqrt(phase_unc));
            h1_rel->SetBinError(i + 1, uncs[i * 2 + 2]);
        }
    }

    TF1* f_amp = new TF1("f_amp", "[2] / sqrt(([0]^2 - x)^2 + [0]^2 * [1]^2)", bins[0], bins.back());
    f_amp->SetParameters(m_f0, w_f0, a_f0);

    TCanvas* C = new TCanvas("C", "C");
    C->Divide(2, 2);

    C->cd(1);
    h1_amp->SetStats(false);
    h1_amp->Draw();
    f_amp->Draw("same");
    h1_amp_tr->SetLineColor(2);
    h1_amp_tr->Draw("same");

    C->cd(2);
    h1_arg->SetStats(false);
    h1_arg->Draw();
    h1_arg_tr->SetLineColor(2);
    h1_arg_tr->Draw("same");

    C->cd(3);
    g_argand->SetMarkerStyle(8);
    g_argand->Draw("ALP");
    TEllipse* e_argand = new TEllipse(0, 3. / m_f0 / w_f0 / 2, 3. / m_f0 / w_f0 / 2);
    e_argand->SetFillStyle(0);
    e_argand->SetLineColor(2);
    e_argand->Draw();
    g_argand_tr->SetMarkerStyle(4);
    g_argand_tr->SetMarkerColor(2);
    g_argand_tr->SetLineColor(2);
    g_argand_tr->Draw("sameP");
    TLine* L = new TLine();
    L->SetLineColor(16);
    L->SetLineStyle(2);
    double* X = g_argand->GetX();
    double* Y = g_argand->GetY();
    double* X_tr = g_argand_tr->GetX();
    double* Y_tr = g_argand_tr->GetY();
    for (int i = 0; i < g_argand->GetN(); ++i)
        L->DrawLine(X[i], Y[i], X_tr[i], Y_tr[i]);

    // h1_re->SetStats(false);
    // h1_re->Draw();

    C->cd(4);
    h1_rel->SetStats(false);
    h1_rel->Draw();
    h1_rel_tr->SetLineColor(2);
    h1_rel_tr->Draw("same");
   // h1_im->SetStats(false);
    // h1_im->Draw();

    C->Print(filename.data());

    delete h1_amp;
    delete f_amp;
    delete h1_arg;
    delete h1_rel;
    delete h1_re;
    delete h1_im;
    delete g_argand;
    delete e_argand;
    delete h1_amp_tr;
    delete h1_arg_tr;
    delete h1_rel_tr;
    delete g_argand;
    delete L;
    delete X;
    delete Y;
    delete X_tr;
    delete Y_tr;
    delete C;
}
