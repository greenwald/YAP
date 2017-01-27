#include "mi_fit.h"

#include "models/d3pi.h"

#include <TCanvas.h>
#include <TEllipse.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLine.h>

#include <BAT/BCParameter.h>
#include <BAT/BCCauchyPrior.h>

#include <DecayingParticle.h>
#include <FreeAmplitude.h>
#include <logging.h>
#include <MeasuredBreakupMomenta.h>
#include <Model.h>
#include <Parameter.h>
#include <StepFunction.h>

// -----------------------
mi_fit::mi_fit(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs)
    : fit_base(name, d3pi_mi(std::move(M), 100), pcs)
{
    
    auto pipi_Swave = std::static_pointer_cast<yap::DecayingParticle>(particle(*model(), is_named("pipi_Swave")));
    auto step_function = std::static_pointer_cast<yap::StepFunction>(pipi_Swave->massShape());
    
    FAtoF2_ = free_amplitude(*model(), yap::to(particle(*model(), is_named("f_2"))));
    
    FAtoPiPiSwave_ = free_amplitude(*model(), to(pipi_Swave));

    /////////////////////////
    // AMP / PHASE MOTION
    // AddParameter("arg_pipi_Swave", -180, 180, "arg(#pi#pi S)");
    // GetParameters().Back().SetPriorConstant();
    // GetParameters().Back().Fix(-44.);

    /////////////////////////
    // REAL / IMAG
    *FAtoPiPiSwave_ = std::polar(1., 0.);
    
    double min_edge = step_function->steps().front()->lowEdge()->value();
    double max_edge = step_function->steps().back()->highEdge()->value();

    // loop over steps
    for (size_t i = 0; i < step_function->steps().size(); ++i) {
        double low_edge = step_function->lowEdges()[i]->value();
        double high_edge = step_function->lowEdges()[i + 1]->value();
        LOG(INFO) << "Adding bin " << i << " = [" << low_edge << ", " << high_edge << ") GeV^2 = ["
                  << static_cast<int>(1.e3 * (low_edge - min_edge) / (max_edge - min_edge)) * 0.1 << ", "
                  << static_cast<int>(1.e3 * (high_edge - min_edge) / (max_edge - min_edge)) * 0.1 << ") %";

        /////////////////////////
        // AMP / PHASE MOTION
        // 
        // AddParameter("amp(bin_" + std::to_string(i) + ")", 0, 25,
        //              "amp [" + std::to_string(low_edge) + ", " + std::to_string(high_edge) + "]");
        // GetParameters().Back().SetPriorConstant();
        // AddParameter("rel_phase(bin_" + std::to_string(i) + ")", 0, 18,
        //              "rel. phase [" + std::to_string(low_edge) + ", " + std::to_string(high_edge) + "]");
        // // GetParameters().Back().SetPriorConstant();
        // GetParameters().Back().SetPrior(new BCCauchyPrior(0, 9.));
        // if (i == 0)
        //     GetParameters().Back().Fix(0);

        /////////////////////////
        // REAL / IMAG
        AddParameter("real(bin_" + std::to_string(i) + ")", -4, 4,
                     "real [" + std::to_string(low_edge) + ", " + std::to_string(high_edge) + "]");
        AddParameter("imag(bin_" + std::to_string(i) + ")", -4, 4,
                     "imag [" + std::to_string(low_edge) + ", " + std::to_string(high_edge) + "]");

    }
    FAfromPiPiSwave_ = step_function->freeAmplitudes();

    AddParameter("real(f_2)", -4, 4);
    GetParameters().Back().Fix(real(FAtoF2_->value()));
    AddParameter("imag(f_2)", -4, 4);
    GetParameters().Back().Fix(imag(FAtoF2_->value()));
    
}

//-------------------------
void mi_fit::setParameters(const std::vector<double>& p)
{
    /////////////////////////
    // AMP / PHASE MOTION
    // set overall phase
    // *FAtoPiPiSwave_ = std::polar(1., yap::rad(p[0]));

    // // set step function amplitudes
    // double phase = 0;
    // for (size_t i = 0; i < FAfromPiPiSwave_.size(); ++i) {
    //     phase += yap::rad(p[i * 2 + 2]);
    //     *FAfromPiPiSwave_[i] = std::polar(p[i * 2 + 1], phase);
    // }

    for (size_t i = 0; i < FAfromPiPiSwave_.size(); ++i)
        *FAfromPiPiSwave_[i] = std::complex<double>(p[i * 2], p[i * 2 + 1]);
    *FAtoF2_ = std::complex<double>(p[FAfromPiPiSwave_.size() * 2], p[FAfromPiPiSwave_.size() + 2 + 1]);
    
    fit_base::setParameters(p);
}

//-------------------------
double mi_fit::LogAPrioriProbability(const std::vector<double>& P)
{
    double L = 0;
    for (size_t i = 0; i < P.size() / 2; ++i) {
        double A2 = pow(P[i * 2], 2) + pow(P[i * 2 + 1], 2);
        if (A2 > 16)
            return -std::numeric_limits<double>::infinity();
        // L -= 0.5 * log(A2);
    }
    return L;
}

//-------------------------
std::vector<double> mi_fit::truth() const
{
    // f_0(980)
    static double m_f0_980 = 0.965;
    static double g_f0_980 = 0.406;
    static double m_pi = 0.13957018;
    static double m_K  = 0.4936770;
    static auto a_f0_980 = std::polar(3., yap::rad(45.));
    static auto mw_f0_980 = 2. * g_f0_980
        * (std::sqrt(std::complex<double>(yap::measured_breakup_momenta::q2(m_f0_980 * m_f0_980, m_pi, m_pi)))
           + 2. * std::sqrt(std::complex<double>(yap::measured_breakup_momenta::q2(m_f0_980 * m_f0_980, m_K, m_K))))
        / m_f0_980;
    // f_0(1500)
    static double m_f0_1500 = 1.505;
    static double w_f0_1500 = 0.109;
    static auto a_f0_1500 = std::polar(3., yap::rad(-44.));
    // f_0(500)
    static std::complex<double> m_f0_500(0.470, -0.220);
    static auto a_f0_500 = std::polar(1., yap::rad(-3.));

    auto pipi_Swave = std::static_pointer_cast<yap::DecayingParticle>(particle(*model(), is_named("pipi_Swave")));
    const auto& low_edges = std::static_pointer_cast<yap::StepFunction>(pipi_Swave->massShape())->lowEdges();
    std::vector<double> P;
    for (size_t i = 0; i < low_edges.size() - 1; ++i) {
        double m2 = 0.5 * (low_edges[i]->value() + low_edges[i + 1]->value());

        std::complex<double> amp = 0;
        // f_0(980)
        auto w_f0_980 = std::sqrt(std::complex<double>(yap::measured_breakup_momenta::q2(m2, m_pi, m_pi)));
        w_f0_980     += 2. * std::sqrt(std::complex<double>(yap::measured_breakup_momenta::q2(m2, m_K, m_K)));
        amp += a_f0_980 * mw_f0_980 / (m_f0_980 * m_f0_980 - m2 - 2. * 1_i * g_f0_980 * w_f0_980 / sqrt(m2));
        // f_0(1500)
        amp += a_f0_1500 * m_f0_1500 * w_f0_1500 / (m_f0_1500 * m_f0_1500 - m2 - 1_i * m_f0_1500 * w_f0_1500);
        // f_0(500)
        amp += a_f0_500 * 2. * real(m_f0_500) * imag(m_f0_500) / (m_f0_500 * m_f0_500 - m2);
        
        P.push_back(real(amp));
        P.push_back(imag(amp));
    }
    return P;
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

    auto par_tr = truth();
    
    TH1D* h1_amp = new TH1D("h1_amp", ";M^{2} [GeV^{2}];amp(A);", bins.size() - 1, &bins[0]);
    TH1D* h1_arg = new TH1D("h1_arg", ";M^{2} [GeV^{2}];arg(A);", bins.size() - 1, &bins[0]);
    TH1D* h1_rel = new TH1D("h1_rel", ";M^{2} [GeV^{2}];relative phase;", bins.size() - 1, &bins[0]);

    TH1D* h1_amp_tr = new TH1D("h1_amp_tr", ";M^{2} [GeV^{2}];amp(A);", bins.size() - 1, &bins[0]);
    TH1D* h1_arg_tr = new TH1D("h1_arg_tr", ";M^{2} [GeV^{2}];arg(A);", bins.size() - 1, &bins[0]);
    TH1D* h1_rel_tr = new TH1D("h1_rel_tr", ";M^{2} [GeV^{2}];relative phase;", bins.size() - 1, &bins[0]);

    TGraphErrors* g_argand = new TGraphErrors();
    g_argand->SetTitle(";re(A);im(A)");

    TGraph* g_argand_tr = new TGraph();
    g_argand_tr->SetTitle(";re(A);im(A)");

    TH1D* h1_re = new TH1D("h1_re",   ";M [GeV];re(A);", bins.size() - 1, &bins[0]);
    TH1D* h1_im = new TH1D("h1_im",   ";M [GeV];im(A);", bins.size() - 1, &bins[0]);

    // double phase = 0;
    std::complex<double> amp_old = 0;
    std::complex<double> amp_tr_old(par_tr[0], par_tr[1]);
    for (size_t i = 0; i < step_function->steps().size(); ++i) {
        // phase += pars[i * 2 + 2];

        // AMP / REL PHASE
        // auto amp = std::polar(pars[i * 2 + 1], yap::rad(phase));
        // REAL / IMAG
        auto amp = std::complex<double>(pars[i * 2], pars[i * 2 + 1]);
        auto amp_tr = std::complex<double>(par_tr[i * 2], par_tr[i * 2 + 1]);
        h1_amp->SetBinContent(i + 1, abs(amp));
        h1_arg->SetBinContent(i + 1, yap::deg(arg(amp)));
        h1_rel->SetBinContent(i + 1, arg(amp - amp_old));
        h1_re->SetBinContent(i + 1, real(amp));
        h1_im->SetBinContent(i + 1, imag(amp));

        g_argand->SetPoint(i, real(amp), imag(amp));
        
        h1_amp_tr->SetBinContent(i + 1, abs(amp_tr));
        h1_arg_tr->SetBinContent(i + 1, yap::deg(arg(amp_tr)));
        h1_rel_tr->SetBinContent(i + 1, yap::deg(arg(amp_tr) - arg(amp_tr_old)));
        g_argand_tr->SetPoint(i, real(amp_tr), imag(amp_tr));

        amp_tr_old = amp_tr;
        amp_old = amp;
        
        if (!uncs.empty()) {
            if (!std::isfinite(uncs[i * 2]) or !std::isfinite(uncs[i * 2 + 1]))
                continue;

            // AMP / REL PHASE
            // phase_unc2 += pow(uncs[i * 2 + 2], 2);
            // amp_unc = uncs[i * 2 + 1];
            // REAL / IMAG
            double amp_unc = sqrt(pow(uncs[i * 2] * pars[i * 2], 2) + pow(uncs[i * 2 + 1] * pars[i * 2 + 1], 2)) / abs(amp);
            double phase_unc = sqrt(pow(cos(arg(amp)) * uncs[i * 2], 2) + pow(sin(arg(amp)) * uncs[i * 2 + 1], 2));

            g_argand->SetPointError(i, uncs[i * 2], uncs[i * 2 + 1]);

            h1_amp->SetBinError(i + 1, amp_unc);
            h1_arg->SetBinError(i + 1, sqrt(phase_unc));
            // h1_rel->SetBinError(i + 1, h1_rel->GetBinError);
            
        }
    }

    TCanvas* C = new TCanvas("C", "C");
    C->Divide(2, 2);

    C->cd(1);
    h1_amp->SetStats(false);
    h1_amp->Draw();
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
    // TEllipse* e_argand = new TEllipse(0, 3. / m_f0 / w_f0 / 2, 3. / m_f0 / w_f0 / 2);
    // e_argand->SetFillStyle(0);
    // e_argand->SetLineColor(2);
    // e_argand->Draw();
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

    // delete h1_amp;
    // delete f_amp;
    // delete h1_arg;
    // delete h1_rel;
    // delete h1_re;
    // delete h1_im;
    // delete g_argand;
    // delete e_argand;
    // delete h1_amp_tr;
    // delete h1_arg_tr;
    // delete h1_rel_tr;
    // delete g_argand;
    // delete L;
    // delete X;
    // delete Y;
    // delete X_tr;
    // delete Y_tr;
    // delete C;
}
