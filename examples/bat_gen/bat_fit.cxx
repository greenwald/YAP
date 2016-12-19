#include "bat_fit.h"

#include "ConstantPrior.h"

#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <SpinAmplitude.h>
#include <VariableStatus.h>

#include <BAT/BCPrior.h>

bool begins_with(const std::string& haystack, const std::string& needle)
{
    return haystack.size() >= needle.size() and std::equal(needle.begin(), needle.end(), haystack.begin());
}

// -----------------------
bat_fit::bat_fit(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs)
    : fit_base(name, std::move(M), pcs)
{
    // loop over all free amplitudes
    for (const auto& fa : free_amplitudes(*model())) {
        // ignore fixed free amplitudes
        if (fa->variableStatus() == yap::VariableStatus::fixed)
            continue;

        auto fa_name = to_string(*fa->decayChannel())
            + " L = " + std::to_string(fa->spinAmplitude()->L())
            + " S = " + yap::spin_to_string(fa->spinAmplitude()->twoS());

        auto re_fa_name = "real(" + fa_name + ")";
        auto im_fa_name = "imag(" + fa_name + ")";
        auto amp_fa_name = "amp(" + fa_name + ")";
        auto arg_fa_name = "arg(" + fa_name + ")";

        unsigned n = 0;
        // check how many are there
        for (unsigned i = 0; i < GetParameters().Size(); ++i)
            if (begins_with(GetParameters()[i].GetName(), re_fa_name) == 0)
                ++n;
        if (n > 0) {
            // if only one found, relabel it
            if (n == 1) {
                GetParameter(re_fa_name).SetName(re_fa_name + "_0");
                GetParameter(im_fa_name).SetName(im_fa_name + "_0");
                GetParameter(amp_fa_name).SetName(amp_fa_name + "_0");
                GetParameter(arg_fa_name).SetName(arg_fa_name + "_0");
            }
            // add index to current parameter
            re_fa_name += "_" + std::to_string(n);
            im_fa_name += "_" + std::to_string(n);
            amp_fa_name += "_" + std::to_string(n);
            arg_fa_name += "_" + std::to_string(n);
        }

        // // add real parameter
        // AddParameter(re_fa_name, -2 * abs(fa->value()), 2 * abs(fa->value()));
        // // add imag parameter
        // AddParameter(im_fa_name, -2 * abs(fa->value()), 2 * abs(fa->value()));

        // AbsPriors_.push_back(new ConstantPrior(0, 2 * abs(fa->value())));
        // ArgPriors_.push_back(new ConstantPrior(-180, 180));

        // // add amplitude observable
        // AddObservable(amp_fa_name, 0, 2 * abs(fa->value()));
        // // add phase observable
        // AddObservable(arg_fa_name, -180, 180);

        AddParameter(amp_fa_name, 0, 2 * abs(fa->value()));
        AddParameter(arg_fa_name, 0, 360);
        
        // add free amplitude to list
        FreeAmplitudes_.push_back(fa);
    }
    
    // // add observables for all fit fractions
    // int N = std::accumulate(Integral_.integrals().begin(), Integral_.integrals().end(), 0,
    //                         [](int n, const yap::IntegralMap::value_type& v)
    //                         {return n + v.second.decayTrees().size();});
    // if (N > 1) {
    // for (const auto& mci : Integral_.integrals())
    //     for (const auto& dt : mci.Integral.decayTrees()) {
    //         DecayTrees_.push_back(dt);

    //         std::string str = "fit_frac(" + to_string(*dt) + ")";
    //         std::replace(str.begin(), str.end(), '-', 'm'); // - will be omitted by BATs "safe name", and when decay channels only differ in some spin projections (-1 vs 1), the safe name would be identical
    //         std::replace(str.begin(), str.end(), '\n', ';'); // make into one line
    //         std::replace(str.begin(), str.end(), '\t', ' ');

    //         unsigned n = 0;
    //         // check how many are there
    //         for (unsigned i = 0; i < GetObservables().Size(); ++i)
    //             if (begins_with(GetObservables()[i].GetName(), str) == 0)
    //                 ++n;
            
    //         if (n > 0) {
    //             if (n == 1)
    //                 GetObservable(str).SetName(str + "_0");
    //             str += "_" + std::to_string(n);
    //         }

    //         AddObservable(str, 0, 1.1);
    //         GetObservables().Back().SetNbins(1000);
    //     }
    // }

    FirstParameter_ = GetParameters().Size();
    FirstObservable_ = GetObservables().Size();
}

//-------------------------
void bat_fit::setParameters(const std::vector<double>& p)
{
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i)
        *FreeAmplitudes_[i] = std::complex<double>(p[i * 2], p[i * 2 + 1]);

    if (!Parameters_.empty())
        yap::set_values(Parameters_.begin(), Parameters_.end(), p.begin() + FirstParameter_, p.end());

    fit_base::setParameters(p);

    // calculate fit fractions
    // if (!CalculatedFitFractions_.empty()) {
    // unsigned c = GetCurrentChain();
    // size_t i = 0;
    // for (const auto& mci : Integral_.integrals()) {
    //     auto ff = fit_fractions(mci.Integral);
    //     for (const auto& f : ff)
    //         CalculatedFitFractions_[c][i++] = f;
    // }
    // }
}

//-------------------------
double bat_fit::LogAPrioriProbability(const std::vector<double>& p)
{
    double logP = 0;
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i) {
        if (GetParameter(i * 2).Fixed() or GetParameter(i * 2 + 1).Fixed())
            continue;
        auto A = std::complex<double>(p[i * 2 + 0], p[i * 2 + 1]);
        logP += AbsPriors_[i]->GetLogPrior(abs(A))
            + ArgPriors_[i]->GetLogPrior(yap::deg(arg(A)))
            - log(abs(A));      // jacobian
    }
    for (size_t i = FreeAmplitudes_.size() * 2; i < GetParameters().Size(); ++i)
        if (GetParameter(i).GetPrior())
            logP += GetParameter(i).GetLogPrior(p[i]);
    return logP;
}

//-------------------------
void bat_fit::CalculateObservables(const std::vector<double>& p)
{
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i) {
        auto A = std::complex<double>(p[i * 2 + 0], p[i * 2 + 1]);
        GetObservables()[i * 2 + 0] = abs(A);
        GetObservables()[i * 2 + 1] = yap::deg(arg(A));
    }
    // unsigned c = GetCurrentChain();
    // for (size_t i = 0; i < CalculatedFitFractions_[c].size(); ++i)
    //     GetObservables()[FreeAmplitudes_.size() * 2 + i] = CalculatedFitFractions_[c][i].value();
}

//-------------------------
void bat_fit::addParameter(std::string name, std::shared_ptr<yap::ComplexParameter> P, std::complex<double> low, std::complex<double> high, std::string latex, std::string units)
{
    if (std::find(Parameters_.begin(), Parameters_.end(), P) != Parameters_.end())
        throw yap::exceptions::Exception("trying to add parameter twice", "bat_fit::addParameter");
    Parameters_.push_back(P);
    AddParameter(name + "_re", real(low), real(high), latex.empty() ? latex : "Re(" + latex + ")", units);
    AddParameter(name + "_im", imag(low), imag(high), latex.empty() ? latex : "Im(" + latex + ")", units);
}

//-------------------------
void bat_fit::addParameter(std::string name, std::shared_ptr<yap::RealParameter> P, double low, double high, std::string latex, std::string units)
{
    if (std::find(Parameters_.begin(), Parameters_.end(), P) != Parameters_.end())
        throw yap::exceptions::Exception("trying to add parameter twice", "bat_fit::addParameter");
    Parameters_.push_back(P);
    AddParameter(name, low, high, latex, units);
}

//-------------------------
size_t bat_fit::findFreeAmplitude(std::shared_ptr<yap::FreeAmplitude> A) const
{
    auto it = std::find(FreeAmplitudes_.begin(), FreeAmplitudes_.end(), A);

    if (it == FreeAmplitudes_.end())
        throw yap::exceptions::Exception("FreeAmplitude not found", "setPrior");

    return (it - FreeAmplitudes_.begin()) * 2;
}

//-------------------------
void bat_fit::setPriors(std::shared_ptr<yap::FreeAmplitude> fa, BCPrior* amp_prior, BCPrior* arg_prior)
{
    if (!amp_prior)
        throw yap::exceptions::Exception("amp_prior is null", "bat_fit::setPrior");
    if (!arg_prior)
        throw yap::exceptions::Exception("phase_prior is null", "bat_fit::setPrior");
    
    auto i = findFreeAmplitude(fa);

    AbsPriors_[i / 2] = amp_prior;
    ArgPriors_[i / 2] = arg_prior;
}

//-------------------------
void bat_fit::setRealImagRanges(std::shared_ptr<yap::FreeAmplitude> fa, double real_low, double real_high, double imag_low, double imag_high)
{
    auto i = findFreeAmplitude(fa);
    GetParameter(i).SetLimits(real_low, real_high);
    GetParameter(i + 1).SetLimits(imag_low, imag_high);
}
    
//-------------------------
void bat_fit::setAbsArgRanges(std::shared_ptr<yap::FreeAmplitude> fa, double abs_low, double abs_high, double arg_low, double arg_high)
{
    auto i = findFreeAmplitude(fa);
    GetObservable(i).SetLimits(abs_low, abs_high);
    GetObservable(i + 1).SetLimits(arg_low, arg_high);
}

//-------------------------
void bat_fit::fix(std::shared_ptr<yap::FreeAmplitude> A, double amp, double phase)
{
    auto i = findFreeAmplitude(A);
    auto a = std::polar(amp, yap::rad(phase));
    GetParameter(i).Fix(real(a));
    GetParameter(i + 1).Fix(imag(a));
}

//-------------------------
void bat_fit::MCMCUserInitialize()
{
    bat_yap_base::MCMCUserInitialize();
    // if (!CalculatedFitFractions_.empty())
    // CalculatedFitFractions_.assign(GetNChains(), yap::RealIntegralElementVector(DecayTrees_.size()));
}

