#include "bat_fit.h"

#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <VariableStatus.h>

#include <BAT/BCPrior.h>

#include <TTree.h>

void unambiguous_importance_sampler_calculate(yap::ModelIntegral& M, yap::DataPartitionVector& D)
{ yap::ImportanceSampler::calculate(M, D); }

// -----------------------
bat_fit::bat_fit(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs)
    : bat_yap_base(name, std::move(M)),
      FitData_(model()->createDataSet()),
      FitPartitions_(1, &FitData_),
      IntegralData_(model()->createDataSet()),
      IntegralPartitions_(1, &IntegralData_),
      Integrator_(integrator_type(unambiguous_importance_sampler_calculate)),
      Integral_(*model()),
      FirstParameter_(-1)
{
    // create mass axes
    axes() = model()->massAxes(pcs);

    // loop over all free amplitudes
    for (const auto& fa : free_amplitudes(*model())) {
        // ignore fixed free amplitudes
        if (fa->variableStatus() == yap::VariableStatus::fixed)
            continue;
        // add amplitude parameter
        AddParameter("amp(" + to_string(*fa->decayChannel()) + " M = " + yap::spin_to_string(fa->twoM()) + ")",
                     1.e-3, 2 * norm(fa->value()));
        /// add phase parameter
        AddParameter("phase(" + to_string(*fa->decayChannel()) + " M = " + yap::spin_to_string(fa->twoM()) + ")",
                     -180, 180);
        /// add free amplitude to list
        FreeAmplitudes_.push_back(fa);
    }

    SetPriorConstantAll();

    // add observables for all fit fractions
    for (const auto& b2_dtvi : Integral_.integrals())
        for (const auto& dt : b2_dtvi.second.decayTrees()) {
            DecayTrees_.push_back(dt);
            AddObservable("fit_frac(" + to_string(*dt->freeAmplitude()->decayChannel()) + " M = " + yap::spin_to_string(dt->freeAmplitude()->twoM()) + ")", 0, 1.1);
        }
    CalculatedFitFractions_.assign(1, yap::RealIntegralElementVector(DecayTrees_.size()));
}

//-------------------------
std::vector<std::vector<unsigned> > find_mass_axes(TTree& t_pars)
{
    std::vector<std::vector<unsigned> > pcs;
    pcs.reserve(t_pars.GetEntries());

    bool parameter;
    char c_parname[10];
    t_pars.SetBranchAddress("parameter", &parameter);
    t_pars.SetBranchAddress("name", &c_parname);
    for (unsigned n = 0; n < t_pars.GetEntries(); ++n) {
        t_pars.GetEntry(n);
        if (!parameter)
            continue;
        std::string parname(c_parname);
        // parameter names should be m2_ij
        if (parname.find("m2_") != 0)
            throw yap::exceptions::Exception("parameter name \"" + parname + "\" does not match \"m2_ij\"",
                                             "find_mass_axes");
        if (parname.size() != 5)
            throw yap::exceptions::Exception("parameter name \"" + parname + "\" is not right length "
                                             + "(" + std::to_string(parname.size()) + " != 5)",
                                             "find_mass_axes");
        // read indices:
        std::vector<unsigned> indices;
        indices.reserve(2);
        std::transform(parname.begin() + parname.rfind("_") + 1, parname.end(), std::back_inserter(indices), [](char c) {return std::atoi(&c);});
        // for (size_t i = parname.rfind("_") + 1; i < parname.size(); ++i)
        //     indices.push_back(std::stoi(parname.substr(i, 1)));
        pcs.push_back(indices);
    }

    return pcs;
}

//-------------------------
void set_address(const yap::MassAxes::value_type& a,
                 std::vector<double>& m2,
                 TTree& t_mcmc)
{
    m2.push_back(0);
    t_mcmc.SetBranchAddress(indices_string(*a, "m2_", "").data(), &m2.back());
}

//-------------------------
size_t load_data(yap::DataSet& data, const yap::Model& M, const yap::MassAxes& A, double initial_mass, TTree& t_mcmc, int N, unsigned lag)
{
    if (A.empty())
        throw yap::exceptions::Exception("mass axes empty", "load_data");

    // set branch addresses
    std::vector<double> m2;
    m2.reserve(A.size());
    std::for_each(A.begin(), A.end(), std::bind(set_address, std::placeholders::_1, std::ref(m2), std::ref(t_mcmc)));
    if (m2.size() != A.size())
        throw yap::exceptions::Exception("not all mass axes loaded from TTree", "load_data");

    //
    // load data
    //
    int Phase = -1;
    t_mcmc.SetBranchAddress("Phase", &Phase);
    unsigned Iteration;
    t_mcmc.SetBranchAddress("Iteration", &Iteration);

    int n_attempted = 0;
    size_t old_size = data.size();

    for (long long n = 0; n < t_mcmc.GetEntries() and (N < 0 or n_attempted < N); ++n) {
        t_mcmc.GetEntry(n);

        if (Phase <= 0)
            continue;

        if (Iteration % lag != 0)
            continue;

        if (abs(m2[0] - 1.35 * 1.35) > 0.1 or m2[1] > 1.55 or m2[1] < 0.58)
            continue;

        ++n_attempted;

        auto P = calculate_four_momenta(initial_mass, M, A, m2);
        if (P.empty())
            std::cout << "point is out of phase space!";
        data.push_back(P);
    }

    if (data.empty())
        LOG(INFO) << "No data loaded.";
    else {
        LOG(INFO) << "Loaded " << data.size() - old_size << " data points (" << ((data.size() - old_size) * data[0].bytes() * 1.e-6) << " MB)";
        if (old_size != 0)
            LOG(INFO) << "Total data size now " << data.size() << " points (" << (data.bytes() * 1.e-6) << " MB)";
    }
    return data.size() - old_size;
}

//-------------------------
void bat_fit::setParameters(const std::vector<double>& p)
{
    for (size_t i = 0; i < FreeAmplitudes_.size(); ++i)
        *FreeAmplitudes_[i] = std::polar(p[i * 2], yap::rad(p[i * 2 + 1]));

    yap::set_values(Parameters_.begin(), Parameters_.end(), p.begin() + FirstParameter_, p.end());

    Integrator_(Integral_, IntegralPartitions_);

    // calculate fit fractions
    unsigned c = GetCurrentChain();
    size_t i = 0;
    for (const auto& b2_dtvi : Integral_.integrals()) {
        auto ff = fit_fractions(b2_dtvi.second);
        for (const auto& f : ff)
            CalculatedFitFractions_[c][i++] = f;
    }
}

// ---------------------------------------------------------
double bat_fit::LogLikelihood(const std::vector<double>& p)
{
    setParameters(p);
    double L = sum_of_log_intensity(*model(), FitPartitions_, log(integral(Integral_).value()));
    model()->setParameterFlagsToUnchanged();
    increaseLikelihoodCalls();
    return L;
}

//-------------------------
void bat_fit::CalculateObservables(const std::vector<double>& p)
{
    unsigned c = GetCurrentChain();
    for (size_t i = 0; i < CalculatedFitFractions_[c].size(); ++i)
        GetObservables()[i] = CalculatedFitFractions_[c][i].value();
}

//-------------------------
void bat_fit::addParameter(std::string name, std::shared_ptr<yap::ComplexParameter> P, std::complex<double> low, std::complex<double> high, std::string latex, std::string units)
{
    if (std::find(Parameters_.begin(), Parameters_.end(), P) != Parameters_.end())
        throw yap::exceptions::Exception("trying to add parameter twice", "bat_fit::addParameter");
    Parameters_.push_back(P);
    if (FirstParameter_ < 0)
        FirstParameter_ = GetParameters().Size();
    AddParameter(name + "_re", real(low), real(high), latex.empty() ? latex : "Re(" + latex + ")", units);
    AddParameter(name + "_im", imag(low), imag(high), latex.empty() ? latex : "Im(" + latex + ")", units);
}

//-------------------------
void bat_fit::addParameter(std::string name, std::shared_ptr<yap::RealParameter> P, double low, double high, std::string latex, std::string units)
{
    if (std::find(Parameters_.begin(), Parameters_.end(), P) != Parameters_.end())
        throw yap::exceptions::Exception("trying to add parameter twice", "bat_fit::addParameter");
    Parameters_.push_back(P);
    if (FirstParameter_ < 0)
        FirstParameter_ = GetParameters().Size();
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
void bat_fit::setPrior(std::shared_ptr<yap::FreeAmplitude> A, BCPrior* amp_prior, BCPrior* phase_prior)
{
    auto i = findFreeAmplitude(A);

    // add amp
    if (amp_prior) {
        double low = amp_prior->GetMean() - 3 * amp_prior->GetStandardDeviation();
        GetParameters()[i].SetLimits(std::max(low, 0.), amp_prior->GetMean() + 3 * amp_prior->GetStandardDeviation());
        GetParameters()[i].SetPrior(amp_prior);
    } else
        throw yap::exceptions::Exception("amp_prior is null", "bat_fit::setPrior");

    if (phase_prior) {
        GetParameters()[i + 1].SetLimits(phase_prior->GetMean() - 3. * phase_prior->GetStandardDeviation(),
                                         phase_prior->GetMean() + 3. * phase_prior->GetStandardDeviation());
        GetParameters()[i + 1].SetPrior(phase_prior);
    } else {
        throw yap::exceptions::Exception("phase_prior is null", "bat_fit::setPrior");
    }
}

//-------------------------
void bat_fit::setPrior(std::shared_ptr<yap::FreeAmplitude> A,  double amp_low, double amp_high, double phase_low, double phase_high)
{
    auto i = findFreeAmplitude(A);

    GetParameters()[i].SetLimits(amp_low, amp_high);
    GetParameters()[i].SetPriorConstant();

    GetParameters()[i + 1].SetLimits(phase_low, phase_high);
    GetParameters()[i + 1].SetPriorConstant();
}

//-------------------------
void bat_fit::fix(std::shared_ptr<yap::FreeAmplitude> A, double amp, double phase)
{
    auto i = findFreeAmplitude(A);
    GetParameters()[i].SetLimits(0, 2 * amp);
    GetParameters()[i].Fix(amp);
    GetParameters()[i + 1].SetLimits(phase - 0.1, phase + 0.1);
    GetParameters()[i + 1].Fix(phase);
}

//-------------------------
void bat_fit::MCMCUserInitialize()
{
    bat_yap_base::MCMCUserInitialize();
    CalculatedFitFractions_.assign(GetNChains(), CalculatedFitFractions_[0]);
}

