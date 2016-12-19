#include "fit_base.h"

#include <DataSet.h>
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

#include <algorithm>

// -----------------------
fit_base::fit_base(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs)
    : bat_yap_base(name, std::move(M)),
      FitData_(model()->createDataSet()),
      FitPartitions_(1, &FitData_),
      IntegralData_(model()->createDataSet()),
      IntegralPartitions_(1, &IntegralData_),
      NIntegrationPoints_(0),
      NIntegrationPointsBatchSize_(0),
      NIntegrationThreads_(1),
      Integral_(*model())
{
    // create mass axes
    axes() = model()->massAxes(pcs);
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
        std::string indices_string = parname.substr(parname.rfind("_") + 1);
        if (indices_string.size() != 2)
            throw yap::exceptions::Exception("parameter name \"" + parname + "\" does not contain two one-digit indices", "find_mass_axes");
        if (!std::all_of(indices_string.begin(), indices_string.end(), isdigit))
            throw yap::exceptions::Exception("parameter name \"" + parname + "\" does not contain two one-digit indices", "find_mass_axes");
        
        // read indices:
        std::vector<unsigned> indices;
        indices.reserve(2);
        std::transform(indices_string.begin(), indices_string.end(), std::back_inserter(indices), [](char c){return c - '0';});
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
size_t load_data(yap::DataSet& data, const yap::Model& M, const yap::MassAxes& A, double initial_mass, TTree& t_mcmc, int N, int lag)
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
    unsigned Chain;
    t_mcmc.SetBranchAddress("Chain", &Chain);

    unsigned long long n_entries = t_mcmc.GetEntries();


    if (N < 0)
        // attempt to load all data
        N = n_entries;

    if (lag < 0)
        // calculate lag
        lag = n_entries / N;
    lag = std::max(lag, 1);

    int n_attempted = 0;
    size_t old_size = data.size();

    for (unsigned long long n = 0; n < n_entries and n_attempted < N; ++n) {
        t_mcmc.GetEntry(n);

        if (Phase <= 0)
            continue;

        if (Iteration % lag != 0)
            continue;

        // if (fabs(m2[0] - 1.35 * 1.35) > 0.1 or m2[1] > 1.55 or m2[1] < 0.58)
        //     continue;

        ++n_attempted;

        auto P = calculate_four_momenta(initial_mass, M, A, m2);
        if (P.empty())
            std::cout << "point is out of phase space!";
        data.push_back(P);
    }

    if (data.empty())
        LOG(INFO) << "No data loaded.";
    else {
        LOG(INFO) << "Loaded " << data.size() - old_size << " data points (" << ((data.size() - old_size) * data[0].bytes() * 1.e-6) << " MB)"
                << " from a tree of size " << n_entries << ", with a lag of " << lag;
        if (old_size != 0)
            LOG(INFO) << "Total data size now " << data.size() << " points (" << (data.bytes() * 1.e-6) << " MB)";
    }

    if (int(data.size() - old_size) < N)
        LOG(WARNING) << "could not load as many data points as requested. Reduce the lag (or set it to -1 to automatically determine the lag).";

    return data.size() - old_size;
}

//-------------------------
void fit_base::setParameters(const std::vector<double>& p)
{
    integrate();
}

//-------------------------
void fit_base::integrate()
{
    if (IntegrationPointGenerator_)
        yap::ImportanceSampler::calculate(Integral_, IntegrationPointGenerator_, NIntegrationPoints_, NIntegrationPointsBatchSize_, NIntegrationThreads_);
    else
        yap::ImportanceSampler::calculate(Integral_, IntegralPartitions_);
}

// ---------------------------------------------------------
double fit_base::LogLikelihood(const std::vector<double>& p)
{
    setParameters(p);
    double L = sum_of_log_intensity(*model(), FitPartitions_, log(integral(Integral_).value()));
    model()->setParameterFlagsToUnchanged();
    increaseLikelihoodCalls();
    return L;
}
