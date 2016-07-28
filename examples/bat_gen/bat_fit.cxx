#include "bat_fit.h"

#include <DataSet.h>
#include <DecayingParticle.h>
#include <FourVector.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <MassAxes.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>
#include <ParticleCombination.h>

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
      Integral_(*model())
{
    // create mass axes
    axes() = model()->massAxes(pcs);
}

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

    long long n = 0;
    // find first entry of main run
    while (Phase <= 0 and n < t_mcmc.GetEntries())
        t_mcmc.GetEntry(n++);

    int n_attempted = 0;
    size_t old_size = data.size();

    for (; n < t_mcmc.GetEntries() and (N < 0 or n_attempted < N); ++n) {
        t_mcmc.GetEntry(n);

        if (Phase <= 0)
            continue;

        if (Iteration % lag != 0)
            continue;

        ++n_attempted;

        auto P = M.calculateFourMomenta(A, m2, initial_mass);
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

// ---------------------------------------------------------
double bat_fit::LogLikelihood(const std::vector<double>& p)
{
    yap::set_values(Parameters_, p);
    Integrator_(Integral_, IntegralPartitions_);
    double L = sum_of_log_intensity(*model(), FitPartitions_, log(integral(Integral_).value()));
    model()->setParameterFlagsToUnchanged();
    increaseLikelihoodCalls();
    return L;
}

//-------------------------
void bat_fit::addParameter(std::string name, std::shared_ptr<yap::ComplexParameter> P, std::complex<double> low, std::complex<double> high, std::string latex, std::string units)
{
    if (std::find(Parameters_.begin(), Parameters_.end(), P) != Parameters_.end())
        throw yap::exceptions::Exception("trying to add parameter twice", "bat_fit::addParameter");
    Parameters_.push_back(P);
    AddParameter(name + "_re", real(low), real(high), "Re(" + latex + ")", units);
    AddParameter(name + "_im", imag(low), imag(high), "Im(" + latex + ")", units);
}

//-------------------------
void bat_fit::addParameter(std::string name, std::shared_ptr<yap::RealParameter> P, double low, double high, std::string latex, std::string units)
{
    if (std::find(Parameters_.begin(), Parameters_.end(), P) != Parameters_.end())
        throw yap::exceptions::Exception("trying to add parameter twice", "bat_fit::addParameter");
    Parameters_.push_back(P);
    AddParameter(name, low, high, latex, units);
}