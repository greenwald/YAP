#ifndef __BAT__FIT_BASE__H
#define __BAT__FIT_BASE__H

#include "bat_yap_base.h"

#include <fwd/DataPartition.h>
#include <fwd/FreeAmplitude.h>
#include <fwd/IntegralElement.h>
#include <fwd/Model.h>
#include <fwd/Parameter.h>

#include <DataSet.h>
#include <FourVector.h>
#include <ModelIntegral.h>

#include <memory>
#include <functional>
#include <string>
#include <random>
#include <vector>

class TTree;

class fit_base : public bat_yap_base
{

public:

    /// constructor
    /// \param name name of bat model
    /// \param M yap::model
    /// \param pcs vector<vector<unsigned>> defining mass axes
    fit_base(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs);

    /// log likelihood
    virtual double LogLikelihood(const std::vector<double>& p) override;

    /// \return FitData_
    yap::DataSet& fitData()
    { return FitData_; }

    /// \return FitPartitions_
    yap::DataPartitionVector& fitPartitions()
    { return FitPartitions_; }

    /// \return IntegralData_
    yap::DataSet& integralData()
    { return IntegralData_; }

    /// \return IntegralPartitions_
    yap::DataPartitionVector& integralPartitions()
    { return IntegralPartitions_; }

    yap::ModelIntegral& modelIntegral()
    { return Integral_; }

    /// set parameters of integration
    /// \param N number of points to use for integration
    /// \param n batch size for integration
    /// \param t number of threads
    void setNIntegrationPoints(unsigned N, unsigned n, unsigned t = 1)
    { NIntegrationPoints_ = N; NIntegrationPointsBatchSize_ = n; NIntegrationThreads_ = t; }

    /// \typedef Generator
    /// function for generating new points for integration
    using Generator = std::function<std::vector<yap::FourVector<double> >()>;

    /// \return IntegrationPointGenerator_
    Generator& integrationPointGenerator()
    { return IntegrationPointGenerator_; }

    /// set parameters into model
    virtual void setParameters(const std::vector<double>& p) = 0;

    /// perform the integration
    virtual void integrate();

protected:

    /// DataSet to fit the model to
    yap::DataSet FitData_;

    /// Partitioning of FitData_
    yap::DataPartitionVector FitPartitions_;

    /// DataSet to fit the model to
    yap::DataSet IntegralData_;

    /// Partitioning of FitData_
    yap::DataPartitionVector IntegralPartitions_;

    /// Number of points to integrate with
    unsigned NIntegrationPoints_;

    /// Batch size for generating integration points
    unsigned NIntegrationPointsBatchSize_;

    /// Number of threads for integration
    unsigned NIntegrationThreads_;

    /// generator for integration
    Generator IntegrationPointGenerator_;

    /// stores integral result
    yap::ModelIntegral Integral_;

    /// list of decay trees integrated over
    yap::DecayTreeVector DecayTrees_;

};

/// load data from a TTree into a DataSet
/// \param data DataSet to load into
/// \param M Model to load with
/// \param A MassAxes to load with
/// \param t_mcmc TTree to load from
/// \param N max number of data points to (attempt to) load
/// \param lag Lag to apply to iterations when reading from TTree
///        per default, data points will be selected uniformly from t_mcmc
/// \param eps Amount to smear momenta by
size_t load_data(yap::DataSet& data, const yap::Model& M,
                 const yap::MassAxes& A, double initial_mass, TTree& t_mcmc,
                 int N = -1, int lag = -1);

/// find mass axes from TTree of parameters
std::vector<std::vector<unsigned> > find_mass_axes(TTree& t_pars);
#endif
