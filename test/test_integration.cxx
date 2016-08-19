#include <catch.hpp>
#include <catch_capprox.hpp>

#include "helperFunctions.h"

#include "BreitWigner.h"
#include "CompensatedSum.h"
#include "DataPartition.h"
#include "DataPoint.h"
#include "DataSet.h"
#include "FinalStateParticle.h"
#include "HelicityFormalism.h"
#include "ImportanceSampler.h"
#include "logging.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "Model.h"
#include "ModelIntegral.h"
#include "Parameter.h"
#include "ParticleFactory.h"
#include "PDL.h"
#include "PHSP.h"
#include "Resonance.h"

#include <future>
#include <memory>

/**
 *  Test the integration
 */

namespace yap {

//-------------------------
// hidden helper function,
// resolves C++ problem related to naming of functions and call to std::async below
const double sum_of_intensities(const Model& M, DataPartition& D, double ped)
{
    // calculate components
    M.calculate(D);

    // sum intensities over data points in partition
    // if pedestal is zero
    if (ped == 0)
        return std::accumulate(D.begin(), D.end(), CompensatedSum<double>(0.),
                               [&](CompensatedSum<double>& l, const DataPoint & d)
                               {return l += intensity(M.initialStateParticles(), d);});
    // else
    return std::accumulate(D.begin(), D.end(), CompensatedSum<double>(0.),
                           [&](CompensatedSum<double>& l, const DataPoint & d)
                           {return l += intensity(M.initialStateParticles(), d) - ped;});
}

//-------------------------
const double sum_of_intensity(const Model& M, DataPartition& D, double ped)
{
    if (M.initialStateParticles().empty())
        throw exceptions::Exception("Model has no initialStateParticles", "sum_of_intensity");

    return sum_of_intensities(M, D, ped);
}

//-------------------------
const double sum_of_intensity(const Model& M, DataPartitionVector& DP, double ped)
{
    // if DataPartitionVector is empty
    if (DP.empty())
        throw exceptions::Exception("DataPartitionVector is empty", "sum_of_intensity");

    // check no partitions are nullptr
    if (std::any_of(DP.begin(), DP.end(), std::logical_not<DataPartitionVector::value_type>()))
        throw exceptions::Exception("DataPartitionVector contains nullptr", "sum_of_intensity");

    // if threading is unnecessary
    if (DP.size() == 1)
        return sum_of_intensity(M, *DP[0], ped);

    if (M.initialStateParticles().empty())
        throw exceptions::Exception("Model has no InitialStateParticles", "sum_of_intensity");

    std::vector<std::future<double> > partial_sums;
    partial_sums.reserve(DP.size());

    // create thread for calculation on each partition
    for (auto& P : DP)
        // since std::async copies its arguments, even if they are supposed to be references, we need to use std::ref and std::cref
        partial_sums.push_back(std::async(std::launch::async, sum_of_intensities, std::cref(M), std::ref(*P), ped));

    // wait for each partition to finish calculating
    return std::accumulate(partial_sums.begin(), partial_sums.end(), 0.,
                           [](double & l, std::future<double>& s) {return l += s.get();});
}

}

TEST_CASE("integration")
{
    // disable debug logs in test
    //yap::disableLogs(el::Level::Debug);
    yap::plainLogs(el::Level::Debug);

    const unsigned nPoints = 1000;

    auto M = create_model();
    auto data = generate_data(*M, nPoints);
    auto partitions = yap::DataPartitionBlock::create(data, 4);

    double bruteForceIntegral = yap::sum_of_intensity(*M, partitions, 0) / nPoints;
    DEBUG("bruteForceIntegral = " << bruteForceIntegral);

    yap::ModelIntegral mi(*M);

    yap::ImportanceSampler::calculate(mi, partitions);
    double smartIntegral = integral(mi).value();
    DEBUG("smartIntegral = " << smartIntegral);

    REQUIRE(bruteForceIntegral == Approx(smartIntegral));
}

