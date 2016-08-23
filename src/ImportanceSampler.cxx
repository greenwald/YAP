#include "ImportanceSampler.h"

#include "DataPartition.h"
#include "DecayTree.h"
#include "DecayTreeVectorIntegral.h"
#include "Exceptions.h"
#include "IntegralElement.h"
#include "Model.h"
#include "ModelIntegral.h"

#include <future>
#include <vector>

namespace yap {

//-------------------------
unsigned ImportanceSampler::partially_calculate_subIntegral(DecayTreeVectorIntegral& I, DataPartition& D)
{
    // create holder for amplitudes
    std::vector<std::complex<double> > A(I.decayTrees().size());

    // counter for number of data points
    unsigned n = 0;

    // loop over data points
    for (const auto& d : D) {

        // calculate the amplitudes of all decay trees
        std::transform(I.decayTrees().begin(), I.decayTrees().end(), A.begin(),
        [&d](const std::shared_ptr<DecayTree>& dt) {return dt->dataDependentAmplitude(d);});

        // increase number of points
        ++n;

        for (size_t i = 0; i < A.size(); ++i) {

            // calculate difference from mean
            double delta_diag = norm(A[i]) - I.diagonals()[i].value();
            // update mean
            diagonals(I)[i].value() += delta_diag / n;

            for (size_t j = i + 1; j < A.size(); ++j) {
                // calculate difference from mean
                auto delta_offdiag = conj(A[i]) * A[j] - I.offDiagonals()[i][j - i - 1].value();
                // update mean
                offDiagonals(I)[i][j - i - 1].value() += delta_offdiag / static_cast<double>(n);
            }
        }
    }

    // return number of data points evaluated
    return n;
}

//-------------------------
unsigned ImportanceSampler::partially_calculate(std::vector<DecayTreeVectorIntegral*>& J, DataPartition& D)
{
    if (J.empty())
        throw exceptions::Exception("vector is empty", "ImportanceSampler::partialCalculate");

    if (!J[0]->model())
        throw exceptions::Exception("Model is nullptr", "ImportanceSampler::partialCalculate");

    // calculate on data partition
    J[0]->model()->calculate(D);

    unsigned n = 0;
    for (auto& j : J)
        n = ImportanceSampler::partially_calculate_subIntegral(*j, D);

    return n;
}

//-------------------------
std::vector<DecayTreeVectorIntegral*> ImportanceSampler::select_changed(ModelIntegral& I)
{
    std::vector<DecayTreeVectorIntegral*> C;
    C.reserve(integrals(I).size());

    for (auto& b_I : integrals(I))
        if (std::any_of(b_I.second.decayTrees().begin(), b_I.second.decayTrees().end(), &has_changed))
            C.push_back(&b_I.second);

    return C;
}

//-------------------------
void ImportanceSampler::calculate(ModelIntegral& I, DataPartition& D)
{
    // get DecayTreeVectorIntegral's for DecayTree's that need to be calculated
    auto J = select_changed(I);

    // if nothing requires recalculation, return
    if (J.empty())
        return;

    // reset those to be recalculated
    for (auto& j : J)
        reset(*j);

    // calculate it
    partially_calculate(J, D);

}

//-------------------------
void ImportanceSampler::calculate(ModelIntegral& I, DataPartitionVector& DPV)
{
    // if only one data partition, don't thread:
    if (DPV.size() == 1) {
        calculate(I, *DPV[0]);
        return;
    }

    // get DecayTreeVectorIntegral's for DecayTree's that need to be calculated
    auto J = select_changed(I);

    // if nothing requires recalculation, return
    if (J.empty())
        return;

    // reset those to be recalculated
    for (auto& j : J)
        reset(*j);

    // create copies for running with in each partition
    std::vector<std::vector<DecayTreeVectorIntegral*> > m(DPV.size());
    for (auto& m_j : m) {
        m_j.reserve(J.size());
        for (const auto& j : J)
            m_j.push_back(new DecayTreeVectorIntegral(*j));
    }

    // run over each partition storing number of events used in each calculation
    std::vector<std::future<unsigned> > n;
    n.reserve(m.size());
    // create thread for each partial calculation
    for (size_t i = 0; i < m.size(); ++i)
        n.push_back(std::async(std::launch::async, &ImportanceSampler::partially_calculate,
                               std::ref(m[i]), std::ref(*DPV[i])));

    // calculate data fractions:
    // also waits for threads to finish calculating
    std::vector<double> f;
    f.reserve(n.size());
    std::transform(n.begin(), n.end(), std::back_inserter(f), std::mem_fn(&std::future<unsigned>::get));
    double N = std::accumulate(f.begin(), f.end(), 0.);
    std::transform(f.begin(), f.end(), f.begin(), std::bind(std::divides<double>(), std::placeholders::_1, N));

    for (size_t i = 0; i < m.size(); ++i)
        for (size_t j = 0; j < J.size(); ++j)
            *J[j] += (*m[i][j] *= f[i]);

}

}
