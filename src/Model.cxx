#include "Model.h"

#include "Attributes.h"
#include "BlattWeisskopf.h"
#include "CalculationStatus.h"
#include "CompensatedSum.h"
#include "DataAccessor.h"
#include "DataPartition.h"
#include "DataPoint.h"
#include "DataSet.h"
#include "DecayChannel.h"
#include "DecayingState.h"
#include "DecayTree.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
#include "FreeAmplitude.h"
#include "Group.h"
#include "HelicityAngles.h"
#include "logging.h"
#include "MassAxes.h"
#include "Parameter.h"
#include "RecalculableDataAccessor.h"
#include "SpinAmplitudeCache.h"
#include "VariableStatus.h"

/// \todo Find better place for this
INITIALIZE_EASYLOGGINGPP

#include <future>

namespace yap {

//-------------------------
Model::Model(std::unique_ptr<SpinAmplitudeCache> SAC) :
    Locked_(false),
    FourMomenta_(std::make_shared<FourMomenta>(*this)),
    HelicityAngles_(*this)
{
    if (!SAC)
        throw exceptions::Exception("SpinAmplitudeCache unset", "Model::Model");
    if (!SAC->empty())
        throw exceptions::Exception("SpinAmplitudeCache not empty", "Model::Model");

    setCoordinateSystem({ThreeVector<double>({1., 0., 0.}),
                         ThreeVector<double>({0., 1., 0.}),
                         ThreeVector<double>({0., 0., 1.})});

    SAC->setModel(*this);
    SpinAmplitudeCache_ = std::move(SAC);
}

//-------------------------
void Model::calculate(DataPartition& D) const
{
    // update calculation statuses
    for (const auto& rda : RecalculableDataAccessors_)
        rda->updateCalculationStatus(D);

    // call calculate on all RecalculableDataAccessors
    for (const auto& rda : RecalculableDataAccessors_)
        rda->calculate(D);
}

//-------------------------
const double intensity(const InitialStateMap::value_type& isp_mix, const DataPoint& d)
{
    return std::accumulate(isp_mix.second.begin(), isp_mix.second.end(), 0.,
                           [&](double& I, const AdmixtureMap::value_type& m_b)
                           // += admixture factory * intensity of decay tree for spin projection
                           { return I += m_b.second->value() * intensity(isp_mix.first->decayTrees().at(m_b.first), d);});
}

//-------------------------
const double intensity(const InitialStateMap& isp_map, const DataPoint& d)
{
    return std::accumulate(isp_map.begin(), isp_map.end(), 0.,
                           [&](double& I, const InitialStateMap::value_type& isp_mix)
                           { return I += intensity(isp_mix, d); });
}

//-------------------------
// hidden helper function,
// resolves C++ problem related to naming of functions and call to std::async below
const double sum_of_logs_of_intensities(const Model& M, DataPartition& D, double ped)
{
    // calculate components
    M.calculate(D);

    // if pedestal is zero
    if (ped == 0)
        return std::accumulate(D.begin(), D.end(), CompensatedSum<double>(0.),
                               [&](CompensatedSum<double>& l, const DataPoint& d)
                               {return l += log(intensity(M.initialStates(), d));});
    // else
    return std::accumulate(D.begin(), D.end(), CompensatedSum<double>(0.),
                           [&](CompensatedSum<double>& l, const DataPoint& d)
                           {return l += (log(intensity(M.initialStates(), d)) - ped);});
}

//-------------------------
const double sum_of_log_intensity(const Model& M, DataPartition& D, double ped)
{
    if (M.initialStates().empty())
        throw exceptions::Exception("Model has no initial states", "sum_of_log_intensity");

    return sum_of_logs_of_intensities(M, D, ped);
}

//-------------------------
const double sum_of_log_intensity(const Model& M, DataPartitionVector& DP, double ped)
{
    // if DataPartitionVector is empty
    if (DP.empty())
        throw exceptions::Exception("DataPartitionVector is empty", "sum_of_log_intensity");

    // check no partitions are nullptr
    if (std::any_of(DP.begin(), DP.end(), std::logical_not<DataPartitionVector::value_type>()))
        throw exceptions::Exception("DataPartitionVector contains nullptr", "sum_of_log_intensity");

    // if threading is unnecessary
    if (DP.size() == 1)
        return sum_of_log_intensity(M, *DP[0], ped);

    if (M.initialStates().empty())
        throw exceptions::Exception("Model has no InitialStates", "sum_of_log_intensity");

    std::vector<std::future<double> > partial_sums;
    partial_sums.reserve(DP.size());

    // create thread for calculation on each partition
    for (auto& P : DP)
        // since std::async copies its arguments, even if they are supposed to be references, we need to use std::ref and std::cref
        partial_sums.push_back(std::async(std::launch::async, sum_of_logs_of_intensities, std::cref(M), std::ref(*P), ped));

    // wait for each partition to finish calculating
    return std::accumulate(partial_sums.begin(), partial_sums.end(), 0.,
                           [](double& l, std::future<double>& s) {return l += s.get();});
}

//-------------------------
bool Model::consistent() const
{
    bool C = true;

    C &= ParticleCombinationCache_.consistent();
    C &= SpinAmplitudeCache_->consistent();

    for (const auto& da : dataAccessors())
        C &= da->consistent();

    for (auto& p : initialStates())
        C &= p.first->consistent();

    return C;
}

//-------------------------
void Model::addParticleCombination(const ParticleCombination& pc)
{
    if (locked())
        throw exceptions::Exception("Model is locked and cannot be modified.", "Model::addParticleCombination");

    // if does not trace up to an ISP, halt
    if (!is_from_initial_state_particle_combination(pc, *this))
        return;

    FourMomenta_->addParticleCombination(pc);

    // call recursively on daughters
    for (auto& d : pc.daughters())
        addParticleCombination(*d);
}

//-------------------------
void Model::setFinalState(const FinalStateParticleVector& FSP)
{
    // check that FinalStateParticles_ is empty
    if (!FinalStateParticles_.empty())
        throw exceptions::Exception("Final-state particles already set", "Model::setFinalState");

    // check that none of the FSP's has yet been used
    // and that FinalStateParticles don't have ISP set to another ISP
    for (auto& fsp : FSP) {
        if (!fsp)
            throw exceptions::Exception("FinalStateParticle empty", "Model::setFinalState");
        if (!fsp->particleCombinations().empty())
            throw exceptions::Exception("FinalStateParticle already used", "Model::setFinalState");
        if (fsp->model() != nullptr)
            throw exceptions::Exception("FinalStateParticle already has Model set", "Model::setFinalState");
    }

    FinalStateParticles_.reserve(FSP.size());

    // set indices by order in vector
    for (auto& fsp : FSP) {
        fsp->addParticleCombination(*ParticleCombinationCache_.fsp(FinalStateParticles_.size()));
        fsp->setModel(this);
        FinalStateParticles_.push_back(fsp);
    }
}

//-------------------------
void Model::setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, StatusManager& sm) const
{
    fourMomenta()->setFinalStateMomenta(d, P, sm);

    // call calculate on all static data accessors in model
    for (const auto& sda : StaticDataAccessors_)
        sda->calculate(d, sm);
}

//-------------------------
AdmixtureMap admixture_map(const DecayTreeVectorMap& dtvm)
{
    AdmixtureMap M;
    // create new (spin projection -> real parameter) entry for each in dtvm
    std::transform(dtvm.begin(), dtvm.end(), std::inserter(M, M.end()),
                   [](const DecayTreeVectorMap::value_type& v)
                   {return AdmixtureMap::value_type(v.first, std::make_shared<RealParameter>(1.));});
    return M;
}

//-------------------------
void Model::addInitialState(std::shared_ptr<DecayingState> p)
{
    if (locked())
        throw exceptions::Exception("Model is locked and cannot be modified.", "Model::addInitialState");

    if (!p)
        throw exceptions::Exception("DecayingState is nullptr", "Model::addInitialState");
    
    // add DataAccessors to model / set model pointers
    p->registerWithModel();
    
    if (p->model() != this)
        throw exceptions::Exception("Initial state does not belong to this model", "Model::addInitialState");

    auto res = InitialStates_.emplace(std::dynamic_pointer_cast<DecayingState>(p->shared_from_this()), admixture_map(p->decayTrees()));

    if (res.second) { // new element was inserted
        for (auto& pc : p->particleCombinations())
            addParticleCombination(*pc);

    // no new element was inserted, check if insertion failed
    } else if (res.first == InitialStates_.end())
        throw exceptions::Exception("Failed to insert initial state", "Model::addInitialState");
}

//-------------------------
void Model::addInitialState(std::shared_ptr<DecayingState> p, double f)
{
    addInitialState(p);
    // retrieve admixture map
    auto it = InitialStates_.find(p);
    if (it == InitialStates_.end())
        throw exceptions::Exception("Initial state not found", "Model::addInitialState");
    if (it->second.size() != 1)
        throw exceptions::Exception("Several spin projections exist", "Model::addInitialState");
    *(it->second.begin()->second) = f;
}

//-------------------------
void fix_solitary_free_amplitudes(Model& m)
{
    // loop over sets of FreeAmplitudes grouped by parent state
    for (const auto fas : group(free_amplitudes(m), by_parent_state<>()))
        // if there is only one free amplitude for the decay of a state
        if (fas.size() == 1) {
            if (!*fas.begin())
                throw exceptions::Exception("FreeAmplitude is nullptr", "fix_solitary_free_amplitudes");
            // fix it
            (*fas.begin())->variableStatus() = VariableStatus::fixed;
        }
}

//-------------------------
std::string to_string(const AdmixtureMap& mix)
{
    return std::accumulate(mix.begin(), mix.end(), std::string(),
                           [](std::string& s, const AdmixtureMap::value_type& m_b)
                           {return s += "; b(m = " + std::to_string(m_b.first) + ") = " + to_string(*m_b.second);}).erase(0, 2);
}

//-------------------------
size_t all_fixed(const AdmixtureMap& mix)
{
    return std::all_of(mix.begin(), mix.end(),
                       [](const AdmixtureMap::value_type& m_b)
                       {return m_b.second->variableStatus() == VariableStatus::fixed;});
}

//-------------------------
std::vector<std::shared_ptr<DecayingState> > full_final_state_isp(const Model& M)
{
    std::vector<std::shared_ptr<DecayingState> > isps;
    isps.reserve(M.initialStates().size());
    // collect first those particles for whom all admixture factors are fixed
    for (const auto& isp_am : M.initialStates())
        if (decays_to_full_final_state(*isp_am.first) and all_fixed(isp_am.second))
            isps.push_back(isp_am.first);
    // then collect the rest
    for (const auto& isp_am : M.initialStates())
        if (decays_to_full_final_state(*isp_am.first) and !all_fixed(isp_am.second))
            isps.push_back(isp_am.first);
    return isps;
}

//-------------------------
void Model::setCoordinateSystem(const CoordinateSystem<double, 3>& cs)
{
    if (!is_right_handed(cs))
        throw exceptions::Exception("Coordinate system not right-handed", "Model::setCoordinateSystem");

    CoordinateSystem_ = unit(cs);
}

//-------------------------
void Model::lock()
{
    // if already locked, do nothing
    if (Locked_)
        return;

    // prune initial state particles
    for (auto& isp_mix : InitialStates_)
        isp_mix.first->pruneParticleCombinations();

    // fix any amplitudes that are the only for a particular
    // particle's decay
    fix_solitary_free_amplitudes(*this);

    // if only one ISP with only one spin projection, fix it's admixture
    if (InitialStates_.size() == 1 and InitialStates_.begin()->second.size() == 1)
        InitialStates_.begin()->second.begin()->second->variableStatus() = VariableStatus::fixed;

    // remove expired elements of DataAccessors_
    remove_expired(DataAccessors_);
    remove_expired(StaticDataAccessors_);

    // prune remaining DataAccessor's
    for (auto& D : DataAccessors_)
        D->pruneSymmetrizationIndices();

    // remove data accessors from list that don't need storage
    for (auto it = DataAccessors_.begin(); it != DataAccessors_.end(); ) {
        if (!(*it)->requiresStorage())
            it = DataAccessors_.erase(it);
        else
            ++it;
    }

    // set DataAccessor indices
    int index = -1;
    for (const auto& da : DataAccessors_)
        da->setIndex(++index);

    Locked_ = true;
}

//-------------------------
const MassAxes Model::massAxes(std::vector<std::vector<unsigned> > pcs)
{
    // if no axes requested, build default:
    if (pcs.empty()) {

        if (finalStateParticles().size() > 4)
            throw exceptions::Exception("Currently only supports final states of 4 or fewer particles for default axes", "Model::massAxes");

        // builds vector down first off diagonal, then second off-diagonal, etc
        for (unsigned i = 1; i < finalStateParticles().size() - 1; ++i)
            for (unsigned j = 0 ; i + j < finalStateParticles().size(); ++j)
                pcs.push_back({j, j + i});
    }

    unsigned n_fsp = finalStateParticles().size();
    unsigned n_axes = 3 * n_fsp - 7;

    // check that number of requested axes == n_axes
    if (pcs.size() != n_axes) {
        if (pcs.size() < n_axes)
            throw exceptions::Exception("too few axes requested ( " + std::to_string(pcs.size()) + " < " + std::to_string(n_axes) + " )",
                                        "Model::massAxes");
        else
            throw exceptions::Exception("too many axes requested ( " + std::to_string(pcs.size()) + " > " + std::to_string(n_axes) + " )",
                                        "Model::massAxes");
    }

    // for the moment, we only support 2-particle axes
    // check that all axes are 2 -particle
    if (std::any_of(pcs.begin(), pcs.end(), [](const std::vector<unsigned>& v) {return v.size() != 2;}))
        throw exceptions::Exception("only 2-particle axes supported currently", "Model::massAxes");

    ParticleCombinationVector M;

    for (auto& v : pcs) {

        // check that all indices are in range
        {
            auto it = std::find_if(v.end(), v.end(), std::bind(std::greater_equal<unsigned>(), std::placeholders::_1, n_fsp));
            if (it != v.end())
                throw exceptions::Exception("particle index out of range (" + std::to_string(*it) + " >= " + std::to_string(n_fsp) + ")",
                                            "Model::massAxes");
        }

        // check for duplicates
        {
            auto it = std::find_if(v.begin(), v.end(), [&](unsigned i) {return std::count(v.begin(), v.end(), i) != 1;});
            if (it != v.end())
                throw exceptions::Exception("duplicate index (" + std::to_string(*it) + ")given", "Model::massAxes");
        }

        // get fsp ParticleCombinations
        ParticleCombinationVector pcv;
        pcv.reserve(v.size());
        std::transform(v.begin(), v.end(), std::back_inserter(pcv), [&](unsigned i) {return particleCombinationCache().fsp(i);});
        auto pc = particleCombinationCache().composite(pcv);

        // check that pc isn't already in M
        if (std::any_of(M.begin(), M.end(), std::bind(&equal_by_orderless_content, pc, std::placeholders::_1)))
            throw exceptions::Exception("axis requested twice", "Model::massAxes");

        M.push_back(pc);
    }

    return MassAxes(M);
}

//-------------------------
DataSet Model::createDataSet(size_t n)
{
    if (!locked())
        lock();

    if (!locked())
        throw exceptions::Exception("data sets cannot be generated from an unlocked model.", "Model::createDataSet");

    // create empty data set
    DataSet D(*this);

    D.addEmptyDataPoints(n);

    return D;
}

//-------------------------
FreeAmplitudeSet free_amplitudes(const Model& M)
{
    FreeAmplitudeSet S;
    for (const auto& ds_am : M.initialStates()) {
        auto s = free_amplitudes(*ds_am.first);
        S.insert(s.begin(), s.end());
    }
    return S;
}

//-------------------------
ParticleSet particles(const Model& M)
{
    ParticleSet S;
    for (const auto& ds_am : M.initialStates()) {
        auto s = particles(*ds_am.first);
        S.insert(s.begin(), s.end());
    }
    return S;
}

//-------------------------
void Model::setParameterFlagsToUnchanged()
{
    for (auto& d : RecalculableDataAccessors_)
        d->setParameterFlagsToUnchanged();
}

}
