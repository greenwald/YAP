/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/// \file

#ifndef yap_Model_h
#define yap_Model_h

#include "fwd/DataAccessor.h"
#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/DataSet.h"
#include "fwd/DecayingParticle.h"
#include "fwd/Filters.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/FourMomenta.h"
#include "fwd/FourVector.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/HelicityAngles.h"
#include "fwd/MassAxes.h"
#include "fwd/MeasuredBreakupMomenta.h"
#include "fwd/Parameter.h"
#include "fwd/Particle.h"
#include "fwd/RecalculableDataAccessor.h"
#include "fwd/StaticDataAccessor.h"
#include "fwd/StatusManager.h"

#include "CoordinateSystem.h"
#include "ParticleCombinationCache.h"
#include "SpinAmplitudeCache.h"

#include <complex>
#include <memory>
#include <vector>

namespace yap {

/// \class Model
/// \brief Class implementing a PWA model
/// \author Johannes Rauch, Daniel Greenwald
class Model
{
public:

    /// Constructor
    /// \param SAC unique_ptr to SpinAmplitudeCache
    Model(std::unique_ptr<SpinAmplitudeCache> SAC);

    /// copy constructor deleted
    /// \todo Implement deep copy.
    Model(const Model&) = delete;

    /// copy assignment operator
    /// \todo Implement deep copy assignment.
    Model& operator=(const Model&) = delete;

    /// move constructor
    /// \todo Implement move construction.
    Model(Model&&) = delete;

    /// move assignment operator
    /// \todo Implement move assignment.
    Model& operator=(Model&&) = delete;

    /// Calculate model for each data point in the data partition
    /// \param D DataPartition to calculate over
    /// \todo This need not be a member function!
    void calculate(DataPartition& D) const;

    /// Check consistency of object
    virtual bool consistent() const;

    /// \return whether model has been locked
    bool locked() const
    { return Locked_; }

    /// prepare model and mark as locked:
    /// removes expired DataAccessor's, prune's remaining, and assigns them indices;
    /// fixes amplitudes that needn't be free
    void lock();

    /// \name Getters
    /// @{

    /// \return coordinate system (const)
    const CoordinateSystem<double, 3>& coordinateSystem() const
    { return CoordinateSystem_; }

    /// \return FourMomenta accessor
    std::shared_ptr<FourMomenta> fourMomenta()
    { return FourMomenta_; }

    /// \return FourMomenta accessor (const)
    const std::shared_ptr<FourMomenta> fourMomenta() const
    { return FourMomenta_; }

    /// \return MeasuredBreakupMomenta accessor
    std::shared_ptr<MeasuredBreakupMomenta> measuredBreakupMomenta()
    { return MeasuredBreakupMomenta_; }

    /// \return MeasuredBreakupMomenta accessor (const)
    const std::shared_ptr<MeasuredBreakupMomenta> measuredBreakupMomenta() const
    { return MeasuredBreakupMomenta_; }

    /// \return HelicityAngles accessor
    /// is nullptr if not using helicity formalism
    std::shared_ptr<HelicityAngles> helicityAngles()
    { return HelicityAngles_; }

    /// \return HelicityAngles accessor (const)
    /// is nullptr if not using helicity formalism
    const std::shared_ptr<HelicityAngles> helicityAngles() const
    { return HelicityAngles_; }

    /// \return ParticleCombinationCache
    ParticleCombinationCache& particleCombinationCache()
    { return ParticleCombinationCache_; }

    /// \return ParticleCombinationCache (const)
    const ParticleCombinationCache& particleCombinationCache() const
    { return ParticleCombinationCache_; }

    /// \return SpinAmplitudeCache
    SpinAmplitudeCache* spinAmplitudeCache()
    { return SpinAmplitudeCache_.get(); }

    /// \return SpinAmplitudeCache (const)
    const SpinAmplitudeCache* spinAmplitudeCache() const
    { return SpinAmplitudeCache_.get(); }

    /// \return vector of shared pointers to final state particles
    const FinalStateParticleVector& finalStateParticles() const
    { return FinalStateParticles_; }

    /// \return InitialStateParticles_
    const InitialStateParticleMap& initialStateParticles() const
    { return InitialStateParticles_; }

    /// \return set of DataAccessors
    const DataAccessorSet& dataAccessors() const
    { return DataAccessors_; }

    /// @}

    /// \name Setters
    /// @{


    /// Set final-state particle content. The order in which particles
    /// are given dictates the order in which four-momenta must be
    /// given in data points. The FinalStateParticle's have their
    /// Model pointer set to this
    /// \param FSP vector of shared pointers to final-state particles
    void setFinalState(const FinalStateParticleVector& FSP);

    /// Set final-state particle content. The order in which particles
    /// are given dictates the order in which four-momenta must be
    /// given in data points. The FinalStateParticle's have their
    /// Model pointer set to this
    /// \param FSPs shared pointers to final-state particles
    template <typename ... Types>
    void setFinalState(Types ... FSPs)
    { FinalStateParticleVector V{FSPs...}; setFinalState(V); }

    /// add an initial state particle
    /// The first particle added that decays to the full final state will become THE initial state particle,
    /// which can be retrieved by calling #initialStateParticle,
    /// and its amplitude will be fixed
    const InitialStateParticleMap::value_type& addInitialStateParticle(std::shared_ptr<DecayingParticle> bg);

    /// set four momenta of data point
    /// \param P vector of FourVectors of final-state momenta
    /// \param sm StatusManager to update
    void setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, StatusManager& sm) const;

    /// set coordinate system
    void setCoordinateSystem(const CoordinateSystem<double, 3>& cs);

    /// @}

    /// \name Monte Carlo Generation
    /// @{

    /// Build vector of mass axes for coordinates in phase space.
    /// Currently only supports two-particle masses.
    /// if argument is left empty, a default set of axes is constructed
    /// \return MassAxes for requested particle combinations
    /// \param pcs vector of vectors of particle indices
    /// \todo Find way to make const.
    const MassAxes massAxes(std::vector<std::vector<unsigned> > pcs = {});

    /// @}

    /// tell model to calculate and store helicity angles
    void requireHelicityAngles();

    /// tell model to calculate and store measured breakup momenta
    void requireMeasuredBreakupMomenta();

    /// create an empty data set
    /// \param n Number of empty data points to place inside data set
    DataSet createDataSet(size_t n = 0);

    /// Set VariableStatus'es of all Parameter's to unchanged, or leave as fixed
    void setParameterFlagsToUnchanged();

    /// Print the list of DataAccessor's
    void printDataAccessors(bool printParticleCombinations = true) const;

    /// Print all VariableStatus'es and CalculationStatus'es
    void printFlags(const StatusManager& sm) const;

    /// grant friend status to DataAccessor to register itself with this
    friend class DataAccessor;

    /// grant friend status to DataAccessor to register itself with this
    friend class RecalculableDataAccessor;

    /// grant friend status to DataAccessor to register itself with this
    friend class StaticDataAccessor;

    /// grant friend status to DecayingParticle to call addParticleCombination
    friend class DecayingParticle;

protected:

    /// add ParticleCombination to to FourMomenta_, HelicityAngles_, and
    /// MeasuredBreakupMomenta_ (along with it's daughters through
    /// recursive calling) if it is NOT for a FSP.
    virtual void addParticleCombination(std::shared_ptr<ParticleCombination> pc);

    /* /// remove a DataAccessor from this Model */
    /* virtual void removeDataAccessor(DataAccessorSet::value_type da); */

private:

    /// stores whether model structure can be modified
    /// (whether DataAccessors can be added or not)
    bool Locked_;

    /// Lab coordinate system to use in calculating helicity angles
    CoordinateSystem<double, 3> CoordinateSystem_;

    /// ParticleCombination cache
    ParticleCombinationCache ParticleCombinationCache_;

    /// SpinAmplitude cache
    std::unique_ptr<SpinAmplitudeCache> SpinAmplitudeCache_;

    /// Set of all DataAccessor's registered to this model
    DataAccessorSet DataAccessors_;

    /// Vector of raw pointers to all StaticDataAccessor's registered
    /// to this model, in order in which they need be calculated
    StaticDataAccessorVector StaticDataAccessors_;

    /// set of pointers to RecalculableDataAccessors
    RecalculableDataAccessorSet RecalculableDataAccessors_;

    /// pointers to initial particles
    /// they will be summed in incoherently
    InitialStateParticleMap InitialStateParticles_;

    /// vector of final state particles
    FinalStateParticleVector FinalStateParticles_;

    /// four momenta manager
    std::shared_ptr<FourMomenta> FourMomenta_;

    /// Breakup momenta manager
    std::shared_ptr<MeasuredBreakupMomenta> MeasuredBreakupMomenta_;

    /// helicity angles manager
    std::shared_ptr<HelicityAngles> HelicityAngles_;

};

/// \return string of AdmixtureMap
std::string to_string(const AdmixtureMap& mix);

/// \return whether a particle decays to its model's full final state
const bool decays_to_full_final_state(const Particle& p);

/// \return vector of shared_ptr to DecayingParticles inside Model
/// that decay to its full final state, sorted such that the first
/// entries have fixed prefactors
std::vector<std::shared_ptr<DecayingParticle> > full_final_state_isp(const Model& M);

/// \return intensity for all spin projections of an ISP
const double intensity(const InitialStateParticleMap::value_type& isp_mix, const DataPoint& d);

/// \return intensity for a data point evaluated over isp_map
const double intensity(const InitialStateParticleMap& isp_map, const DataPoint& d);

/// \return The sum of the logs of squared amplitudes evaluated over the data partition
/// \param M Model to evaluate
/// \param D DataPartition to evalue over
/// \param ped Pedestal to substract from each term in the sum
const double sum_of_log_intensity(const Model& M, DataPartition& D, double ped = 0);

/// \return The sum of the logs of squared amplitudes evaluated over the data partitions
/// \param DP DataPartitionVector of partitions to use
/// \param ped Pedestal to substract from each term in the sum
const double sum_of_log_intensity(const Model& M, DataPartitionVector& DP, double ped = 0);

/// \return all free amplitudes in a model
FreeAmplitudeSet free_amplitudes(const Model& M);

/// \return free amplitude in a model from decay trees evaluating to true
template <typename Last, typename ... UnaryPredicates>
FreeAmplitudeSet free_amplitudes(const Model& M, Last p, UnaryPredicates ... P)
{ return filter(free_amplitudes(M), p, P...); }

/// \return lone free amplitude in a model passing predicates
/// throws if unique amplitude is not found
template <typename Last, typename ... UnaryPredicates>
FreeAmplitudeSet::value_type free_amplitude(const Model& M, Last p, UnaryPredicates ... P)
{ return lone_elt(filter(free_amplitudes(M), p, P...)); }

/// \return set of all particles in model
ParticleSet particles(const Model& M);

/// \return Set of particles in model for which all predicates evaluate true
template <typename Last, typename ... UnaryPredicates>
ParticleSet particles(const Model& M, Last p, UnaryPredicates ... P)
{ return filter(particles(M), p, P...); }

/// \return lone particle in model for which all predicates evaluate true
/// throws if unique particle is not found
template <typename Last, typename ... UnaryPredicates>
ParticleSet::value_type particle(const Model& M, Last p, UnaryPredicates ... P)
{ return lone_elt(filter(particles(M), p, P...)); }

}

#endif
