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

#include "fwd/DataPoint.h"

#include "CoordinateSystem.h"
#include "DataPartition.h"
#include "DataSet.h"
#include "DecayTree.h"
#include "FourVector.h"
#include "ParticleCombinationCache.h"
#include "StaticDataAccessor.h"

#include <complex>
#include <memory>
#include <vector>

namespace yap {

class DecayingParticle;
class FinalStateParticle;
class FourMomenta;
class HelicityAngles;
class MassAxes;
class MeasuredBreakupMomenta;
class SpinAmplitudeCache;
class StatusManager;

/// \class Model
/// \brief Class implementing a PWA model
/// \author Johannes Rauch, Daniel Greenwald

class Model
{
public:

    /// Constructor
    /// \param SAC unique_ptr to SpinAmplitudeCache
    Model(std::unique_ptr<SpinAmplitudeCache> SAC);

    /// \name Amplitude-related
    /// @{

    /// \return amplitude with a sum over all particle combinations of ISP
    /// \param d DataPoint to calculate with
    /// \param two_m 2 * the spin projection of ISP to calculate for
    /// \param sm StatusManager to update
    std::complex<double> amplitude(DataPoint& d, int two_m, StatusManager& sm) const;

    /// \return amplitude with a sum over all particle combinations and spin projections of ISP
    /// \param d DataPoint to calculate with
    /// \param sm StatusManager to update
    std::complex<double> amplitude(DataPoint& d, StatusManager& sm) const;

    /// \return The sum of the logs of squared amplitudes evaluated over the data partition
    /// \param D DataPartition to evalue over
    /// \param global StatusManager to reset partition's statuses to
    double partialSumOfLogsOfSquaredAmplitudes(DataPartition& D, const StatusManager& global) const;

    /// Calculate the sum of the logs of the squared amplitudes evaluated over all partitions
    /// \param DS DataSet to evaluate over
    /// \param DP DataPartitionVector of partitions to use
    double sumOfLogsOfSquaredAmplitudes(DataSet& DS, DataPartitionVector& DP) const;

    /// Calculate the sum of the logs of the squared amplitudes evaluated over the whole data set
    /// \param DS DataSet to evaluate over
    double sumOfLogsOfSquaredAmplitudes(DataSet& DS) const;

    /// @}

    /// Check consistency of object
    virtual bool consistent() const;

    bool locked() const
    { return Locked_; }

    /// removes expired DataAccessor's, prune's remaining, and assigns them indices
    void prepareDataAccessors();

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

    /// \return Initial-state particle
    std::shared_ptr<DecayingParticle> initialStateParticle()
    { return InitialStateParticle_; }

    /// \return Initial-state particle (const)
    std::shared_ptr<DecayingParticle> initialStateParticle() const
    { return InitialStateParticle_; }

    /// \return vector of shared pointers to final state particles
    const std::vector<std::shared_ptr<FinalStateParticle> >& finalStateParticles() const
    { return FinalStateParticles_; }

    /// \return set of DataAccessors
    const DataAccessorSet& dataAccessors() const
    { return DataAccessors_; }

    /// \return set of DataAccessors
    const StaticDataAccessorVector& staticDataAccessors() const
    { return StaticDataAccessors_; }

    /// \return (min, max) array[2] of mass range for particle combination
    /// \param pc shared pointer to ParticleCombination to get mass range of
    std::array<double, 2> massRange(const std::shared_ptr<ParticleCombination>& pc) const;

    /// \return free amplitudes of DecayChannels_
    ComplexParameterVector freeAmplitudes() const;

    /// @}

    /// \name Setters
    /// @{

    /// Set initial-state particle
    /// \param isp shared pointer to initial-state particle
    void setInitialStateParticle(std::shared_ptr<DecayingParticle> isp);

    /// Set final-state particle content. The order in which particles
    /// are given dictates the order in which four-momenta must be
    /// given in data points. The FinalStateParticle's have their
    /// Model pointer set to this
    /// \param FSP list of shared pointers to final-state particles
    void  setFinalState(std::initializer_list<std::shared_ptr<FinalStateParticle> > FSP);

    /// set coordinate system
    void setCoordinateSystem(const CoordinateSystem<double, 3>& cs);

    /// @}

    /// \name Monte Carlo Generation
    /// @{

    /// Build vector of mass axes for coordinates in phase space.
    /// Currently only supports two-particle masses; the PCs put into
    /// the returned MassAxes will have their daughters sorted (i.e. (10) will become (01)).
    /// \return MassAxes for requested particle combinations
    /// \param pcs vector of vectors of particle indices
    const MassAxes massAxes(std::vector<std::vector<unsigned> > pcs);

    /// Calculate four-momenta for final-state particles for phase-space coordinate
    /// \param axes phase-space axes
    /// \param squared_masses phase-space coordinate
    std::vector<FourVector<double> > calculateFourMomenta(const MassAxes& axes, const std::vector<double>& squared_masses) const;

    /// @}

    /// create an empty data set
    /// \param n Number of empty data points to place inside data set
    DataSet createDataSet(size_t n = 0);

    /// create a vector of all possible DecayTree's
    DecayTreeVector decayTrees();

    /// Set VariableStatus'es of all Parameter's to calculated, or leave as fixed
    void setParameterFlagsToUnchanged();

    /// Print the list of DataAccessor's
    void printDataAccessors(bool printParticleCombinations = true) const;

    /// Print all VariableStatus'es and CalculationStatus'es
    void printFlags(const StatusManager& sm) const;

    /// grant friend status to DataAccessor to register itself with this
    friend class DataAccessor;

    /// grant friend status to DecayingParticle to call addParticleCombination
    friend class DecayingParticle;

protected:

    /// add ParticleCombination to to FourMomenta_, HelicityAngles_, and
    /// MeasuredBreakupMomenta_ (along with it's daughters through
    /// recursive calling) if it is NOT for a FSP.
    virtual void addParticleCombination(std::shared_ptr<ParticleCombination> pc);

    /// register a DataAccessor with this Model
    virtual void addDataAccessor(DataAccessorSet::value_type da);

    /* /// remove a DataAccessor from this Model */
    /* virtual void removeDataAccessor(DataAccessorSet::value_type da); */

private:

    // check if the fourMomenta produce the given invariant masses
    bool checkInvariantMasses(const MassAxes& axes, const std::vector<double>& squared_masses, const std::vector<FourVector<double> >& fourMomenta) const;

    /// lock model
    void lock()
    { Locked_ = true; }

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

    /// Raw pointer to initial-state particle
    std::shared_ptr<DecayingParticle> InitialStateParticle_;

    /// vector of final state particles
    std::vector<std::shared_ptr<FinalStateParticle> > FinalStateParticles_;

    /// four momenta manager
    std::shared_ptr<FourMomenta> FourMomenta_;

    /// Breakup momenta manager
    std::shared_ptr<MeasuredBreakupMomenta> MeasuredBreakupMomenta_;

    /// helicity angles manager
    std::shared_ptr<HelicityAngles> HelicityAngles_;

};

}

#endif
