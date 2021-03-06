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

#ifndef yap_InitialStateParticle_h
#define yap_InitialStateParticle_h

#include "CoordinateSystem.h"
#include "DataSet.h"
#include "DecayingParticle.h"
#include "FourMomenta.h"
#include "FourVector.h"
#include "HelicityAngles.h"
#include "MeasuredBreakupMomenta.h"

#include <complex>
#include <memory>
#include <vector>

namespace yap {

class DataPartitionBase;
class FinalStateParticle;

/// \class InitialStateParticle
/// \brief Class implementing an initial state particle.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class InitialStateParticle : public DecayingParticle
{
public:

    /// \name Constructor, destructor, clone
    /// @{

    /// Constructor
    InitialStateParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize);

    /// Destructor
    ~InitialStateParticle();

    /// @}

    /// \name Amplitude-related
    /// @{

    /// inherit std::complex<double> DecayingParticle::amplitude(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, unsigned dataPartitionIndex) const
    using DecayingParticle::amplitude;

    /// \return amplitude with a sum over all particle combinations
    std::complex<double> amplitude(DataPoint& d, unsigned dataPartitionIndex) const;

    /// \return ln(|amplitude|^2), with sum over all particle combinations in amp. calculation
    /// calls resetCalculationStatuses before calculation
    double logOfSquaredAmplitude(DataPoint& d, unsigned dataPartitionIndex)
    {
        resetCalculationStatuses(dataPartitionIndex);
        std::complex<double> a = amplitude(d, dataPartitionIndex);

        if (isnan(a.real() or isnan(a.imag())))
            return -std::numeric_limits<double>::infinity();

        return log(norm(a));
    }

    /// \return The sum of the logs of squared amplitudes evaluated over the data partition
    /// \param D Pointer to a #DataPartitionBase object
    double partialSumOfLogsOfSquaredAmplitudes(DataPartitionBase* D);

    /// Calculate the sum of the logs of the squared amplitudes evaluated over all partitions
    double sumOfLogsOfSquaredAmplitudes();

    /// Stores the sum of the lqogs of the squared amplitudes evaluated over the data partition
    /// \param D Pointer to a #DataPartitionBase
    /// \param s Double to store sum of log of squared amplitudes in
    void storeSumOfLogsOfSquaredAmplitudes(DataPartitionBase* D, double& s)
    { s = partialSumOfLogsOfSquaredAmplitudes(D); }

    /// @}

    /// call before looping over all DataPartitions
    void updateGlobalCalculationStatuses();

    /// calculate FourMomenta_, MeasuredBreakupMomenta_ and HelicityAngles_
    void calculate(DataPoint& d);

    /// set all parameter flags to kUnchanged (or leave at kFixed)
    /// call after looping over ALL DataPartitions
    void setParameterFlagsToUnchanged();

    /// Check consistency of object
    virtual bool consistent() const override;

    /// you MUST call this function after you have added all decay channels and before adding DataPoints
    bool prepare();

    /// \name Getters
    /// @{

    /// \return coordinate system
    CoordinateSystem<double, 3>& coordinateSystem()
    { return CoordinateSystem_; }

    /// \return coordinate system (const)
    const CoordinateSystem<double, 3>& coordinateSystem() const
    { return CoordinateSystem_; }

    /// \return FourMomenta accessor
    FourMomenta& fourMomenta()
    { return FourMomenta_; }

    /// \return FourMomenta accessor (const)
    const FourMomenta& fourMomenta() const
    { return FourMomenta_; }

    /// \return MeasuredBreakupMomenta accessor
    MeasuredBreakupMomenta& measuredBreakupMomenta()
    { return MeasuredBreakupMomenta_; }

    /// \return MeasuredBreakupMomenta accessor (const)
    const MeasuredBreakupMomenta& measuredBreakupMomenta() const
    { return MeasuredBreakupMomenta_; }

    /// \return HelicityAngles accessor
    HelicityAngles& helicityAngles()
    { return HelicityAngles_; }

    /// \return HelicityAngles accessor (const)
    const HelicityAngles& helicityAngles() const
    { return HelicityAngles_; }

    /// \return vector of shared pointers to final state particles
    const std::vector<std::shared_ptr<FinalStateParticle> >& finalStateParticles() const
    { return FinalStateParticles_; }

    /// \return (min, max) array[2] of mass range for particle combination
    /// \param pc shared pointer to ParticleCombination to get mass range of
    std::array<double, 2> getMassRange(const std::shared_ptr<const ParticleCombination>& pc) const;

    /// \return if prepare() has been called for this InitialStateParticle
    bool prepared() const
    {return Prepared_; }

    /// \return free amplitudes of DecayChannels_
    std::vector<std::shared_ptr<ComplexParameter> > freeAmplitudes() const;

    /// @}

    /// \name Setters
    /// @{

    /// Set final-state particle content. The order in which particles
    /// are given dictates the order in which four-momenta must be
    /// given in data points
    /// \param FSP list of shared pointers to final-state particles
    /// \return Success of action
    bool setFinalStateParticles(std::initializer_list<std::shared_ptr<FinalStateParticle> > FSP);

    /// @}

    /// \name Data set and partitions
    /// @{

    /// \return DataSet
    DataSet& dataSet()
    { return DataSet_; }

    /// \return the data partitions
    std::vector<DataPartitionBase*> dataPartitions();

    /// set data partitions
    /// ownership over DataPartitionBase objects will be taken
    void setDataPartitions(std::vector<std::unique_ptr<DataPartitionBase> > partitions);

    /// @}

    /// \name DataPoints
    /// @{

    /// Add data point via four-momenta
    /// This method is faster since it avoids unneccessary copying of objects
    /// and resizing of the DataPoint's storage
    bool addDataPoint(const std::vector<FourVector<double> >& fourMomenta);

    /// Add data point via move
    /// \param d DataPoint to move into DataSet
    /// \return Success of action
    bool addDataPoint(DataPoint&& d);

    /// Add data point via copy
    /// \param d DataPoint to copy into DataSet
    /// \return Success of action
    bool addDataPoint(const DataPoint& d);

    /// @}

    /// \name Monte Carlo Generation
    /// @{

    /// Initialize DataSet for MC Generation
    /// \param n Number of simultaneous streams for MC generation
    bool initializeForMonteCarloGeneration(unsigned n);

    /// @}

    /// Print the list of DataAccessor's
    void printDataAccessors(bool printParticleCombinations = true);

private:

    /// \name InitialStateParticle friends
    /// @{

    friend class DataAccessor;
    friend class AmplitudeComponent;
    friend void DecayingParticle::optimizeSpinAmplitudeSharing();

    /// @}

    /// check if d is in DataPartitions_
    bool hasDataPartition(DataPartitionBase* d);

    /// set number of data partitions of all #CachedDataValue's
    void setNumberOfDataPartitions(unsigned n);

    /// Set parents of symmetrization indices (recursively)
    virtual void setSymmetrizationIndexParents() override;

    /// reset all CalculationStatus'es for the dataPartitionIndex to the GlobalCalculationStatus_
    /// call before calculating the amplitude for a new dataPoint
    void resetCalculationStatuses(unsigned dataPartitionIndex);

    /// set all parameter flags to kUnchanged (or leave at kFixed)
    /// call after looping over a DataPartition
    void setCachedDataValueFlagsToUnchanged(unsigned dataPartitionIndex);

    /// add DataAccessor to set
    void addDataAccessor(DataAccessor* d);

    /// remove DataAccessor from set
    void removeDataAccessor(DataAccessor* d);

    /// add AmplitudeComponent to set
    void addAmplitudeComponent(AmplitudeComponent* d)
    { AmplitudeComponents_.insert(d); }

    /// remove AmplitudeComponent from set
    void removeAmplitudeComponent(AmplitudeComponent* d);

    void setDataAcessorIndices();

    /// Stores whether prepare() has been called
    bool Prepared_;

    /// Lab coordinate system to use in calculating helicity angles
    CoordinateSystem<double, 3> CoordinateSystem_;

    /// List of all DataAccessor objects in the InitialsStateParticle and below
    std::set<DataAccessor*> DataAccessors_;

    /// List of all AmplitudeComponent objects in the InitialsStateParticle and below
    std::set<AmplitudeComponent*> AmplitudeComponents_;

    /// List of all DecayChannel objects in the InitialsStateParticle and below
    std::vector<DecayChannel*> DecayChannels_;

    /// vector of final state particles
    std::vector<std::shared_ptr<FinalStateParticle> > FinalStateParticles_;

    /// four momenta manager
    FourMomenta FourMomenta_;

    /// Breakup momenta manager
    MeasuredBreakupMomenta MeasuredBreakupMomenta_;

    /// helicity angles manager
    HelicityAngles HelicityAngles_;

    /// Data set
    DataSet DataSet_;

    /// Data partitions
    std::vector<std::unique_ptr<DataPartitionBase> > DataPartitions_;

};

}

#endif
