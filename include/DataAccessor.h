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

#ifndef yap_DataAccessor_h
#define yap_DataAccessor_h

#include "Amp.h"
#include "CalculationStatus.h"
#include "ParticleCombination.h"

#include <map>
#include <vector>

namespace yap {

class DataPoint;
class InitialStateParticle;

/// \name DataAccessor
/// \brief Base class for all objects accessing DataPoint's
/// \author Johannes Rauch, Daniel Greenwald

class DataAccessor
{
public:

    /// \todo Make a function to reserve space in the DataPoint?

    /// \name Constructors, destructor, & operators
    /// @{

    /// Constructor
    /// \param equiv ParticleCombination equivalence struct for determining index assignments
    DataAccessor(InitialStateParticle* isp, ParticleCombination::Equiv equiv = ParticleCombination::equivBySharedPointer);

    /// Copy constructor
    DataAccessor(const DataAccessor& other);

    // Defaulted move constructor
    // Defaulted destructor
    // Defaulted move assignment operator

    /// @}

    /// \name Access to indices
    /// @{

    /// \return index inside DataPoint structure that this DataAccessor accesses
    unsigned index() const
    { return Index_; }

    /// \return index inside row of DataPoint for the requested symmetrization
    unsigned symmetrizationIndex(std::shared_ptr<ParticleCombination> c) const
    { return SymmetrizationIndices_.at(c); }

    /// \return list of all ParticleCombinations
    std::vector<std::shared_ptr<ParticleCombination> > particleCombinations() const;

    /// \return calculation statuses
    std::vector<CalculationStatus>& CalculationStatuses()
    { return CalculationStatuses_; }

    /// \return calculation statuses (const)
    const std::vector<CalculationStatus>& CalculationStatuses() const
    { return CalculationStatuses_; }

    /// @}

    /// Check consistency of object
    bool consistent() const;

    /// \name Symmetrization functions
    /// @{

    /// add symmetrizationIndex to SymmetrizationIndices_
    virtual void addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c);

    /// @}

    /// \name Data access
    /// @{

    /// Access a data point's data (by friendship)
    std::vector<double>& data(DataPoint& d, unsigned i) const;

    /// Access a data point's data (by friendship) (const)
    const std::vector<double>& data(const DataPoint& d, unsigned i) const;

    /// Get pointer to the initial state particle
    InitialStateParticle* initialStateParticle() const
    { return InitialStateParticle_; }

    /// @}

protected:

    /// pointer to the initial state particle for access to FourMomenta, HelicityAngles etc.
    InitialStateParticle* InitialStateParticle_;

    /// Object to check equality of symmetrizations for determining storage indices
    ParticleCombination::Equiv Equiv_;

    /// vector of calculation statuses for each index in the data (as needed by #SymmetrizationIndices_
    std::vector<CalculationStatus> CalculationStatuses_;

    /// Map of indices for each used symmetrization stored with key = shared_ptr<ParticleCombination>
    std::map<std::shared_ptr<ParticleCombination>, unsigned, std::owner_less<std::shared_ptr<ParticleCombination> > > SymmetrizationIndices_;

private:

    /// storage index used in DataPoint. Must be unique.
    unsigned Index_;

    /// static counter for setting indices
    static unsigned GlobalIndex;
};

}

#endif
