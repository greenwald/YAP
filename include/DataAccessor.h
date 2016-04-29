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

#include "CachedDataValue.h"
#include "ReportsParticleCombinations.h"
#include "ParticleCombination.h"

#include <memory>
#include <set>

namespace yap {

class Model;

/// \name DataAccessor
/// \brief Abstract base class for all objects accessing DataPoint's
/// \author Johannes Rauch, Daniel Greenwald

class DataAccessor :
    public virtual ReportsParticleCombinations
{
public:

    /// Constructor
    /// \param equiv ParticleCombination equivalence struct for determining index assignments
    DataAccessor(ParticleCombination::Equiv* equiv = &ParticleCombination::equivBySharedPointer);

    /// \name Access to indices
    /// @{

    /// \return index inside DataPoint structure that this DataAccessor accesses
    int index() const
    { return Index_; }

    /// \return if the given ParticleCombination is in SymmetrizationIndices_
    bool hasParticleCombination(const std::shared_ptr<ParticleCombination>& c) const
    { return SymmetrizationIndices_.find(c) != SymmetrizationIndices_.end(); }

    /// \return if the given ParticleCombination is in SymmetrizationIndices_
    /// \param c ParticleCombination to look for equivalent of
    /// \param equiv ParticleCombination::Equiv object for checking equivalence
    bool hasParticleCombination(const std::shared_ptr<ParticleCombination>& c,
                                const ParticleCombination::Equiv& equiv) const;

    /// \return index inside row of DataPoint for the requested ParticleCombination
    unsigned symmetrizationIndex(const std::shared_ptr<ParticleCombination>& c) const
    { return SymmetrizationIndices_.at(c); }

    /// \return SymmetrizationIndices_
    const ParticleCombinationMap<unsigned>& symmetrizationIndices() const
    { return SymmetrizationIndices_; }

    /// \return maximum index of SymmetrizationIndices_
    /// -1 means empty
    int maxSymmetrizationIndex() const;

    /// \return vector of ParticleCombination's
    const ParticleCombinationVector& particleCombinations() const override
    { return ParticleCombinations_; }

    /// print ParticleCombination map
    void printParticleCombinations() const;

    /// @}

    /// \return CachedDataValueSet
    const CachedDataValueSet& cachedDataValues() const
    { return CachedDataValues_; }

    /// \return size of storage in data point (number of real values)
    unsigned size() const
    { return Size_; }

    /// Check consistency of object
    bool consistent() const;

    /// get raw pointer to Model (const)
    virtual const Model* model() const = 0;

    /// \return string denoting DataAccessor type
    /// \todo REMOVE
    virtual std::string data_accessor_type() const = 0;

    /// grant friend status to Model to access CachedDataValues_
    friend class Model;

    /// grant friend status to CachedDataValue to call addCachedDataValue
    friend class CachedDataValue;

protected:

    /// register with Model
    void virtual addToModel();

    /// add CachedDataValue
    void addCachedDataValue(std::shared_ptr<CachedDataValue> c);

    /// Increase storage
    /// \param n number of elements to increase by
    void increaseSize(unsigned n)
    { Size_ += n; }

    /// add ParticleCombination to SymmetrizationIndices_
    virtual unsigned addParticleCombination(std::shared_ptr<ParticleCombination> pc) override;

    /// prune SymmetrizationIndices_ to only contain ParticleCombination's tracing back up the ISP
    virtual void pruneSymmetrizationIndices();

    /// set storage index used in DataPoint. Must be unique.
    void setIndex(size_t i)
    { Index_ = i; }

private:

    /// Object to check equality of symmetrizations for determining storage indices
    ParticleCombination::Equiv* Equiv_;

    /// Map of indices for each used symmetrization stored with key = shared_ptr<ParticleCombination>
    ParticleCombinationMap<unsigned> SymmetrizationIndices_;

    /// Vector of particle combinations. This is a cache kept for performance reasons,
    /// it must always be in sync with SymmetrizationIndices_
    ParticleCombinationVector ParticleCombinations_;

    /// Set of CachedDataValues that have this DataAccessor as an owner
    CachedDataValueSet CachedDataValues_;

    /// number of real values stored per symm. index
    unsigned Size_;

    /// storage index used in DataPoint. Must be unique.
    int Index_;

};

/// \typedef DataAccessorSet
using DataAccessorSet = std::set<DataAccessor*>;
// using DataAccessorSet = std::set<std::shared_ptr<DataAccessor>, std::owner_less<std::shared_ptr<DataAccessor> > >;

/// remove expired elements of set
void removeExpired(DataAccessorSet& S);

}

#endif
