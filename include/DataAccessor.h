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

#include "fwd/DataAccessor.h"

#include "fwd/CachedValue.h"
#include "fwd/DataPoint.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"

#include <memory>

namespace yap {

/// \class DataAccessor
/// \brief Abstract base class for all objects accessing DataPoint's
/// \author Johannes Rauch, Daniel Greenwald
class DataAccessor
{
public:

    /// Constructor
    /// \param equal ParticleCombination equality struct for determining index assignments
    DataAccessor(const ParticleCombinationEqualTo& equal);

    /// \return Equality function
    const ParticleCombinationEqualTo& equal() const
    { return Equal_; }

    /// \return index inside DataPoint structure that this DataAccessor accesses
    int index() const
    { return Index_; }

    /// \return index inside row of DataPoint for the requested ParticleCombination
    unsigned symmetrizationIndex(const std::shared_ptr<const ParticleCombination>& c) const
    { return SymmetrizationIndices_.at(c); }

    /// \return SymmetrizationIndices_
    const ParticleCombinationMap<unsigned>& symmetrizationIndices() const
    { return SymmetrizationIndices_; }

    /// \return number of independent indices in SymmetrizationIndices_
    const unsigned nSymmetrizationIndices() const
    { return NIndices_; }

    /// \return CachedValueSet
    const CachedValueSet& CachedValues() const
    { return CachedValues_; }

    /// \return size of storage in data point (number of real values pr symm Index)
    const unsigned size() const
    { return Size_; }

    /// \return whether DataAccessor stores any data
    const bool requiresStorage() const
    { return size() > 0 and nSymmetrizationIndices() > 0; }

    /// Check consistency of object
    bool consistent() const;

    /// get raw pointer to Model (const)
    virtual const Model* model() const = 0;

    /// grant friend status to Model to access CachedValues_
    friend Model;

    /// grant friend status to CachedValue to call addCachedValue
    friend CachedValue;

protected:

    /// register with Model
    void virtual registerWithModel();

    /// add CachedValue
    void addCachedValue(CachedValue& c);

    /// add ParticleCombination to SymmetrizationIndices_
    virtual void addParticleCombination(const ParticleCombination& pc);

    /// prune SymmetrizationIndices_ to only contain ParticleCombination's tracing back up the ISP
    virtual void pruneSymmetrizationIndices();

    /// set storage index used in DataPoint. Must be unique.
    void setIndex(size_t i)
    { Index_ = i; }

private:

    /// Increase storage
    /// \param n number of elements to increase by
    void increaseSize(unsigned n)
    { Size_ += n; }

    /// Object to check equality of symmetrizations for determining storage indices
    ParticleCombinationEqualTo Equal_;

    /// Map of indices for each used symmetrization stored with key = shared_ptr<ParticleCombination>
    ParticleCombinationMap<unsigned> SymmetrizationIndices_;

    /// Number of independent indices stored in SymmetrizationIndices_
    unsigned NIndices_;

    /// Set of CachedValues that have this DataAccessor as an owner
    CachedValueSet CachedValues_;

    /// number of real values stored per symm. index
    unsigned Size_;

    /// storage index used in DataPoint. Must be unique.
    int Index_;

};

/// remove expired elements of set
void remove_expired(DataAccessorSet& S);

/// \todo REMOVE
std::string to_string(const DataAccessorSet& S, bool print_particle_combinations = true);

}

#endif
