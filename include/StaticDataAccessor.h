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

#ifndef yap_StaticDataAccessor_h
#define yap_StaticDataAccessor_h

#include "fwd/DataPointBase.h"

#include "DataAccessor.h"
#include "Exceptions.h"
#include "ParticleCombination.h"

namespace yap {

class StatusManager;

class Model;

/// \name StaticDataAccessor
/// \brief Base class for all data accessors that will only write to DataPointBase once at initial data loading
/// \author Johannes Rauch, Daniel Greenwald
class StaticDataAccessor : public DataAccessor
{
public:

    /// Constructor
    /// \param equiv ParticleCombination equivalence struct for determining index assignments
    StaticDataAccessor(ParticleCombination::Equiv* equiv = &ParticleCombination::equivBySharedPointer)
        : DataAccessor(equiv), Model_(nullptr) {}

    /// Constructor
    /// \param model Raw pointer to owning Model
    /// \param equiv ParticleCombination equivalence struct for determining index assignments
    StaticDataAccessor(Model* m, ParticleCombination::Equiv* equiv = &ParticleCombination::equivBySharedPointer)
        : DataAccessor(equiv), Model_(nullptr)
    {
        setModel(m);
    }

    /// Set the Model and add this DataAccessor to the Model's DataAccessorSet
    virtual void setModel(Model* m)
    {
        Model_ = m;

        if (!model())
            throw exceptions::Exception("Model unset", "StaticDataAccessor::StaticDataAccessor");

        // register with Model
        addToModel();
    }

    /// calculate cachedDataValues, store to DataPointBase, and update StatusManager.
    /// Must be overriden in derived classes.
    virtual void calculate(DataPointBase& d, StatusManager& sm) const = 0;

    /// \return Raw pointer to owning Model
    const Model* model() const override
    { return Model_; }

private:

    Model* Model_;

};

/// \typedef StaticDataAccessorVector
using StaticDataAccessorVector = std::vector<StaticDataAccessor*>;

/// remove expired elements of set
inline void removeExpiredStatic(StaticDataAccessorVector& S)
{
    for (auto it = S.begin(); it != S.end(); )
        if (!*it) it = S.erase(it);
        else ++it;
}

}

#endif
