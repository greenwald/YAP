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

#ifndef yap_BlattWeisskopf_h
#define yap_BlattWeisskopf_h

#include "fwd/DataPointBase.h"

#include "AmplitudeComponent.h"
#include "DataAccessor.h"
#include "RequiresMeasuredBreakupMomenta.h"

#include <memory>
#include <string>

namespace yap {

class DataPartition;
class DecayingParticle;
class Model;
class ParticleCombination;
class RealCachedDataValue;
class StatusManager;

/// \class BlattWeisskopf
/// \brief Class implementing BlattWeisskopf barrier factors
/// \author Johannes Rauch, Daniel Greenwald

class BlattWeisskopf : public AmplitudeComponent, public DataAccessor, public RequiresMeasuredBreakupMomenta
{
public:

    /// Constructor
    /// \param L angular momentum of Blatt-Weisskopf barrier factor
    /// \param dp raw pointer to owning DecayingParticle
    BlattWeisskopf(unsigned L, DecayingParticle* dp);

    /// \return angular momentum
    unsigned L() const
    { return L_; }

    /// Calculate amplitude
    /// \param d pointer to DataPointBase to calculate with
    /// \param pc (shared_ptr to) ParticleCombination to calculate for
    /// \param sm StatusManager to update
    virtual double amplitude(DataPointBase& d, const std::shared_ptr<ParticleCombination>& pc, StatusManager& sm) const;

    /// functor
    /// \return Blatt-Weisskopf barrier factor for data point and particle combination
    /// \param d pointer to DataPointBase
    /// \param pc shared_ptr to ParticleCombination
    double operator()(DataPointBase& d, const std::shared_ptr<ParticleCombination>& pc) const;

    /// Calculate barrier factor for and store into each data point in a partition
    /// \param D DataPartition to calculate on
    /// \param pc (shared_ptr to) ParticleCombination to calculate for
    virtual void calculate(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc) const;

    /// check consistency of object
    virtual bool consistent() const override
    { return DataAccessor::consistent(); }

    virtual CachedDataValueSet cachedDataValuesItDependsOn() override;

    /// \return raw pointer to Model through owning DecayingParticle
    const Model* model() const override;

    virtual std::string data_accessor_type() const override
    {return "BlattWeisskopf"; }

    /// grant friend status to DecayingParticle to call addParticleCombination
    friend class DecayingParticle;

private:

    /// raw pointer to owning DecayingParticle
    DecayingParticle* DecayingParticle_;

    /// angular momentum
    unsigned L_;

    /// Blatt-Weisskopf barrier factor
    std::shared_ptr<RealCachedDataValue> BarrierFactor_;

};

}

#endif
