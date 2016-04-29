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

#ifndef yap_MeasuredBreakupMomenta_h
#define yap_MeasuredBreakupMomenta_h

#include "fwd/DataPointBase.h"

#include "StaticDataAccessor.h"

#include <memory>
#include <string>

namespace yap {

class Model;
class ParticleCombination;
class RealCachedDataValue;
class StatusManager;

/// \class MeasuredBreakupMomenta
/// \brief Calculates, stores and gives access to breakup momenta (using measured masses)
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup SpinAmplitude
class MeasuredBreakupMomenta : public StaticDataAccessor
{
public:

    /// Constructor
    /// \param m Raw pointer to owning Model
    MeasuredBreakupMomenta(Model* m);

    /// Calculate breakup momenta for all possible symmetrization indices
    /// \param d DataPointBase to caluclate into
    /// \param sm StatusManager to update
    virtual void calculate(DataPointBase& d, StatusManager& sm) const override;

    /// Access squared breakup momentum
    /// \param d DataPointBase to get data from
    /// \param pc ParticleCombination to return breakup momentum of
    double q2(const DataPointBase& d, const std::shared_ptr<ParticleCombination>& pc) const
    { return Q2_->value(d, symmetrizationIndex(pc)); }

    /// Access breakup momentum
    /// \param d DataPointBase to get data from
    /// \param pc ParticleCombination to return breakup momentum of
    double q(const DataPointBase& d, const std::shared_ptr<ParticleCombination>& pc) const
    { return sqrt(q2(d, pc)); }

    /// Calculate breakup momentum from parent and daughter masses
    static double calcQ2(double m2_R, double m_a, double m_b);

    /// \return Breakup Momentum
    std::shared_ptr<RealCachedDataValue> breakupMomenta()
    { return Q2_; }

    /// \return Breakup Momentum (const)
    std::shared_ptr<RealCachedDataValue> breakupMomenta() const
    { return Q2_; }

    virtual std::string data_accessor_type() const override
    {return "MeasuredBreakupMomenta"; }

    /// grant friend status to Model to call addParticleCombination
    friend class Model;

protected:

    /// override to throw on adding final-state PC
    unsigned addParticleCombination(std::shared_ptr<ParticleCombination> pc) override;

    /// squared breakup momentum [GeV^2]
    std::shared_ptr<RealCachedDataValue> Q2_;

};

}

#endif
