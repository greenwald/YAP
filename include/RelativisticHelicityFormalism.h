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

#ifndef yap_RelativisticHelicityFormalism_h
#define yap_RelativisticHelicityFormalism_h

#include "fwd/DataPoint.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/Spin.h"

#include "RequiresHelicityAngles.h"
#include "SpinAmplitude.h"
#include "SpinAmplitudeCache.h"

#include <complex>
#include <map>
#include <memory>

namespace yap {

/// \class RelativisticHelicitySpinAmplitude
/// \brief Class implementing a canonical spin amplitude, i.e. with defined relative angular momentum.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup SpinAmplitude
class RelativisticHelicitySpinAmplitude : public SpinAmplitude, public RequiresHelicityAngles
{
protected:

    /// Constructor
    /// \param m Owning model
    /// \param two_J twice the spin of initial state
    /// \param two_j SpinVector of daughters
    /// \param l orbital angular momentum
    /// \param two_s twice the total spin angular momentum
    RelativisticHelicitySpinAmplitude(Model& m, unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s);

public:

    /// \return precalculated complex amplitude
    /// \param d DataPoint to retrieve value from
    /// \param pc ParticleCombination to retrieve value for
    /// \param two_M 2 * spin projection of parent
    /// \param two_m SpinProjectionVector of daughters
    virtual const std::complex<double> amplitude(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc,
                                         int two_M, const SpinProjectionVector& two_m) const override
    { return (initialTwoJ() == 0) ? Coefficients_.at(two_m) : SpinAmplitude::amplitude(d, pc, two_M, two_m); }

    /// Overrides SpinAmplitude::calculate to do nothing if initialTwoJ() == 0
    /// \param d DataPoint to calculate into
    /// \param sm StatusManager to update
    virtual void calculate(DataPoint& d, StatusManager& sm) const override
    { if (initialTwoJ() != 0) SpinAmplitude::calculate(d, sm); }

    /// Calculate spin amplitude for given ParticleCombination and spin projections
    /// \param two_M 2 * spin projection of parent
    /// \param two_m SpinProjectionVector of daughters
    /// \param d DataPoint to retrieve data from for calculation
    /// \param pc ParticleCombination to calculate for
    virtual const std::complex<double> calc(int two_M, const SpinProjectionVector& two_m,
                                            const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const override;

    /// check equality
    virtual bool equals(const SpinAmplitude& other) const override
    { return dynamic_cast<const RelativisticHelicitySpinAmplitude*>(&other) and SpinAmplitude::equals(other); }

    /// \return "helicity formalism"
    virtual std::string formalism() const override
    { return "helicity formalism"; }

    /// grant RelativisticHelicityFormalism friend status to call constructor
    friend class RelativisticHelicityFormalism;

private:

    /// map from (m1,m2) -> L-S * S-S coupling coefficients
    /// value is sqrt((2L+1)/4pi) * (L 0 S m1-m2 | J m1-m2) * (j1 m1 j2 m2 | S m1-m2);
    std::map<SpinProjectionVector, double> Coefficients_;

};

/// \class RelativisticHelicityFormalism
/// \brief Caches RelativisticHelicitySpinAmplitude's
/// \author Daniel Greenwald
class RelativisticHelicityFormalism : public SpinAmplitudeCache
{
public:

    /// Constructor
    RelativisticHelicityFormalism() : SpinAmplitudeCache() {}

private:

    /// override in inherting classes
    /// \return shared_ptr to SpinAmplitude object
    /// \param m Owning model
    /// \param two_J twice the spin of initial state
    /// \param two_j SpinVector for daughters
    /// \param L orbital angular momentum
    /// \param two_S 2 * the total spin angular momentum
    virtual std::shared_ptr<SpinAmplitude> create(Model& m, unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s) const override
    { return std::shared_ptr<SpinAmplitude>(new RelativisticHelicitySpinAmplitude(m, two_J, two_j, l, two_s)); }

};

}

#endif
