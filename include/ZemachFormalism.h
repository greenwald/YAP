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

#ifndef yap_ZemachFormalism_h
#define yap_ZemachFormalism_h

#include "fwd/DataPoint.h"
#include "fwd/ParticleCombination.h"
#include "fwd/Spin.h"

#include "SpinAmplitude.h"
#include "SpinAmplitudeCache.h"
#include "UnitSpinAmplitude.h"

#include <complex>
#include <memory>

namespace yap {

/// \class ZemachSpinAmplitude
/// \brief Class implementing Zemach tensors
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup SpinAmplitude
class ZemachSpinAmplitude : public SpinAmplitude
{
public:

    /// Calculate spin amplitude for given ParticleCombination and spin projections
    /// \param two_M 2 * spin projection of parent
    /// \param two_m2 SpinProjectionVector of daughters
    /// \param d DataPoint to retrieve data from for calculation
    /// \param pc ParticleCombination to calculate for
    virtual const std::complex<double> calc(int two_M, const SpinProjectionVector& two_m,
                                            const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const override;

    /// \return "Zemach formalism"
    virtual std::string formalism() const override
    { return "Zemach formalism"; }

    /// grant SpinAmplitudeCache friend status to call constructor
    friend class ZemachFormalism;

protected:
    /// Constructor
    /// \param two_J twice the spin of initial state
    /// \param two_j SpinVector of daughters
    /// \param l orbital angular momentum
    /// \param two_s twice the total spin angular momentum
    ZemachSpinAmplitude(unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s);

private:
    /// check equality
    virtual bool equals(const SpinAmplitude& other) const override
    { return dynamic_cast<const ZemachSpinAmplitude*>(&other) and SpinAmplitude::equals(other); }

};

/// \class ZemachFormalism
/// \brief Caches ZemachSpinAmplitude's
/// \author Daniel Greenwald
class ZemachFormalism : public SpinAmplitudeCache
{
public:

    /// Constructor
    ZemachFormalism() : SpinAmplitudeCache() {}

private:

    /// override in inherting classes
    /// \return shared_ptr to SpinAmplitude object
    /// \param two_J twice the spin of initial state
    /// \param two_j SpinVector of daughters
    /// \param L orbital angular momentum
    /// \param two_S 2 * the total spin angular momentum
    virtual std::shared_ptr<SpinAmplitude> create(unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s) const override
    { return (two_s == 0) ? unit(two_J, two_j, l, two_s) : std::shared_ptr<SpinAmplitude>(new ZemachSpinAmplitude(two_J, two_j, l, two_s)); }

};


}

#endif
