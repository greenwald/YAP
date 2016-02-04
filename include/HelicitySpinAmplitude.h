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

#ifndef yap_HelicitySpinAmplitude_h
#define yap_HelicitySpinAmplitude_h

#include "DataPoint.h"
#include "SpinAmplitude.h"

#include <complex>
#include <map>
#include <memory>

namespace yap {

class ParticleCombination;

template <class T>
class SpinAmplitudeCache;

/// \class HelicitySpinAmplitude
/// \brief Class implementing a canonical spin amplitude, i.e. with defined relative angular momentum.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup SpinAmplitude
class HelicitySpinAmplitude : public SpinAmplitude
{
public:

    /// Calculate spin amplitude for given ParticleCombination and spin projections
    /// \param two_M 2 * spin projection of parent
    /// \param two_m1 2 * spin projection of first daughter
    /// \param two_m2 2 * spin projection of second daughter
    /// \param d DataPoint to retrieve data from for calculation
    /// \param pc ParticleCombination to calculate for
    virtual std::complex<double> calc(int two_M, int two_m1, int two_m2,
                                      const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const override;

    /// \return "helicity formalism"
    virtual std::string formalism() const override
    { return "helicity formalism"; }

    /// grant SpinAmplitudeCache friend status to call constructor
    friend class SpinAmplitudeCache<HelicitySpinAmplitude>;

protected:

    /// Constructor
    /// \param two_J  twice the spin of Initial-state
    /// \param two_j1 twice the spin of first daughter
    /// \param two_j2 twice the spin of second daughter
    /// \param l orbital angular momentum
    /// \param two_s twice the total spin angular momentum
    /// \param isp raw pointer to owning InitialStateParticle
    HelicitySpinAmplitude(unsigned two_J, unsigned two_j1, unsigned two_j2, unsigned l, unsigned two_s,
                          InitialStateParticle* isp);

private:
    /// check equality
    virtual bool equals(const SpinAmplitude& other) const override
    { return dynamic_cast<const HelicitySpinAmplitude*>(&other) and SpinAmplitude::equals(other); }

    /// L-S * S-S coupling coefficients
    /// first map key is m1;
    /// second map key is m2;
    /// value is sqrt((2L+1)/4pi) * (L 0 S m1-m2 | J m1-m2) * (j1 m1 j2 m2 | S m1-m2);
    std::map<int, std::map<int, double> > Coefficients_;

};

}

#endif
