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

#ifndef yap_RelativisticBreitWigner_h
#define yap_RelativisticBreitWigner_h

#include "fwd/BlattWeisskopf.h"
#include "fwd/DataPartition.h"
#include "fwd/DecayChannel.h"
#include "fwd/ParticleCombination.h"

#include "BreitWigner.h"
#include "RequiresMeasuredBreakupMomenta.h"

#include <memory>

namespace yap {

/// \class RelativisticBreitWigner
/// \brief Class for Relativistic Breit-Wigner resonance shape
/// \author Daniel Greenwald
/// \ingroup MassShapes
///
/// Amplitude is\n
/// \f$\frac{1}{M_{R}^2 - s - i*M_{R}*\Gamma(s)}\f$\n
/// with\n
/// \f$\Gamma(s) = \Gamma_{R}
///                \left(\frac{p^*(s)}{p^*(M_{R}^2)}\right)^{2J_{R}+1}
///                \frac{M_{R}}{\sqrt{s}} F^{2}_{R}\f$\n
/// with \f$ F^{2}_{R} \f$ is the Blatt-Weisskopf barrier factor
class RelativisticBreitWigner :
    public BreitWigner,
    public RequiresMeasuredBreakupMomenta

{
public:

    /// Constructor
    /// \param m Mass of resonance [GeV]
    /// \param w Width of resonance [GeV]
    RelativisticBreitWigner(double m = -1, double w = -1) :
    BreitWigner(m, w), RequiresMeasuredBreakupMomenta(true) {}

    /// Check if a DecayChannel is valid for this MassShape; will throw if invalid.
    /// Cheks that decay is to two spin-zero particles
    virtual void checkDecayChannel(const DecayChannel& c) const override;

protected:

    /// Retrieve BlattWeisskopf object from Resonance now that it is added to the Model
    virtual void addDecayChannel(std::shared_ptr<DecayChannel> c) override;

    /// Calculate dynamic amplitude T for and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculateT(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc, unsigned si) const override;

private:

    /// BlattWeisskopf object needed for width calculation
    std::shared_ptr<BlattWeisskopf> BlattWeisskopf_;

};

}

#endif
