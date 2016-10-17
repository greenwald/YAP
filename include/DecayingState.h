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

#ifndef yap_DecayingState_h
#define yap_DecayingState_h

#include "fwd/DecayingState.h"

#include "fwd/BlattWeisskopf.h"
#include "fwd/Parameter.h"
#include "fwd/QuantumNumbers.h"

#include "Particle.h"

#include <memory>
#include <string>

namespace yap {

/// \class DecayingState
/// \brief Represents a named decaying state with quantum numbers and a radial size
/// \ingroup Particle
/// \author Daniel Greenwald
class DecayingState : public Particle
{
protected:
    /// Constructor
    DecayingState(const std::string& name, const QuantumNumbers& q, double radial_size);

public:
    /// \return radial size [GeV^-1]
    std::shared_ptr<PositiveRealParameter>& radialSize()
    { return RadialSize_; }

    /// \return radial size [GeV^-1]
    const std::shared_ptr<PositiveRealParameter>& radialSize() const
    { return RadialSize_; }

protected:

    using Particle::addParticleCombination;

    /// add ParticleCombination to BlattWeisskopf object
    void addParticleCombination(BlattWeisskopf& bw, const ParticleCombination& pc) const;

    using Particle::registerWithModel;

    /// call BlattWeisskopf's registerWithModel
    void registerWithModel(BlattWeisskopf& bw) const;

private:
    /// Radial size parameter [GeV^-1]
    std::shared_ptr<PositiveRealParameter> RadialSize_;
};

}

#endif
