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

#ifndef yap_DecayChannel_h
#define yap_DecayChannel_h

#include "fwd/DecayChannel.h"

#include "fwd/DecayingState.h"
#include "fwd/Model.h"
#include "fwd/Particle.h"
#include "fwd/ParticleCombination.h"
#include "fwd/Spin.h"
#include "fwd/SpinAmplitude.h"

#include <memory>
#include <string>

namespace yap {

/// \class DecayChannel
/// \brief Class implementing a decay channel.
/// \author Johannes Rauch, Daniel Greenwald
class DecayChannel
{
public:

    /// \name Constructors
    /// @{

    /// N-particle Constructor (at the moment only valid for 2 particles).
    /// DecayChannel inherits ISP from daughters.
    /// \param daughters Vector of shared_ptr's to daughter Particle's
    DecayChannel(const ParticleVector& daughters);

    /// @}

    /// check consistency of object
    virtual bool consistent() const;

    /// \name Getters
    /// @{

    /// Get Daughters
    const ParticleVector& daughters() const
    { return Daughters_; }

    /// Get SpinAmplitude objects
    const SpinAmplitudeVector& spinAmplitudes() const
    { return SpinAmplitudes_; }

    /// \return vector of ParticleCombinations
    const ParticleCombinationSet& particleCombinations() const
    { return ParticleCombinations_; }

    /// \return raw pointer to model through first Daughter
    const Model* model() const;

    /// @}

    /// add a spin amplitude
    void addSpinAmplitude(std::shared_ptr<SpinAmplitude> sa);

    /// automatically create all possible spin amplitudes given initial spin J
    /// \param two_J, (twice) decaying particle's spin
    void addAllPossibleSpinAmplitudes(unsigned two_J);

    /// grant friend status to DecayingState to call
    /// addParticleCombination, pruneParticleCombinations,
    /// and registerWithModel
    friend class DecayingState;

protected:

    /// register any necessary DataAccessor's with model
    virtual void registerWithModel();

    /// Add particle combination
    virtual void addParticleCombination(const ParticleCombination& c);

    /// prune ParticleCombinations_ to only contain ParticleCombination's tracing back up the ISP
    virtual void pruneParticleCombinations();

private:

    /// daughters of the decay
    ParticleVector Daughters_;

    /// Vector of SpinAmplitudes
    SpinAmplitudeVector SpinAmplitudes_;

    /// vector of shared_ptr<ParticleCombination>
    ParticleCombinationSet ParticleCombinations_;

};

/// convert to string
std::string to_string(const DecayChannel& dc);

/// convert to string
std::string to_string(const DecayChannel& dc, int two_M, const SpinProjectionVector& two_m);

/// << operator
inline std::ostream& operator<< (std::ostream& os, const DecayChannel& dc)
{ os << to_string(dc); return os; }

}

#endif
