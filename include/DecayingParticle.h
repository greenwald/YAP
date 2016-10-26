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

#ifndef yap_DecayingParticle_h
#define yap_DecayingParticle_h

#include "fwd/DecayChannel.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/QuantumNumbers.h"

#include "DecayingState.h"

#include <memory>
#include <string>

namespace yap {

/// Class for a particle that will decay
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle
///
/// The amplitude function returns a sum over the amplitudes of all
/// #DecayChannel's for the decay of the particle (denoted P, with
/// daughters in channel c denoted D1, D2; and amplitude A_c):\n
/// A_c = a_c * Blatt-Weisskopf(P->D1+D2) * SpinAmplitude(P->D1+D2) * A(D1->xx) * A(D2->xx)\n
/// with free amplitude a_c.
class DecayingParticle : public DecayingState
{
protected:

    /// Constructor
    /// see #create
    DecayingParticle(const std::string& name, const QuantumNumbers& q, double radial_size)
        : DecayingState(name, q, radial_size) {}

public:

    /// \note inherits addChannel
    using DecayingState::addChannel;
    
    /// create
    /// \param name Name of decaying particle
    /// \param q QuantumNumbers of decaying particle
    /// \param radialSize radial size of decaying particle
    static std::shared_ptr<DecayingParticle> create(const std::string& name, const QuantumNumbers& q, double radialSize)
    { return std::shared_ptr<DecayingParticle>(new DecayingParticle(name, q, radialSize)); }

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param daughters ParticleVector of daughters to create DecayChannel object from
    /// \return shared_ptr to DecayChannel that has been added
    std::shared_ptr<DecayChannel> addChannel(const ParticleVector& daughters);

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param A shared_ptr to a daughter
    /// \param B shared_ptr to a daughter
    /// \param other_daughters... other daughters
    /// \return shared_ptr to DecayChannel that has been added
    template <typename ... Types>
    std::shared_ptr<DecayChannel> addChannel(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, Types ... other_daughters)
    { ParticleVector V{A, B, other_daughters...}; return addChannel(V); }

    /// grant friend status to Model to call fixSolitaryFreeAmplitudes()
    friend Model;

protected:

    /// if only one decay channel is available, fix its free amplitude to the current value
    void fixSolitaryFreeAmplitudes();

};

/// \return whether Particle is a DecayingParticle
inline bool is_decaying_particle(const Particle& p)
{ return dynamic_cast<const Particle*>(&p) != nullptr; }

}

#endif
