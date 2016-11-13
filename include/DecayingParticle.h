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

#include "fwd/DecayingParticle.h"

#include "fwd/DecayChannel.h"
#include "fwd/Particle.h"
#include "fwd/ParticleFactory.h"
#include "fwd/QuantumNumbers.h"

#include "AttributeUtilities.h"
#include "DecayingState.h"

#include <memory>

namespace yap {

/// A DecayingState with one or more possible DecayChannels
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle
class DecayingParticle : public DecayingState
{
protected:

    /// Constructor
    /// see #create
    DecayingParticle(const std::string& name, const QuantumNumbers& q, double radial_size, std::shared_ptr<MassShape> mass_shape)
        : DecayingState(name, q, radial_size, mass_shape) {}

    /// Constructor
    /// see #create
    DecayingParticle(const ParticleTableEntry& pte, double radial_size, std::shared_ptr<MassShape> mass_shape)
        : DecayingState(pte, radial_size, mass_shape) {}
    
public:

    /// create
    /// \param name Name of decaying particle
    /// \param q QuantumNumbers of decaying particle
    /// \param radialSize radial size of decaying particle
    /// \param mass_shape Dynamic ampltiude object
    static std::shared_ptr<DecayingParticle> create(const std::string& name, const QuantumNumbers& q, double radial_size, std::shared_ptr<MassShape> mass_shape = nullptr)
    { return std::shared_ptr<DecayingParticle>(new DecayingParticle(name, q, radial_size, mass_shape)); }

    /// create
    /// \param pte ParticleTableEntry
    /// \param radialSize radial size of decaying particle
    /// \param mass_shape Dynamic ampltiude object
    static std::shared_ptr<DecayingParticle> create(const ParticleTableEntry& pte, double radial_size, std::shared_ptr<MassShape> mass_shape = nullptr)
    { return std::shared_ptr<DecayingParticle>(new DecayingParticle(pte, radial_size, mass_shape)); }

    /// automatically create all possible spin amplitudes
    /// \param dc DecayChannel to add to
    /// \param conserve_parity whether to conserve parity
    void addAllPossibleSpinAmplitudes(DecayChannel& dc, bool conserve_parity) const;
    
    /// Create a DecayChannel and add it to this DecayingParticle
    /// \param daughters ParticleVector of daughters to create DecayChannel object from
    /// \param conserve_parity whether to conserve parity in decay, when adding spin amplitudes automatically
    /// \return shared_ptr to DecayChannel that has been added
    std::shared_ptr<DecayChannel> addDecay(const ParticleVector& daughters, bool conserve_parity);

    /// Create and add a DecayChannel, creating SpinAmplitudes with parity conservation ignored
    /// \param daughters ParticleVector of daughters to create DecayChannel object from
    /// \return shared_ptr to DecayChannel that has been added
    std::shared_ptr<DecayChannel> addWeakDecay(const ParticleVector& daughters)
    { return addDecay(daughters, false); }

    /// Create and add a DecayChannel, creating SpinAmplitudes with parity conserved
    /// \param daughters ParticleVector of daughters to create DecayChannel object from
    /// \return shared_ptr to DecayChannel that has been added
    std::shared_ptr<DecayChannel> addStrongDecay(const ParticleVector& daughters)
    { return addDecay(daughters, true); }

    /// Create and add a DecayChannel, creating SpinAmplitudes with parity conservation ignored
    /// \param A shared_ptr to a daughter
    /// \param B shared_ptr to a daughter
    /// \param other_daughters... other daughters
    /// \return shared_ptr to DecayChannel that has been added
    template <typename ... Types>
    std::shared_ptr<DecayChannel> addWeakDecay(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, Types ... other_daughters)
    { ParticleVector V{A, B, other_daughters...}; return addWeakDecay(V); }

    /// Create and add a DecayChannel, creating SpinAmplitudes with parity conserved
    /// \param A shared_ptr to a daughter
    /// \param B shared_ptr to a daughter
    /// \param other_daughters... other daughters
    /// \return shared_ptr to DecayChannel that has been added
    template <typename ... Types>
    std::shared_ptr<DecayChannel> addStrongDecay(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, Types ... other_daughters)
    { ParticleVector V{A, B, other_daughters...}; return addStrongDecay(V); }

};

/// checks if something inherits from DecayingParticle
extern const is_of_type<DecayingParticle> is_decaying_particle;
 
}

#endif
