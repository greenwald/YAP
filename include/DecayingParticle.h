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

#include "fwd/AmplitudeComponent.h"
#include "fwd/BlattWeisskopf.h"
#include "fwd/DecayChannel.h"
#include "fwd/DecayTree.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/MassShape.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/QuantumNumbers.h"

#include "AttributeUtilities.h"
#include "Particle.h"

#include <memory>

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
class DecayingParticle : public Particle
{
protected:

    /// Constructor
    /// see #create
    DecayingParticle(const std::string& name, const QuantumNumbers& q, double radial_size, std::shared_ptr<MassShape> mass_shape);

    /// Constructor
    /// see #create
    DecayingParticle(const ParticleTableEntry& pde, double radial_size, std::shared_ptr<MassShape> mass_shape);

public:

    /// create
    /// \param name Name of decaying particle
    /// \param q QuantumNumbers of decaying particle
    /// \param radial_size radial size of decaying particle
    /// \param mass_shape shared_ptr to dynamic amplitude component
    static std::shared_ptr<DecayingParticle> create(const std::string& name, const QuantumNumbers& q, double radial_size, std::shared_ptr<MassShape> mass_shape = nullptr)
    { return std::shared_ptr<DecayingParticle>(new DecayingParticle(name, q, radial_size, mass_shape)); }

    /// create
    /// \param pde ParticleTableEntry to take name and quantum numbers from
    /// \param radial_size radial size of decaying particle
    /// \param mass_shape shared_ptr to dynamic amplitude component
    static std::shared_ptr<DecayingParticle> create(const ParticleTableEntry& pde, double radial_size, std::shared_ptr<MassShape> mass_shape = nullptr)
    { return std::shared_ptr<DecayingParticle>(new DecayingParticle(pde, radial_size, mass_shape)); }

    /// access MassShape
    std::shared_ptr<MassShape> massShape()
    { return MassShape_; }

    /// access MassShape (const)
    std::shared_ptr<const MassShape> massShape() const
    { return MassShape_; }

    /// \return DecayTrees
    const DecayTreeVector& decayTrees() const
    { return DecayTrees_; }

    /// Check consistency of object
    virtual bool consistent() const override;

    /// Add a DecayChannel and set its parent to this DecayingState.
    /// \param c shared_ptr to DecayChannel
    virtual void addDecayChannel(std::shared_ptr<DecayChannel> c);

    /// automatically create all possible spin amplitudes
    /// \param dc DecayChannel to add to
    /// \param conserve_parity whether to conserve parity
    virtual void addAllPossibleSpinAmplitudes(DecayChannel& dc, bool conserve_parity) const;
    
    /// Create a DecayChannel and add it to this DecayingParticle
    /// \param daughters ParticleVector of daughters to create DecayChannel object from
    /// \param conserve_parity whether to conserve parity in decay, when adding spin amplitudes automatically
    /// \return shared_ptr to DecayChannel that has been added
    std::shared_ptr<DecayChannel> addDecay(const ParticleVector& daughters, bool conserve_parity);
    
    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// Parity is _not_ converved
    /// \param daughters ParticleVector of daughters to create DecayChannel object from
    /// \param conserve_parity whether to conserve parity in decay, when adding spin amplitudes automatically
    /// \return shared_ptr to DecayChannel that has been added
    std::shared_ptr<DecayChannel> addWeakDecay(const ParticleVector& daughters)
    { return addDecay(daughters, false); }

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// Parity _is_ converved
    /// \param daughters ParticleVector of daughters to create DecayChannel object from
    /// \param conserve_parity whether to conserve parity in decay, when adding spin amplitudes automatically
    /// \return shared_ptr to DecayChannel that has been added
    std::shared_ptr<DecayChannel> addStrongDecay(const ParticleVector& daughters)
    { return addDecay(daughters, true); }
    
    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// Parity is _not_ converved
    /// \param A shared_ptr to a daughter
    /// \param B shared_ptr to a daughter
    /// \param other_daughters... other daughters
    /// \return shared_ptr to DecayChannel that has been added
    template <typename ... Types>
    std::shared_ptr<DecayChannel> addWeakDecay(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, Types ... other_daughters)
    { ParticleVector V{A, B, other_daughters...}; return addWeakDecay(V); }

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// Parity _is_ conserved
    /// \param A shared_ptr to a daughter
    /// \param B shared_ptr to a daughter
    /// \param other_daughters... other daughters
    /// \return shared_ptr to DecayChannel that has been added
    template <typename ... Types>
    std::shared_ptr<DecayChannel> addStrongDecay(std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, Types ... other_daughters)
    { ParticleVector V{A, B, other_daughters...}; return addStrongDecay(V); }

    /// \name Getters
    /// @{

    /// \return channels
    const DecayChannelVector& channels() const
    { return DecayChannels_;}

    /// \return Radial size [GeV^-1]
    std::shared_ptr<PositiveRealParameter> radialSize()
    { return RadialSize_; }

    /// \return Blatt-Weisskopf factors
    const BlattWeisskopfMap& blattWeisskopfs() const
    { return BlattWeisskopfs_; }

    /// @}

    /// \return raw pointer to Model through first DecayChannel
    const Model* model() const override;

    /// grant friend status to Model to call registerWithModel() and prune()
    friend Model;
    
    /// grant friend status to DecayChannel to call
    /// registerWithModel()
    friend DecayChannel;

    /// grant friend status to MassShape to manipulate DecayTrees_
    friend MassShape;
    
protected:

    /// add ParticleCombination to SymmetrizationIndices_ and BlattWeisskopfs_
    virtual void addParticleCombination(const ParticleCombination& c) override;

    /// calls Particle::prune() and then prunes DecayTrees_ to only
    /// contain trees used by the model; then calls prune on all DecayChannels
    virtual void prune() override;

    /// register any necessary DataAccessor's with model
    virtual void registerWithModel() override;

    /// modify a DecayTree
    /// \param dt DecayTree to modify
    virtual void modifyDecayTree(DecayTree& dt);

    /// add an AmplitudeComponent to a DecayTree
    void addAmplitudeComponent(const AmplitudeComponent& ac, DecayTree& dt) const;

    /// remove an AmplitudeCompoent from a DecayTree
    void removeAmplitudeComponent(const AmplitudeComponent& ac, DecayTree& dt) const;

    /// replace FreeAmplitude in a DecayTree
    void replaceFreeAmplitude(std::shared_ptr<FreeAmplitude> fa, DecayTree& dt) const;
    
private:

    /// MassShape object
    std::shared_ptr<MassShape> MassShape_;
    
    /// vector of decay channel objects
    DecayChannelVector DecayChannels_;

    /// map of Blatt-Weisskopf barrier factors, key = angular momentum
    BlattWeisskopfMap BlattWeisskopfs_;

    /// Radial size parameter [GeV^-1]
    std::shared_ptr<PositiveRealParameter> RadialSize_;

    /// All DecayTree's of this DecayingParticle
    DecayTreeVector DecayTrees_;

};

/// checks if something inherits from DecayingParticle
extern const is_of_type<DecayingParticle> is_decaying_particle;
 
/// convert to (multiline) string
std::string to_decay_string(const DecayingParticle& dp, unsigned level = 0);
 
/// \return Set of all decay trees in provided DecayingParticle
/// \todo Have it recursively travel down DecayChannels?
inline DecayTreeVector decay_trees(const DecayingParticle& dp)
{ return dp.decayTrees(); }

/// \return DecayTreeVector for decays passing provided predicates
/// \param p last predicate to apply in filtering DecayTree's
/// \param P... predicates to apply in filtering DecayTree's
template <typename Last, typename ... Predicates>
DecayTreeVector decay_trees(const DecayingParticle& dp, Last p, Predicates ... P)
{ return filter(decay_trees(dp), p, P...); }

/// \return lone DecayTree passing provided predicates
/// \param p last predicate to apply in filtering DecayTree's
/// \param P... predicates to apply in filtering DecayTree's
/// throws if no unique DecayTree is found
template <typename Last, typename ... Predicates>
DecayTreeVector::value_type decay_tree(const DecayingParticle& dp, Last p, Predicates ... P)
{ return lone_elt(filter(decay_trees(dp), p, P...)); }

/// \return all the free amplitudes under a decaying particle
FreeAmplitudeSet free_amplitudes(const DecayingParticle& dp);

/// \return free amplitude in a model from decay trees evaluating to true
template <typename Last, typename ... UnaryPredicates>
FreeAmplitudeSet free_amplitudes(const DecayingParticle& dp, Last p, UnaryPredicates ... P)
{ return filter(free_amplitudes(dp), p, P...); }

/// \return lone free amplitude in a model from decay trees evaluating to true
/// throws if no unique free amplitude is found
template <typename Last, typename ... UnaryPredicates>
FreeAmplitudeSet::value_type free_amplitude(const DecayingParticle& dp, Last p, UnaryPredicates ... P)
{ return lone_elt(free_amplitudes(dp, p, P...)); }

/// \return set of all particles below given DecayingParticle, including itself
ParticleSet particles(DecayingParticle& dp);

/// \return DecayChannel's of DecayingParticle matching predicates
template <typename Last, typename ... UnaryPredicates>
DecayChannelSet decay_channels(const DecayingParticle& dp, Last p, UnaryPredicates ... P)
{ return filter(DecayChannelSet(dp.channels().begin(), dp.channels().end()), p, P...); }

/// \return lone DecayChannel of DecayingParticle matching predicates
/// throws if no unique DecayChannel is found
template <typename Last, typename ... UnaryPredicates>
DecayChannelSet::value_type decay_channel(const DecayingParticle& dp, Last p, UnaryPredicates ... P)
{ return lone_elt(decay_channels(dp, p, P...)); }

}

#endif
