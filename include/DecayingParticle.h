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

#include "fwd/BlattWeisskopf.h"
#include "fwd/DecayChannel.h"
#include "fwd/DecayTree.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/QuantumNumbers.h"

#include "Particle.h"

#include <memory>

namespace yap {

/// \class DecayingParticle
/// \brief Class for a particle that will decay
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
    DecayingParticle(const QuantumNumbers& q, std::string name, double radialSize);

public:

    /// create
    /// \param q QuantumNumbers of decaying particle
    /// \param radialSize radial size of decaying particle
    static std::shared_ptr<DecayingParticle> create(const QuantumNumbers& q, std::string name, double radialSize)
    { return std::shared_ptr<DecayingParticle>(new DecayingParticle(q, name, radialSize)); }

    /// \return DecayTrees
    /// map key is spin projection
    const DecayTreeVectorMap& decayTrees() const
    { return DecayTrees_; }

    /// Check consistency of object
    virtual bool consistent() const override;

    /// Check if a DecayChannel is valid for DecayingParticle;
    /// will throw if invalid
    virtual void checkDecayChannel(const DecayChannel& c) const;

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param c unique_ptr to DecayChannel, should be constructed in function call, or use std::move(c)
    /// \return shared_ptr to DecayChannel that has been added
    virtual std::shared_ptr<DecayChannel> addChannel(std::shared_ptr<DecayChannel> c);

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

    /// \name Getters
    /// @{

    /// \return channels
    const DecayChannelVector& channels() const
    { return Channels_;}

    /// \return Radial size [GeV^-1]
    std::shared_ptr<RealParameter> radialSize()
    { return RadialSize_; }

    /// \return Blatt-Weisskopf factors
    const BlattWeisskopfMap& blattWeisskopfs() const
    { return BlattWeisskopfs_; }

    /// @}

    /// Print complete decay chain
    void printDecayChain() const
    { printDecayChainLevel(0); }

    /// \return raw pointer to Model through first DecayChannel
    const Model* model() const override;

    /// grant friend status to Model to call fixSolitaryFreeAmplitudes()
    friend Model;

    /// grant friend status to DecayChannel to call registerWithModel()
    friend DecayChannel;

protected:

    /// add ParticleCombination to SymmetrizationIndices_ and BlattWeisskopfs_
    virtual void addParticleCombination(const std::shared_ptr<ParticleCombination>& c) override;

    /// prune ParticleCombinations_ to only contain ParticleCombination's tracing back up the ISP
    virtual void pruneParticleCombinations() override;

    /// register any necessary DataAccessor's with model
    virtual void registerWithModel() override;

    /// if only one decay channel is available, fix its free amplitude to the current value
    void fixSolitaryFreeAmplitudes();

    void printDecayChainLevel(int level) const;

    /// modify a DecayTree
    /// \param dt DecayTree to modify
    virtual void modifyDecayTree(DecayTree& dt) const;

private:

    /// vector of decay channel objects
    DecayChannelVector Channels_;

    /// map of Blatt-Weisskopf barrier factors, key = angular momentum
    BlattWeisskopfMap BlattWeisskopfs_;

    /// Radial size parameter [GeV^-1]
    std::shared_ptr<RealParameter> RadialSize_;

    /// Map of spin projection to DecayTreeVector
    DecayTreeVectorMap DecayTrees_;

};

/// convert to (multiline) string
std::string to_string(const DecayTreeVectorMap& m_dtv_map);

/// \return Set of all decay trees in provided DecayingParticle
/// \todo Have it recursively travel down DecayChannels?
DecayTreeSet decay_trees(const DecayingParticle& dp);

/// \return DecayTreeSet for decays passing provided predicates
/// \param p last predicate to apply in filtering DecayTree's
/// \param P... predicates to apply in filtering DecayTree's
template <typename Last, typename ... Predicates>
DecayTreeSet decay_trees(const DecayingParticle& dp, Last p, Predicates ... P)
{ return filter(decay_trees(dp), p, P...); }

/// \return lone DecayTree passing provided predicates
/// \param p last predicate to apply in filtering DecayTree's
/// \param P... predicates to apply in filtering DecayTree's
/// throws if no unique DecayTree is found
template <typename Last, typename ... Predicates>
DecayTreeSet::value_type decay_tree(const DecayingParticle& dp, Last p, Predicates ... P)
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
