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

#include "fwd/AmplitudeComponent.h"
#include "fwd/BlattWeisskopf.h"
#include "fwd/DecayChannel.h"
#include "fwd/DecayTree.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/QuantumNumbers.h"

#include "AttributeUtilities.h"
#include "Particle.h"

#include <memory>

namespace yap {

/// Class for a named state with definitive quantum numbers that will decay
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle
///
/// The amplitude function returns a sum over the amplitudes of all
/// #DecayChannel's for the decay of the state (denoted P, with
/// daughters in channel c denoted D1, D2; and amplitude A_c):\n
/// A_c = a_c * Blatt-Weisskopf(P->D1+D2) * SpinAmplitude(P->D1+D2) * A(D1->xx) * A(D2->xx)\n
/// with free amplitude a_c.
class DecayingState : public Particle
{
protected:

    /// Constructor
    /// see #create
    DecayingState(const std::string& name, const QuantumNumbers& q, double radialSize);

public:

    /// create
    /// \param name Name of decaying state
    /// \param q QuantumNumbers of decaying state
    /// \param radialSize radial size of decaying state
    static std::shared_ptr<DecayingState> create(const std::string& name, const QuantumNumbers& q, double radialSize)
    { return std::shared_ptr<DecayingState>(new DecayingState(name, q, radialSize)); }

    /// \return DecayTrees
    /// map key is spin projection
    const DecayTreeVectorMap& decayTrees() const
    { return DecayTrees_; }

    /// Check consistency of object
    virtual bool consistent() const override;

    /// \return channels
    const DecayChannelVector& channels() const
    { return Channels_;}

    /// \return Radial size [GeV^-1]
    std::shared_ptr<RealParameter>& radialSize()
    { return RadialSize_; }

    /// \return Radial size [GeV^-1]
    const std::shared_ptr<RealParameter>& radialSize() const
    { return RadialSize_; }

    /// \return Blatt-Weisskopf factors
    const BlattWeisskopfMap& blattWeisskopfs() const
    { return BlattWeisskopfs_; }

    /// \return raw pointer to Model through first DecayChannel
    const Model* model() const override;

    /// grant friend status to DecayChannel to call registerWithModel()
    friend DecayChannel;

    /// grant friend status to Model to call registerWithModel(),
    /// pruneParticleCombinations()
    friend Model;
    
protected:

    /// automatically create all possible spin amplitudes given initial spin J
    /// parity conservation is ignored if parity is set 0
    /// \param dc DecayChannel to add to
    /// \param conserve_parity whether to conserve parity
    virtual void addAllPossibleSpinAmplitudes(DecayChannel& dc, bool conserve_parity) const;

    /// Add a DecayChannel and set its parent to this DecayingState.
    /// \param c unique_ptr to DecayChannel, should be constructed in function call, or use std::move(c)
    /// \param conserve_parity whether to conserve parity in decay, when adding spin amplitudes automatically
    /// \return shared_ptr to DecayChannel that has been added
    virtual std::shared_ptr<DecayChannel> addDecayChannel(std::shared_ptr<DecayChannel> c, bool conserve_parity = false);

    /// add ParticleCombination to SymmetrizationIndices_ and BlattWeisskopfs_
    virtual void addParticleCombination(const ParticleCombination& c) override;

    /// prune ParticleCombinations_ to only contain ParticleCombination's tracing back up the ISP
    virtual void pruneParticleCombinations() override;

    /// register any necessary DataAccessor's with model
    virtual void registerWithModel() override;

    /// modify a DecayTree
    /// \param dt DecayTree to modify
    virtual void modifyDecayTree(DecayTree& dt);

    /// add an AmplitudeComponent to a DecayTree
    void addAmplitudeComponent(DecayTree& dt, const AmplitudeComponent& da) const;

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

/// checks if something inherits from DecayingState
extern const is_of_type<DecayingState> is_decaying_state;
 
/// convert to (multiline) string
std::string to_decay_string(const DecayingState& ds, unsigned level = 0);
 
/// convert to (multiline) string
std::string to_string(const DecayTreeVectorMap& m_dtv_map);

/// \return Set of all decay trees in provided DecayingState
/// \todo Have it recursively travel down DecayChannels?
DecayTreeSet decay_trees(const DecayingState& ds);

/// \return DecayTreeSet for decays passing provided predicates
/// \param p last predicate to apply in filtering DecayTree's
/// \param P... predicates to apply in filtering DecayTree's
template <typename Last, typename ... Predicates>
DecayTreeSet decay_trees(const DecayingState& ds, Last p, Predicates ... P)
{ return filter(decay_trees(ds), p, P...); }

/// \return lone DecayTree passing provided predicates
/// \param p last predicate to apply in filtering DecayTree's
/// \param P... predicates to apply in filtering DecayTree's
/// throws if no unique DecayTree is found
template <typename Last, typename ... Predicates>
DecayTreeSet::value_type decay_tree(const DecayingState& ds, Last p, Predicates ... P)
{ return lone_elt(filter(decay_trees(ds), p, P...)); }

/// \return all the free amplitudes under a decaying state
FreeAmplitudeSet free_amplitudes(const DecayingState& ds);

/// \return free amplitude in a model from decay trees evaluating to true
template <typename Last, typename ... UnaryPredicates>
FreeAmplitudeSet free_amplitudes(const DecayingState& ds, Last p, UnaryPredicates ... P)
{ return filter(free_amplitudes(ds), p, P...); }

/// \return lone free amplitude in a model from decay trees evaluating to true
/// throws if no unique free amplitude is found
template <typename Last, typename ... UnaryPredicates>
FreeAmplitudeSet::value_type free_amplitude(const DecayingState& ds, Last p, UnaryPredicates ... P)
{ return lone_elt(free_amplitudes(ds, p, P...)); }

/// \return set of all particles below given DecayingState, including itself
ParticleSet particles(const DecayingState& ds);

/// \return DecayChannel's of DecayingState matching predicates
template <typename Last, typename ... UnaryPredicates>
DecayChannelSet decay_channels(const DecayingState& ds, Last p, UnaryPredicates ... P)
{ return filter(DecayChannelSet(ds.channels().begin(), ds.channels().end()), p, P...); }

/// \return lone DecayChannel of DecayingState matching predicates
/// throws if no unique DecayChannel is found
template <typename Last, typename ... UnaryPredicates>
DecayChannelSet::value_type decay_channel(const DecayingState& ds, Last p, UnaryPredicates ... P)
{ return lone_elt(decay_channels(ds, p, P...)); }

}

#endif
