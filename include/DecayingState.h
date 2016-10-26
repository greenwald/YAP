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
#include "fwd/DecayChannel.h"
#include "fwd/DecayTree.h"
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

    /// \return channels
    const DecayChannelVector& channels() const
    { return DecayChannels_;}

    /// \return Blatt-Weisskopf factors
    const BlattWeisskopfMap& blattWeisskopfs() const
    { return BlattWeisskopfs_; }

    /// \return DecayTrees
    /// map key is spin projection
    const DecayTreeVector& decayTrees() const
    { return DecayTrees_; }

    /// \return raw pointer to Model through first DecayChannel
    const Model* model() const override;

    /// Check if a DecayChannel is valid for DecayingState;
    /// will throw if invalid
    virtual void checkDecayChannel(const DecayChannel& c) const;

    /// check consistency of object
    virtual bool consistent() const override;
    
protected:

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param c unique_ptr to DecayChannel, should be constructed in function call, or use std::move(c)
    /// \return shared_ptr to DecayChannel that has been added
    virtual std::shared_ptr<DecayChannel> addChannel(std::shared_ptr<DecayChannel> c);

    /// calls Particle::addParticleCombination and then calls on BlattWeisskopf's and DecayChannel's
    virtual void addParticleCombination(const ParticleCombination& c) override;

    /// calls registerWithModel of BlattWeisskopf's and DecayChannel's
    virtual void registerWithModel() override;
    
    /// calls Particle::pruneParticleCombinations, then calls on each DecayChannel
    virtual void pruneParticleCombinations() override;
    
    /// modify a DecayTree
    /// \param dt DecayTree to modify
    virtual void modifyDecayTree(DecayTree& dt) const;

private:

    /// Radial size parameter [GeV^-1]
    std::shared_ptr<PositiveRealParameter> RadialSize_;

    /// vector of decay channel objects
    DecayChannelVector DecayChannels_;

    /// map of Blatt-Weisskopf barrier factors, key = angular momentum
    BlattWeisskopfMap BlattWeisskopfs_;

    /// vector of DecayTrees
    DecayTreeVector DecayTrees_;

};

/// \return whether Particle is a DecayingState
inline bool is_decaying_state(const Particle& p)
{ return dynamic_cast<const DecayingState*>(&p) != nullptr; }

/// \return multiline string displaying decay
std::string to_decay_string(const DecayingState& ds, unsigned level = 0);

/// \return set of all particles below given DecayingState, including itself
ParticleSet particles(DecayingState& ds);

/// \return lone DecayTree passing provided predicates
/// \param p last predicate to apply in filtering DecayTree's
/// \param P... predicates to apply in filtering DecayTree's
/// throws if no unique DecayTree is found
template <typename Last, typename ... Predicates>
DecayTreeSet::value_type decay_tree(const DecayingState& ds, Last p, Predicates ... P)
{ return lone_elt(filter(ds.decayTrees(), p, P...)); }

/// \return all the free amplitudes under a DecayingState
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
