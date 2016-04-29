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

#include "fwd/DataPointBase.h"

#include "DataAccessor.h"
#include "DecayChannel.h"
#include "DecayTree.h"
#include "make_unique.h"
#include "Particle.h"

#include <complex>
#include <map>
#include <memory>
#include <vector>

namespace yap {

class BlattWeisskopf;
class FinalStateParticle;
class Model;
class ParticleCombination;
class QuantumNumbers;
class StatusManager;

/// \class DecayingParticle
/// \brief Class for a particle that will decay
/// \authors Johannes Rauch, Daniel Greenwald
/// \ingroup Particle
///
/// The amplitude function returns a sum over the amplitudes of all
/// #DecayChannel's for the decay of the particle (denoted P, with
/// daughters in channel c denoted D1, D2; and amplitude A_c):\n
/// A_c = a_c * Blatt-Weisskopf(P->D1+D2) * SpinAmplitude(P->D1+D2) * A(D1->xx) * A(D2->xx)\n
/// with free amplitude a_c.

class DecayingParticle : public Particle, public DataAccessor
{
protected:

    /// Constructor
    /// see #create
    DecayingParticle(const QuantumNumbers& q, double mass, std::string name, double radialSize);

public:

    /// create
    /// \param q QuantumNumbers of decaying particle
    /// \param mass mass of decaying particle
    /// \param radialSize radial size of decaying particle
    static std::shared_ptr<DecayingParticle> create(const QuantumNumbers& q, double mass, std::string name, double radialSize)
    { return std::shared_ptr<DecayingParticle>(new DecayingParticle(q, mass, name, radialSize)); }

    /// Calculate complex amplitude
    /// \param d DataPointBase to calculate with
    /// \param pc (shared_ptr to) ParticleCombination to calculate for
    /// \param two_m 2 * the spin projection to calculate for
    /// \param sm StatusManager to update
    virtual std::complex<double> amplitude(DataPointBase& d, const std::shared_ptr<ParticleCombination>& pc, int two_m, StatusManager& sm) const override;

    /// Check consistency of object
    virtual bool consistent() const override;

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param c unique_ptr to DecayChannel, should be constructed in function call, or use std::move(c)
    /// \return shared_ptr to DecayChannel that has been added
    virtual std::shared_ptr<DecayChannel> addChannel(std::shared_ptr<DecayChannel> c);

    /// Add a DecayChannel and set its parent to this DecayingParticle.
    /// \param daughters ParticleVector of daughters to create DecayChannel object from
    /// \return shared_ptr to DecayChannel that has been added
    virtual std::shared_ptr<DecayChannel> addChannel(const ParticleVector& daughters)
    { return addChannel(std::make_shared<DecayChannel>(daughters)); }

    /// Return final state particles of a channel (vector should be identical for all channels)
    /// \return vector of shared_ptr's to FinalStateParticles of this decaying particle (in channel i)
    /// \param i index of DecayChannel to return FinalStateParticles of.
    std::vector< std::shared_ptr<FinalStateParticle> > finalStateParticles(unsigned i = 0) const;

    /// \name Getters
    /// @{

    /// \return channel
    std::shared_ptr<DecayChannel> channel(const ParticleVector& daughters);

    /// \return channels
    const DecayChannelVector& channels() const
    { return Channels_;}

    /// \return Number of decay channels for this object
    unsigned int nChannels() const
    { return Channels_.size(); }

    /// Return Channel i
    std::shared_ptr<DecayChannel> channel(unsigned i) const
    { return Channels_.at(i); }

    /// Return Channel i
    std::shared_ptr<DecayChannel> channel(unsigned i)
    { return Channels_.at(i); }

    /// \return Radial size [GeV^-1]
    std::shared_ptr<RealParameter> radialSize()
    { return RadialSize_; }

    /// @}

    /// Print complete decay chain
    void printDecayChain() const
    { printDecayChainLevel(0); }

    virtual CachedDataValueSet cachedDataValuesItDependsOn() override;

    /// \return raw pointer to Model through first DecayChannel
    const Model* model() const override
    { return Channels_.empty() ? nullptr : Channels_[0]->model(); }

    /// \return string denoting DataAccessor type
    virtual std::string data_accessor_type() const
    { return "DecayingParticle"; }

    /// grant friend status to DecayChannel to see BlattWeiskopffs_
    friend DecayChannel;

    /// grant friend status to Model to see freeAmplitudes
    friend class Model;

protected:

    /// add ParticleCombination to SymmetrizationIndices_ and BlattWeisskopfs_
    virtual unsigned addParticleCombination(std::shared_ptr<ParticleCombination> c) override;

    void printDecayChainLevel(int level) const;

    /// \return vector of shared_ptr's to all free amplitudes from this point in decay tree and down
    ComplexParameterVector freeAmplitudes() const;

    /// \return vector od DecayTree's from this point in decay tree and down
    DecayTreeVector decayTrees() const;

private:

    /// vector of decay channel objects
    DecayChannelVector Channels_;

    /// map of Blatt-Weisskopf barrier factors, key = angular momentum
    std::map<unsigned, std::shared_ptr<BlattWeisskopf> > BlattWeisskopfs_;

    /// Radial size parameter [GeV^-1]
    std::shared_ptr<RealParameter> RadialSize_;

    /// Cached amplitudes for each spin projection
    /// key = 2 * spin projection
    std::map<int, std::shared_ptr<ComplexCachedDataValue> > Amplitudes_;

};

}

#endif
