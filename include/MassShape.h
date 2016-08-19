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

#ifndef yap_MassShape_h
#define yap_MassShape_h

#include "fwd/CachedValue.h"
#include "fwd/DataPartition.h"
#include "fwd/DecayChannel.h"
#include "fwd/Model.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"
#include "fwd/ParticleFactory.h"
#include "fwd/Resonance.h"
#include "fwd/StatusManager.h"

#include "RecalculableDataAccessor.h"

#include <complex>
#include <memory>
#include <string>

namespace yap {

/// \class MassShape
/// \brief Abstract base class for all mass shapes
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup MassShapes Mass Shapes
class MassShape :
    public RecalculableDataAccessor
{
public:

    /// Constructor
    MassShape();

    /// \return dynamic amplitude for data point and particle combination
    /// \param d DataPoint
    /// \param pc shared_ptr to ParticleCombination
    virtual std::complex<double> value(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const override final;

    /// Calculate complex amplitudes for and store in each DataPoint in DataPartition;
    /// calls calculateT, which must be overrided in derived classes
    /// \param D DataPartition to calculate on
    virtual void calculate(DataPartition& D) const override final;

    /// Set parameters from ParticleTableEntry
    /// Can be overloaded in inheriting classes
    /// \param entry ParticleTableEntry containing information to create mass shape object
    virtual void setParameters(const ParticleTableEntry& entry)
    { }

    /// Check consistency of object
    virtual bool consistent() const;

    /// get raw pointer to owning resonance
    Resonance* resonance() const
    { return Resonance_; }

    /// update the calculationStatus for a DataPartition
    virtual void updateCalculationStatus(StatusManager& D) const override final;

    /// get raw pointer to Model through resonance
    const Model* model() const override;

    /// Check if a DecayChannel is valid for this MassShape; will throw if invalid.
    virtual void checkDecayChannel(const DecayChannel& c) const
    {}

    /// Grant Resonance friendship, so it can set itself as owner
    /// and call addDecayChannel
    friend class Resonance;

protected:

    /// Set raw pointer to owning Resonance.
    virtual void setResonance(Resonance* r);

    /// Give MassShape chance to perform operations based on the
    /// addition of a DecayChannel to its Resonance
    virtual void addDecayChannel(std::shared_ptr<DecayChannel> c)
    {}

    /// access cached dynamic amplitude
    std::shared_ptr<ComplexCachedValue> T()
    { return T_; }

    /// access cached dynamic amplitude (const)
    const std::shared_ptr<ComplexCachedValue> T() const
    { return const_cast<MassShape*>(this)->T(); }

    /// Calculate dynamic amplitude T for and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculateT(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc, unsigned si) const = 0;

private:

    /// raw pointer to resonance that owns this mass shape
    Resonance* Resonance_;

    /// cached dynamic amplitude
    std::shared_ptr<ComplexCachedValue> T_;

};

}

#endif
