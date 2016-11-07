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

#ifndef yap_Resonance_h
#define yap_Resonance_h

#include "fwd/Resonance.h"

#include "fwd/DataPartition.h"
#include "fwd/DataPoint.h"
#include "fwd/DecayChannel.h"
#include "fwd/MassShape.h"
#include "fwd/ParticleCombination.h"
#include "fwd/QuantumNumbers.h"
#include "fwd/SpinAmplitude.h"
#include "fwd/StatusManager.h"

#include "AttributeUtilities.h"
#include "DecayingParticle.h"

#include <complex>
#include <memory>
#include <string>

namespace yap {

/// \class Resonance
/// \brief Class for a particle that will decay and has a mass shape
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class Resonance : public DecayingParticle
{
protected:

    /// Constructor
    /// see #create
    Resonance(const std::string& name, const QuantumNumbers& q, double radialSize, std::shared_ptr<MassShape> massShape);

public:

    /// create
    /// \param name Name of resonance
    /// \param q QuantumNumbers of resonance
    /// \param radialSize radial size of resonance
    /// \param massShape shared_ptr to MassShape of resonance
    static std::shared_ptr<Resonance> create(const std::string& name, const QuantumNumbers& q, double radialSize, std::shared_ptr<MassShape> massShape)
    { return std::shared_ptr<Resonance>(new Resonance(name, q, radialSize, massShape)); }

    /// Add a DecayChannel to this Resonance
    /// \param c unique_ptr to DecayChannel, should be constructed in function call, or use std::move(c)
    /// \param conserve_parity whether to conserve parity in decay, when adding spin amplitudes automatically
    /// \return shared_ptr to DecayChannel that has been added
    virtual void addDecayChannel(std::shared_ptr<DecayChannel> c) override;

    /// Check consistency of object
    virtual bool consistent() const override;

    /// \name Getters
    /// @{

    /// access MassShape
    std::shared_ptr<MassShape> massShape()
    { return MassShape_; }

    /// access MassShape (const)
    std::shared_ptr<const MassShape> massShape() const
    { return MassShape_; }

    /// @}

protected:

    /// overrides DecayingParticle's function to register MassShape_ with Model
    virtual void registerWithModel() override;

    /// add ParticleCombination to ParticleCombinations_,
    /// also add to MassShape_
    virtual void addParticleCombination(const ParticleCombination& c) override;

    /// modify a DecayTree
    /// \param dt DecayTree to modify
    virtual void modifyDecayTree(DecayTree& dt) override;

private:

    /// MassShape object
    /// \todo Perhaps replace with reference to MassShape
    std::shared_ptr<MassShape> MassShape_;

};

/// checks if something inherits from Resonance
extern const is_of_type<Resonance> is_resonance;

}

#endif
