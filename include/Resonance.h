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

#include "fwd/DataPointBase.h"

#include "DecayingParticle.h"

#include <complex>
#include <memory>
#include <string>

namespace yap {

class MassShape;
class ParticleCombination;
class QuantumNumbers;
class StatusManager;

/// \class Resonance
/// \brief Class for a particle that will decay and has a mass shape
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class Resonance : public DecayingParticle
{
protected:

    /// Constructor
    /// see #create
    Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, std::shared_ptr<MassShape> massShape);

public:

    /// create
    /// \param q QuantumNumbers of resonance
    /// \param mass mass of resonance
    /// \param radialSize radial size of resonance
    /// \param massShape shared_ptr to MassShape of resonance
    static std::shared_ptr<Resonance> create(const QuantumNumbers& q, double mass, std::string name, double radialSize, std::shared_ptr<MassShape> massShape)
    { return std::shared_ptr<Resonance>(new Resonance(q, mass, name, radialSize, massShape)); }

    /// Calculate complex amplitude
    /// \param d DataPointBase to calculate with
    /// \param pc (shared_ptr to) ParticleCombination to calculate for
    /// \param two_m 2 * the spin projection to calculate for
    /// \param sm StatusManager to update
    virtual std::complex<double> amplitude(DataPointBase& d, const std::shared_ptr<ParticleCombination>& pc, int two_m, StatusManager& sm) const override;

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

    virtual std::string data_accessor_type() const override
    {return "Resonance"; }

    /// grant friend status to MassShape to call setMass
    friend MassShape;

protected:

    /// overrides DataAccessor's function to also register MassShape_ with Model
    virtual void addToModel() override;

    /// add ParticleCombination to ParticleCombinations_,
    /// also add to MassShape_
    virtual unsigned addParticleCombination(std::shared_ptr<ParticleCombination> c) override;

private:

    /// MassShape object
    std::shared_ptr<MassShape> MassShape_;

};

}

#endif
