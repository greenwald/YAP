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

#ifndef yap_FinalStateParticle_h
#define yap_FinalStateParticle_h

#include "fwd/FinalStateParticle.h"

#include "fwd/DataPoint.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/ParticleFactory.h"

#include "AttributeUtilities.h"
#include "Particle.h"

#include <complex>
#include <memory>

namespace yap {

/// \class FinalStateParticle
/// \brief Class representing a final-state particle
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle

class FinalStateParticle : public Particle
{
protected:

    /// Constructor
    /// see #create
    FinalStateParticle(const std::string& name, const QuantumNumbers& q, double m);

    /// Constructor
    /// see #create
    FinalStateParticle(const ParticleTableEntry& pte);

public:

    /// create
    /// \param name Name of particle
    /// \param q Quantum numbers of particle
    /// \param m Mass of particle
    static std::shared_ptr<FinalStateParticle> create(const std::string& name, const QuantumNumbers& q, double m)
    { return std::shared_ptr<FinalStateParticle>(new FinalStateParticle(name, q, m)); }

    /// create
    /// \param pte ParticleTableEntry
    static std::shared_ptr<FinalStateParticle> create(const ParticleTableEntry& pte)
    { return std::shared_ptr<FinalStateParticle>(new FinalStateParticle(pte)); }
    
    /// Get mass [GeV]
    double mass() const
    { return Mass_; }

    /// \return owning Model
    const Model* model() const override
    { return Model_; }

    /// Grant friend status to Model to set FSP's indices
    friend class Model;

protected:

    /// set raw pointer to owning Model
    void setModel(Model* m)
    { Model_ = m; }

    /// add ParticleCombination to ParticleCombinations
    virtual void addParticleCombination(const ParticleCombination& pc) override;

    /// register any necessary DataAccessor's with model
    virtual void registerWithModel() override
    { }

private:

    /// raw pointer to Model owning this final state particle
    Model* Model_;

    /// Mass [GeV]
    double Mass_;

};

/// checks if something inherits from FinalStateParticle
extern const is_of_type<FinalStateParticle> is_final_state_particle;

}

#endif
