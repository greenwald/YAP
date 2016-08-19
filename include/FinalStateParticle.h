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

#include "Constants.h"
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
    FinalStateParticle(const QuantumNumbers& q, double m, std::string name)
        : Particle(q, name), Model_(nullptr), Mass_(m) {}

public:

    /// create
    /// \param q Quantum numbers of particle
    /// \param m Mass of particle
    /// \param name Name of particle
    static std::shared_ptr<FinalStateParticle> create(const QuantumNumbers& q, double m, std::string name)
    { return std::shared_ptr<FinalStateParticle>(new FinalStateParticle(q, m, name)); }

    /// check mass is positive semi-definite
    bool consistent() const;

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
    virtual void addParticleCombination(const std::shared_ptr<ParticleCombination>& pc) override;

    /// register any necessary DataAccessor's with model
    virtual void registerWithModel() override
    { }

private:

    /// raw pointer to Model owning this final state particle
    Model* Model_;

    /// Mass [GeV]
    double Mass_;

};

/// \return whether pc and FSPs are for the same final state
bool valid_final_state(const std::shared_ptr<ParticleCombination>& pc, const FinalStateParticleVector& FSPs);

}

#endif
