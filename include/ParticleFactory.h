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

#ifndef yap_ParticleFactory_h
#define yap_ParticleFactory_h

#include "fwd/DecayingParticle.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/MassShape.h"
#include "fwd/ParticleFactory.h"
#include "fwd/Resonance.h"

#include "QuantumNumbers.h"

#include <memory>
#include <string>
#include <vector>

namespace yap {

/// \struct ParticleTableEntry
/// \brief Data container for storing particle information in database
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup ParticleFactory
struct ParticleTableEntry : public QuantumNumbers {
    /// constructor
    /// \param pdg PDG code
    /// \param name Particle name
    /// \param q QuantumNumbers of particle
    /// \param mass Mass of particle
    /// \param parameters Further parameters of particle (implementation dependent)
    ParticleTableEntry(int pdg = 0, std::string name = "", QuantumNumbers q = QuantumNumbers(), double mass = -1, std::vector<double> parameters = {});

    /// \return consistency of entry
    bool consistent() const override;

    /// PDG code of particle
    int PDG;

    /// Name of particle
    std::string Name;

    /// Mass of particle
    double Mass;

    /// further parameters of particle (implementation dependent)
    std::vector<double> MassShapeParameters;
};

/// \class ParticleFactory
/// \brief Factory class for easy creation of Particle objects from PDG codes.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle
/// \defgroup ParticleFactory
class ParticleFactory
{
public:

    /// \typedef ParticleFactory::value_type
    /// Define this to allow `std::inserter` to use `insert`
    using value_type = ParticleTableEntry;

    /// \typedef ParticleFactory::iterator
    /// Define this to allow `std::inserter` to use `insert`
    using iterator = ParticleTableMap::iterator;

    /// Create a FinalStateParticle from a PDG code
    /// \param PDG PDG code of particle to create
    /// \return shared pointer to new final state particle
    std::shared_ptr<FinalStateParticle> fsp(int PDG) const;

    /// Create an decayingParticle from a PDG code
    /// \param PDG PDG code of particle to create
    /// \param radialSize radial size of particle to create [GeV^-1]
    /// \return shared pointer to new DecayingParticle object
    std::shared_ptr<DecayingParticle> decayingParticle(int PDG, double radialSize) const;

    /// Create a Resonance from a PDG code and a MassShape
    /// \param PDG PDG code of particle to create
    /// \param radialSize Radial size of particle to create [GeV^-1]
    /// \param massShape Pointer to MassShape object describing resonance
    /// \return shared pointer to new Resonance object
    std::shared_ptr<Resonance> resonance(int PDG, double radialSize, std::shared_ptr<MassShape> massShape) const;

    /// Adds content of rhs to this
    /// \param rhs ParticleFactory to add into this
    ParticleFactory& operator+=(const ParticleFactory& rhs);

    /// \name Particle table access
    /// @{

    /// get ParticleTableEntry from #ParticleTable_ with safety checks
    /// \param PDG pdg code labeling particle table entry
    const ParticleTableEntry& operator[](int PDG) const;

    /// get ParticleTableEntry from #ParticleTable_ with safety checks
    /// \param name Name of particle in table
    const ParticleTableEntry& operator[](std::string name) const
    { return (*this)[pdgCode(name)]; }

    /// get #QuantumNumbers from #ParticleTable_ with safety checks
    /// \param PDG pdg code labeling particle table entry
    const QuantumNumbers& quantumNumbers(int PDG) const
    { return static_cast<const QuantumNumbers&>((*this)[PDG]); }

    /// get #QuantumNumbers from #ParticleTable_ with safety checks
    /// \param name Name of particle in table
    const QuantumNumbers& quantumNumbers(std::string name) const
    { return static_cast<const QuantumNumbers&>((*this)[name]); }

    /// inserts the pair `ParticleTableEntry::PDG` and `ParticleTableEntry` to #ParticleTable_
    /// \param entry a A ParticleTableEntry to add to #ParticleTable_
    std::pair<ParticleTableMap::iterator, bool> insert(const ParticleTableEntry& entry);

    /// convenience function to allow `inserter()` to be used in the `std::copy` algorithm
    ParticleTableMap::iterator insert(ParticleTableMap::iterator hint, const ParticleTableEntry& entry);

    /// #ParticleFactory's own inserter
    friend std::insert_iterator<ParticleFactory> inserter(ParticleFactory& F)
    { return std::insert_iterator<ParticleFactory>(F, F.ParticleTable_.end()); }

    // find PDG number by particle name
    // \return PDG code number
    // \param name Particle name as listed in particle table
    int pdgCode(std::string name) const;

    /// @}

private:

    /// maps PDGCodes to ParticleTableEntry's
    ParticleTableMap ParticleTable_;
};

}

#endif
