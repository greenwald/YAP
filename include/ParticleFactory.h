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

#include "QuantumNumbers.h"

#include <memory>
#include <string>
#include <vector>

namespace yap {

/// \struct ParticleTableEntry
/// \brief Data container for storing particle information in database
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup ParticleFactory
class ParticleTableEntry
{
public:
    /// constructor
    /// \param pdg PDG code
    /// \param name Particle name
    /// \param q QuantumNumbers of particle
    /// \param mass Mass of particle
    /// \param parameters Further parameters of particle (implementation dependent)
    ParticleTableEntry(int pdg , const std::string& name , QuantumNumbers q, double mass, std::vector<double> parameters = {});

    /// \return PDG_
    const int pdg() const
    { return PDG_; }

    /// \return Name_
    std::string name() const
    { return Name_; }

    /// \return QuantumNumbers_
    const QuantumNumbers& quantumNumbers() const
    { return QuantumNumbers_; }

    /// \return Mass_
    const double mass() const
    { return Mass_; }

    /// \return MassShapeParameters_
    const std::vector<double>& massShapeParameters() const
    { return MassShapeParameters_; }

    /// set name
    void setName(const std::string& name)
    { Name_ = name; }

    /// set QuantumNumbers_
    void setQuantumNumbers(const QuantumNumbers& qn)
    { QuantumNumbers_ = qn; }

    /// set Mass_
    void mass(double m)
    { Mass_ = m; }

    /// set MassShapeParameters_
    void setMassShapeParameter(const std::vector<double>& pars)
    { MassShapeParameters_ = pars; }

private:
    /// PDG code of particle
    int PDG_;

    /// Name of particle
    std::string Name_;

    /// Quantum numbers of particle
    QuantumNumbers QuantumNumbers_;
    
    /// Mass of particle
    double Mass_;

    /// further parameters of particle (implementation dependent)
    std::vector<double> MassShapeParameters_;
};

/// retrieves element from MassShapeParameters_ and throws with
/// appropriate message if cannot be retrieved.
/// Useful for delegating constructors in MassShape's
/// \param pte ParticleTableEntry to retrieve from
/// \param n elt in vector to retrieve
/// \param where function name to pass to exception thrown if MassShapeParameters_ does not contain n elts
/// \return n-th element in MassShapeParameters_ of entry
const double get_nth_entry(const ParticleTableEntry& pte, size_t n, const std::string& where);

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

    /// Adds content of rhs to this
    /// \param rhs ParticleFactory to add into this
    ParticleFactory& operator+=(const ParticleFactory& rhs);

    /// \name Particle table access
    /// @{

    /// get ParticleTableEntry from #ParticleTable_ with safety checks
    /// \param PDG pdg code labeling particle table entry
    const ParticleTableEntry& operator[](int PDG) const;

    /// get ParticleTableEntry from #ParticleTable_ with safety checks
    /// \param PDG pdg code labeling particle table entry
    ParticleTableEntry& operator[](int PDG)
    { return const_cast<ParticleTableEntry&>(static_cast<const ParticleFactory*>(this)->operator[](PDG)); }

    /// get ParticleTableEntry from #ParticleTable_ with safety checks
    /// \param name Name of particle in table
    const ParticleTableEntry& operator[](const std::string& name) const;

    /// get ParticleTableEntry from #ParticleTable_ with safety checks
    /// \param name Name of particle in table
    ParticleTableEntry& operator[](const std::string& name)
    { return const_cast<ParticleTableEntry&>(static_cast<const ParticleFactory*>(this)->operator[](name)); }

    /// inserts the pair `ParticleTableEntry::PDG` and `ParticleTableEntry` to #ParticleTable_
    /// \param entry a A ParticleTableEntry to add to #ParticleTable_
    std::pair<ParticleTableMap::iterator, bool> insert(const ParticleTableEntry& entry);

    /// convenience function to allow `inserter()` to be used in the `std::copy` algorithm
    ParticleTableMap::iterator insert(ParticleTableMap::iterator hint, const ParticleTableEntry& entry);

    /// #ParticleFactory's own inserter
    friend std::insert_iterator<ParticleFactory> inserter(ParticleFactory& F)
    { return std::insert_iterator<ParticleFactory>(F, F.ParticleTable_.end()); }

    /// @}

private:

    /// maps PDGCodes to ParticleTableEntry's
    ParticleTableMap ParticleTable_;
};

}

#endif
