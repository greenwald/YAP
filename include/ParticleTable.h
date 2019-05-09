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

#ifndef yap_ParticleTable_h
#define yap_ParticleTable_h

#include "fwd/ParticleTable.h"

#include "fwd/DecayingParticle.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/MassShape.h"

#include "QuantumNumbers.h"

#include <memory>
#include <string>
#include <vector>

namespace yap {

/// \struct ParticleTableEntry
/// \brief Data container for storing particle information in database
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup ParticleTable
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

    /// \return QuantumNumbers_ (non const)
    QuantumNumbers& quantumNumbers()
    { return QuantumNumbers_; }

    /// \return Mass_
    const double mass() const
    { return Mass_; }

    /// \return MassShapeParameters_
    const std::vector<double>& massShapeParameters() const
    { return MassShapeParameters_; }

    /// set pdg
    void pdg(int p)
    { PDG_ = p; }
    
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

/// \return cp conjugate copy of a particle table entry
ParticleTableEntry cp_conjugate(ParticleTableEntry t, const std::string& prefix = "anti-");

/// \return n'th element in ParticleTableEntry's MassShapeParameters_, throwing error if necessary
/// \param pde ParticleTableEntry to access
/// \param n element to retrieve
/// \param where function name to point to when throwing
double get_nth_element(const ParticleTableEntry& pde, size_t n, const std::string& where);

/// \class ParticleTable
/// \brief Factory class for easy creation of Particle objects from PDG codes.
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Particle
/// \defgroup ParticleTable
class ParticleTable
{
public:

    /// \typedef ParticleTable::value_type
    /// Define this to allow `std::inserter` to use `insert`
    using value_type = ParticleTableEntry;

    /// \typedef ParticleTable::iterator
    /// Define this to allow `std::inserter` to use `insert`
    using iterator = ParticleTableMap::iterator;

    /// \typedef ParticleTable::const_iterator
    using const_iterator = ParticleTableMap::const_iterator;

    /// Adds content of rhs to this
    /// \param rhs ParticleTable to add into this
    ParticleTable& operator+=(const ParticleTable& rhs);

    /// get ParticleTableEntry from #ParticleTableMap_ with safety checks
    /// \param PDG pdg code labeling particle table entry
    const ParticleTableEntry& operator[](int PDG) const;

    /// get ParticleTableEntry from #ParticleTableMap_ with safety checks
    /// \param PDG pdg code labeling particle table entry
    ParticleTableEntry& operator[](int PDG)
    { return const_cast<ParticleTableEntry&>(static_cast<const ParticleTable*>(this)->operator[](PDG)); }

    /// get ParticleTableEntry from #ParticleTableMap_ with safety checks
    /// \param name Name of particle in table
    const ParticleTableEntry& operator[](const std::string& name) const;

    /// get ParticleTableEntry from #ParticleTableMap_ with safety checks
    /// \param name Name of particle in table
    ParticleTableEntry& operator[](const std::string& name)
    { return const_cast<ParticleTableEntry&>(static_cast<const ParticleTable*>(this)->operator[](name)); }

    /// inserts the pair `ParticleTableEntry::PDG` and `ParticleTableEntry` to #ParticleTableMap_
    /// \param entry a A ParticleTableEntry to add to #ParticleTableMap_
    std::pair<iterator, bool> insert(const ParticleTableEntry& entry);

    /// convenience function to allow `inserter()` to be used in the `std::copy` algorithm
    iterator insert(iterator hint, const ParticleTableEntry& entry);

    /// #ParticleTable's own inserter
    friend std::insert_iterator<ParticleTable> inserter(ParticleTable& F)
    { return std::insert_iterator<ParticleTable>(F, F.ParticleTableMap_.end()); }

    /// \return ParticleTableMap_.begin()
    iterator begin()
    { return ParticleTableMap_.begin(); }

    /// \return ParticleTableMap_.begin()
    const_iterator begin() const
    { return ParticleTableMap_.begin(); }

    /// \return ParticleTableMap_.end()
    iterator end()
    { return ParticleTableMap_.end(); }

    /// \return ParticleTableMap_.end()
    const_iterator end() const
    { return ParticleTableMap_.end(); }

    // \return ParticleTableMap_.empty()
    const bool empty() const
    { return ParticleTableMap_.empty(); }
    
    // \return ParticleTableMap_.size()
    const ParticleTableMap::size_type size() const
    { return ParticleTableMap_.size(); }

private:

    /// maps PDGCodes to ParticleTableEntry's
    ParticleTableMap ParticleTableMap_;
};

/// \return cp conjugate copy of a particle table entry
ParticleTableEntry cp_conjugate(const ParticleTableMap::value_type& vt, const std::string& prefix = "anti-");

}

#endif
