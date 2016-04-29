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

#ifndef yap_MassShapeWithNominalMass_h
#define yap_MassShapeWithNominalMass_h

#include "MassShape.h"

#include <memory>

namespace yap {

class ParticleCombination;

/// \class MassShapeWithNominalMass
/// \brief Class for MassShape that gets its nominal mass from its owning resonance
/// \author Daniel Greenwald
/// \ingroup MassShapes
class MassShapeWithNominalMass : public MassShape
{
public:

    /// Constructor
    MassShapeWithNominalMass() : MassShape() {}

    /// Set parameters from ParticleTableEntry
    /// \param entry ParticleTableEntry containing information to create mass shape object
    virtual void setParameters(const ParticleTableEntry& entry) override;

    /// Get mass
    std::shared_ptr<RealParameter> mass();

    /// Get mass (const)
    const std::shared_ptr<RealParameter> mass() const
    { return const_cast<MassShapeWithNominalMass*>(this)->mass(); }

    virtual std::string data_accessor_type() const override
    {return "MassShapeWithNominalMass"; }

protected:

    /// set dependency on resonance's mass
    virtual void setDependenciesFromResonance() override;

    /// set dependency on masses from model
    virtual void setDependenciesFromModel() override;

};

}

#endif
