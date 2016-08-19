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

#ifndef yap_PoleMass_h
#define yap_PoleMass_h

#include "fwd/DataPoint.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"

#include "MassShape.h"

#include <complex>
#include <memory>

namespace yap {

/// \class PoleMass
/// \brief Class for pole-mass resonance shape
/// \author Daniel Greenwald
/// \ingroup MassShapes
///
/// Amplitude is 1 / (mass^2 - s)\n\n
/// mass is complex

class PoleMass : public MassShape
{
public:

    /// Constructor
    /// \param mass Mass of resonance [GeV]
    PoleMass(std::complex<double> mass = std::complex<double>(-1, -1));

    /// Set parameters from ParticleTableEntry;
    /// If width is available, sets M = mass + i/2 * width
    /// \param entry ParticleTableEntry containing information to create mass shape object
    virtual void setParameters(const ParticleTableEntry& entry) override;

    /// Get mass
    std::shared_ptr<ComplexParameter> mass() const
    { return Mass_; }

    /// Check consistency of object
    virtual bool consistent() const override;

protected:

    /// Calculate dynamic amplitude T for and store in each DataPoint in DataPartition
    /// \param D DataPartition to calculate on
    /// \param pc ParticleCombination to calculate for
    /// \param si SymmetrizationIndec to calculate for
    virtual void calculateT(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc, unsigned si) const override;

    /// Complex mass [GeV]
    std::shared_ptr<ComplexParameter> Mass_;

};

}

#endif
