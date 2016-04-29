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

#ifndef yap_BreitWigner_h
#define yap_BreitWigner_h

#include "fwd/DataPointBase.h"

#include "MassShapeWithNominalMass.h"

#include <complex>
#include <memory>

namespace yap {

class DataPartition;
class ParticleCombination;

/// \class BreitWigner
/// \brief Class for Breit-Wigner resonance shape
/// \author Daniel Greenwald
/// \ingroup MassShapes
///
/// Amplitude is 1 / (mass^2 - s - i*mass*width)\n\n

class BreitWigner : public MassShapeWithNominalMass
{
public:

    /// Constructor
    /// \param width Width of resonance [GeV]
    BreitWigner(double w = -1);

    /// Calculate complex amplitude
    /// \param d DataPointBase to calculate with
    /// \param pc (shared_ptr to) ParticleCombination to calculate for
    /// \param sm StatusManager to update
    virtual std::complex<double> amplitude(DataPointBase& d, const std::shared_ptr<ParticleCombination>& pc, StatusManager& sm) const override;

    /// Calculate complex amplitudes for and store in each DataPointBase in DataPartition
    /// Must be overrided in derived classes.
    /// \param pc (shared_ptr to) ParticleCombination to calculate for
    virtual void calculate(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc) const override;

    /// Set parameters from ParticleTableEntry
    /// \param entry ParticleTableEntry containing information to create mass shape object
    virtual void setParameters(const ParticleTableEntry& entry) override;

    /// Get width
    std::shared_ptr<RealParameter> width()
    { return Width_; }

    /// Get width (const)
    const std::shared_ptr<RealParameter> width() const
    { return const_cast<BreitWigner*>(this)->width(); }

    virtual bool consistent() const override;

    virtual std::string data_accessor_type() const override
    {return "BreitWigner"; }

private:

    std::shared_ptr<RealParameter> Width_; ///< [GeV]

};

}

#endif
