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

#ifndef yap_MassRange_h
#define yap_MassRange_h

#include "fwd/MassRange.h"

#include "fwd/FinalStateParticle.h"
#include "fwd/MassAxes.h"
#include "fwd/ParticleCombination.h"

#include <memory>

namespace yap {

/// \return mass range for particle combination inside decay of ISP to FSPs
/// \param isp_mass initial-state mass
/// \param pc shared_ptr to ParticleCombination to get mass range for
/// \param FSP FinalStateParticles of the ISP to calculate with
const MassRange mass_range(double isp_mass, const std::shared_ptr<ParticleCombination>& pc, const FinalStateParticleVector& FSPs);

/// \return mass range for MassAxes inside decay of ISP to FSPs
/// \param isp_mass initial-state mass
/// \param A MassAxes to get mass ranges for
/// \param FSPs FinalStateParticles of the ISP to calculate with
const std::vector<MassRange> mass_range(double isp_mass, const MassAxes& A, const FinalStateParticleVector& FSPs);

/// \return squared masses of mass range
const MassRange squared(MassRange mr);

/// \return squared masses of vector of mass ranges
const std::vector<MassRange> squared(std::vector<MassRange> mr);

}

#endif
