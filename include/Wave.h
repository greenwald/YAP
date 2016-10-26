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

#ifndef yap_Wave_h
#define yap_Wave_h

#include "DecayingState.h"

/// Class for a unified description of decay to a particular state
/// across the full dynamic range
/// \author Daniel Greenwald
/// \ingroup Particle
class Wave : public DecayingState
{
protected:
    
    /// Constructor
    /// see #create
    Wave(const std::string& name, const QuantumNumbers& q, unsigned l, unsigned two_s,
         const ParticleVector& daughters, double radial_size);

public:
    
    /// create
    /// \param name Name of wave
    /// \param q QuantumNumbers of wave
    /// \param l orbital angular momentum of wave
    /// \param two_s twice the spin angular momentum of the wave
    /// \param daughters ParticleVector of daughters of wave
    /// \param radialSize radial size of wave
    static std::shared_ptr<Wave> create(const std::string& name, const QuantumNumbers& q, unsigned l, unsigned two_s,
                                        const ParticleVector& daughters, double radialSize)
    { return std::shared_ptr<Wave>(new Wave(name, q, l, two_s, daughters, radialSize)); }

    
    
private:

    std::shared_ptr<DecayChannel> DecayChannel_;

    std::shared_ptr<BlattWeisskopf> BlattWeisskopf_;
    
};

#endif

