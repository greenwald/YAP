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

#include "fwd/Wave.h"

#include "fwd/Particle.h"
#include "fwd/QuantumNumbers.h"

#include "DecayingState.h"

#include <memory>

namespace yap {

/// Class for a complete description of the full dynamic range of
/// decay to a particular final state
class Wave : public DecayingState
{
protected:

    /// Constructor
    /// see #create
    Wave(const std::string& name, const QuantumNumbers& q, unsigned l, unsigned two_s, double radial_size, const ParticleVector& daughters);

public:

    /// create
    /// \param name Name of wave
    /// \param q QuantumNumbers of wave
    /// \param l orbital angular momentum of wave
    /// \param s spin angular momentum of wave
    /// \param radial_size radial size of wave
    /// \param daughters daughters of wave
    static std::shared_ptr<Wave> create(const std::string& name, const QuantumNumbers& q, unsigned l, unsigned two_s, double radial_size, const ParticleVector& daughters)
    { return std::shared_ptr<Wave>(new Wave(name, q, l, two_s, radial_size, daughters)); }

    /// create, l & s provided explicitly
    template <typename ... Types>
    std::shared_ptr<Wave> create(const std::string& name, const QuantumNumbers& q, unsigned l, unsigned two_s, double radial_size, std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, Types ... other_daughters)
    { return create(name, q, l, two_s, radial_size, ParticleVector({A, B, other_daughters...})); }

    /// create
    /// spin angular momentum is calculated; if ambiguous an exception is thrown
    /// \param name Name of wave
    /// \param q QuantumNumbers of wave
    /// \param l orbital angular momentum of wave
    /// \param radial_size radial size of wave
    /// \param daughters daughters of wave
    static std::shared_ptr<Wave> create(const std::string& name, const QuantumNumbers& q, unsigned l, double radial_size, const ParticleVector& daughters);

    /// create
    template <typename ... Types>
    std::shared_ptr<Wave> create(const std::string& name, const QuantumNumbers& q, unsigned l, double radial_size, std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, Types ... other_daughters)
    { return create(name, q, l, radial_size, ParticleVector({A, B, other_daughters...})); }

    /// create
    /// spin and orbital angular momenta are calculated; if ambiguous an exception is thrown
    /// \param name Name of wave
    /// \param q QuantumNumbers of wave
    /// \param l orbital angular momentum of wave
    /// \param radial_size radial size of wave
    /// \param daughters daughters of wave
    static std::shared_ptr<Wave> create(const std::string& name, const QuantumNumbers& q, double radial_size, const ParticleVector& daughters);

    /// create
    template <typename ... Types>
    std::shared_ptr<Wave> create(const std::string& name, const QuantumNumbers& q, double radial_size, std::shared_ptr<Particle> A, std::shared_ptr<Particle> B, Types ... other_daughters)
    { return create(name, q, radial_size, ParticleVector({A, B, other_daughters...})); }

private:
    
};

/// checks if something inherits from Wave
extern const is_of_type<Wave> is_wave;
 
}

#endif
