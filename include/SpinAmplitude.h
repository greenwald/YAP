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

#ifndef yap_SpinAmplitude_h
#define yap_SpinAmplitude_h

#include "fwd/DataPointBase.h"

#include "MathUtilities.h"
#include "StaticDataAccessor.h"

#include <array>
#include <cstdlib>
#include <iostream>
#include <memory>

namespace yap {

class ComplexCachedDataValue;
class StatusManager;

/// \class SpinAmplitude
/// \brief Abstract base class implementing a spin amplitude.
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup SpinAmplitude Spin Amplitudes
class SpinAmplitude : public StaticDataAccessor
{
public:

    /// \typedef SpinProjectionPair
    using SpinProjectionPair = std::array<int, 2>;

    /// \typedef AmplitudeSubmap
    /// \brief maps SpinProjectionPair to ComplexCachesDataValue
    using AmplitudeSubmap = std::map<SpinProjectionPair, std::shared_ptr<ComplexCachedDataValue> >;

    /// \typedef AmplitudeMap
    /// \brief maps parent spin projectin to AmplitudeSubmap
    using AmplitudeMap = std::map<int, AmplitudeSubmap>;

    /// \return whether three spins fulfill the triangle relationship
    /// \param two_a 2 * spin a
    /// \param two_b 2 * spin b
    /// \param two_c 2 * spin c
    /// \return \f$ \Delta(abc) \f$
    static constexpr bool triangle(unsigned two_a, unsigned two_b, unsigned two_c)
    { return is_even(two_a + two_b + two_c) and std::abs<int>(two_a - two_b) <= two_c and two_c <= (two_a + two_b); }

    /// \return Whether angular momentum is conserved in J -> j1 + j2 with orbital angular momentum l
    /// \param two_J 2 * spin of initial state
    /// \param two_j1 2 * spin of first daughter
    /// \param two_j2 2 * spin of second daughter
    /// \param l orbital angular momentum
    static constexpr bool conserves(unsigned two_J, unsigned two_j1, unsigned two_j2, int l)
    {
        // check that the spins are consistent, and that the triangle requirements for (Jlj) and (j1j2j) can be met simultaneously
        return is_even(two_J + two_j1 + two_j2)
               and (std::min(two_j1 + two_j2, two_J + 2 * l) >= std::max(std::abs<int>(two_j1 - two_j2), std::abs<int>(two_J - 2 * l)));
    }

    /// cast into string
    virtual operator std::string() const;

    /// check consistency of object
    virtual bool consistent() const
    { return true; }

    /// \name Getters
    /// @{

    /// \return initial spin * 2 (const)
    unsigned initialTwoJ() const
    { return InitialTwoJ_; }

    /// \return array of spin * 2 of daughters (const)
    const std::array<unsigned, 2>& finalTwoJ() const
    { return FinalTwoJ_; }

    /// \return orbital angular momentum
    unsigned L() const
    { return L_; }

    /// \return total spin angular momentum */
    unsigned twoS() const
    { return TwoS_; }

    /// \return set of (twice the) spin projections of initial state
    std::set<int> twoM() const;

    /// @}

    /// Calculate spin amplitude for caching.
    /// Must be overrided in derived classes.
    /// \param two_M 2 * spin projection of parent
    /// \param two_m1 2 * spin projection of first daughter
    /// \param two_m2 2 * spin projection of second daughter
    /// \param d DataPointBase to retrieve data from for calculation
    /// \param pc ParticleCombination to calculate for
    virtual std::complex<double> calc(int two_M, int two_m1, int two_m2,
                                      const DataPointBase& d, const std::shared_ptr<ParticleCombination>& pc) const = 0;

    /// Loops over particle combinations (pc) and all (M, m1, m2) combinations
    /// and call calc(M, m1, m2, d, pc) when necessary
    /// \param d DataPointBase to calculate into
    /// \param sm StatusManager to update
    void calculate(DataPointBase& d, StatusManager& sm) const override;

    /// \return precalculated complex amplitude
    /// \param d DataPointBase to retrieve value from
    /// \param pc ParticleCombination to retrieve value for
    /// \param two_M 2 * spin projection of parent
    /// \param two_m1 2 * spin projection of first daughter
    /// \param two_m2 2 * spin projection of second daughter
    std::complex<double> amplitude(const DataPointBase& d, const std::shared_ptr<ParticleCombination>& pc,
                                   int two_M, int two_m1, int two_m2) const
    { return amplitude(two_M, two_m1, two_m2)->value(d, symmetrizationIndex(pc)); }

    /// access cached spin amplitude
    /// \param two_M 2 * spin projection of parent
    /// \param two_m1 2 * spin projection of first daughter
    /// \param two_m2 2 * spin projection of second daughter
    std::shared_ptr<ComplexCachedDataValue>& amplitude(int two_M, int two_m1, int two_m2)
    { return Amplitudes_.at(two_M).at({two_m1, two_m2}); }

    /// access cached spin amplitude (const)
    /// \param two_M 2 * spin projection of parent
    /// \param two_m1 2 * spin projection of first daughter
    /// \param two_m2 2 * spin projection of second daughter
    const std::shared_ptr<ComplexCachedDataValue>& amplitude(int two_M, int two_m1, int two_m2) const
    { return Amplitudes_.at(two_M).at({two_m1, two_m2}); }

    /// \return set of cached spin amplitudes
    CachedDataValueSet amplitudeSet();

    /// \return AmplitudeMap Amplitudes_
    const AmplitudeMap& amplitudes() const
    { return Amplitudes_; }

    virtual std::string data_accessor_type() const override
    {return "SpinAmplitude"; }

    /// \return a string naming the formalism used for the SpinAmplitude calculation
    virtual std::string formalism() const = 0;

    /// grant friend access to SpinAmplitudeCache to set Model
    friend class SpinAmplitudeCache;

    /// grant friend access to DecayChannel to call addParticleCombination
    friend class DecayChannel;

protected:

    /// add spin amplitude for transition from state with initial
    /// projection to states with final projections.
    /// \param two_M twice the spin projection of the initial state
    /// \param two_m1 twice the spin projection of the first daughter
    /// \param two_m2 twice the spin projection of the second daughter
    virtual void addAmplitude(int two_M, int two_m1, int two_m2);

    /// check equivalence: only check spins and angular momenta
    virtual bool equiv(const SpinAmplitude& other) const;

    /// check equality: calls #equiv() and checks symmetrizationIndices
    virtual bool equals(const SpinAmplitude& other) const;

    /// Constructor
    /// declared private to ensure SpinAmplitude's are only created by a SpinAmplitudeCache
    /// \param two_J  twice the spin of Initial-state
    /// \param two_j1 twice the spin of first daughter
    /// \param two_j2 twice the spin of second daughter
    /// \param l orbital angular momentum
    /// \param two_s twice the total spin angular momentum
    /// \param equiv ParticleCombination equivalence struct for determining index assignments
    SpinAmplitude(unsigned two_J, unsigned two_j1, unsigned two_j2, unsigned l, unsigned two_s,
                  ParticleCombination::Equiv* equiv = &ParticleCombination::equivBySharedPointer);

    /// set raw pointer to owning model; overrides to call setDependencies
    virtual void setModel(Model* m) override;

private:

    /// set dependencies for amplitudes
    virtual void setDependencies(std::shared_ptr<CachedDataValue> a) = 0;

    /// Initial-state spin * 2
    unsigned InitialTwoJ_;

    /// array of final-state spin * 2
    std::array<unsigned, 2> FinalTwoJ_;

    /// orbital angular momentum
    unsigned L_;

    /// twice the total spin angular momentum
    unsigned TwoS_;

    /// Cached complex spin amplitude
    AmplitudeMap Amplitudes_;

    /// equality operator
    friend bool operator==(const SpinAmplitude& A, const SpinAmplitude& B)
    { return typeid(A) == typeid(B) and A.equals(B); }

};

/// convert to string
inline std::string to_string(const SpinAmplitude& sa)
{ return (std::string)sa; }

/// << operator
inline std::ostream& operator<< (std::ostream& os, const SpinAmplitude& sa)
{ os << to_string(sa); return os; }

/// \typedef SpinAmplitudeVector
using SpinAmplitudeVector = std::vector<std::shared_ptr<SpinAmplitude> >;

/// convert to string
std::string to_string(const SpinAmplitudeVector& saV);

/// \typedef SpinAmplitudeMap
/// \tparam T Object to store in map, with shared_ptr to SpinAmplitude as key
template<typename T>
using SpinAmplitudeMap = std::map<std::shared_ptr<SpinAmplitude>, T,
      std::owner_less<std::shared_ptr<SpinAmplitude> > >;

}

#endif
