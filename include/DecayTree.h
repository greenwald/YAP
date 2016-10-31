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

#ifndef yap_DecayTree_h
#define yap_DecayTree_h

#include "fwd/DecayTree.h"

#include "fwd/AmplitudeComponent.h"
#include "fwd/DataPoint.h"
#include "fwd/DecayChannel.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/Spin.h"
#include "fwd/SpinAmplitude.h"
#include "fwd/StatusManager.h"
#include "fwd/VariableStatus.h"

#include <array>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <vector>

namespace yap {

/// \class DecayTree
/// \brief Class holding vectors of fixed and free amplitudes that define a decay tree
/// \author Johannes Rauch, Daniel Greenwald
class DecayTree
{
public:

    /// \typedef DaughterDecayTreeMap
    /// key is daughter index, value is decay tree of daughter
    using DaughterDecayTreeMap = std::map<unsigned, std::shared_ptr<DecayTree> >;

    /// Constructor
    /// \param free_amp shared_ptr to ComplexParameter for the free amplitude
    /// \param two_M (twice) the spin projection of the parent particle
    /// \param two_m (twice) the spin projections of the daughter particles
    explicit DecayTree(std::shared_ptr<FreeAmplitude> free_amp, int two_M, const SpinProjectionVector& two_m);

    /// \return FreeAmplitude_
    const std::shared_ptr<FreeAmplitude>& freeAmplitude() const
    { return FreeAmplitude_; }

    /// \return InitialTwoM_
    const int initialTwoM() const
    { return InitialTwoM_; }

    /// \return FinalTwoM_
    const SpinProjectionVector& finalTwoM() const
    { return FinalTwoM_; }

    /// \return DaughterDecayTrees_
    const DaughterDecayTreeMap daughterDecayTrees() const
    { return DaughterDecayTrees_; }

    /// \return Model this DecayTree belongs to (via FreeAmplitude)
    const Model* model() const;

    /// \return product of all free amplitudes in this decay tree
    const std::complex<double> dataIndependentAmplitude() const;

    /// \return product of all particle-combination-dependent amplitudes in this tree
    /// summing over particle combinations of DecayChannel inside FreeAmplitude
    /// \param d DataPoint
    const std::complex<double> dataDependentAmplitude(const DataPoint& d) const;

    /// \return product of all particle-combination-dependent amplitudes in this tree
    /// \param d DataPoint
    /// \param pc shared_ptr<const ParticleCombination>
    const std::complex<double> dataDependentAmplitude(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const;

    /// \return VariableStatus of dataDependentAmplitude
    const VariableStatus dataDependentAmplitudeStatus() const;

    /// grant friend status to DecayChannel to call addDataAccessor
    friend class DecayChannel;

    /// grant friend status to DecayingState to call addDataAccessor
    friend class DecayingState;

protected:

    /// Set the DecayTree of the i'th daughter
    /// \param i index of daughter to set decay tree for
    /// \param dt shared_ptr to DecayTree to set
    void setDaughterDecayTree(unsigned i, std::shared_ptr<DecayTree> dt);

    /// Add an AmplitudeComponent
    void addAmplitudeComponent(const AmplitudeComponent& rda);

private:

    /// ComplexParameter of the free amplitude for the decay
    std::shared_ptr<FreeAmplitude> FreeAmplitude_;

    /// parent spin projection
    int InitialTwoM_;

    /// daughter spin projections
    SpinProjectionVector FinalTwoM_;

    /// vector of AmplitudeComponent's
    std::vector<const AmplitudeComponent*> AmplitudeComponents_;

    /// map of daughter index -> daughter DecayTree
    DaughterDecayTreeMap DaughterDecayTrees_;

};

/// convert to (multiline) string
std::string to_string(const DecayTree& dt, std::string offset = "");

/// convert to (mutliline string)
std::string to_string(const DecayTreeVector& dtv);

/// equality operator
inline bool operator==(const DecayTree& lhs, const DecayTree& rhs)
{ return lhs.freeAmplitude() == rhs.freeAmplitude() and lhs.daughterDecayTrees() == rhs.daughterDecayTrees(); }

/// \return Depth of DecayTree
unsigned depth(const DecayTree& DT);

/// \return amplitude evaluated for DataPoint over all symmetrizations
/// \param dt DecayTree to operate on
/// \param d DataPoint to evaluate on
inline const std::complex<double> amplitude(const DecayTree& dt, const DataPoint& d)
{ return dt.dataIndependentAmplitude() * dt.dataDependentAmplitude(d); }

/// \return amplitude evaluated for DataPoint for particular symmetrization
/// \param dt DecayTree to operate on
/// \param d DataPoint to evaluate on
/// \param pc ParticleCombination
inline const std::complex<double> amplitude(const DecayTree& dt, const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc)
{ return dt.dataIndependentAmplitude() * dt.dataDependentAmplitude(d, pc); }

/// \return sum of amplitudes of decay trees in a vector
/// \param dtv DecayTreeVector to sum over
/// \param d DataPoint to evaluate on
const std::complex<double> amplitude(const DecayTreeVector& dtv, const DataPoint& d);

/// \return square of sum of amplitudes of decay trees in a vector
/// \param dtv DecayTreeVector to sum over
/// \param d DataPoint to evaluate on
inline const double intensity(const DecayTreeVector& dtv, const DataPoint& d)
{ return norm(amplitude(dtv, d)); }

/// \return set of all free amplitudes in a DecayTree
FreeAmplitudeSet free_amplitudes(const DecayTree& DT);

/// \return set of all free amplitudes in a DecayTreeVector
FreeAmplitudeSet free_amplitudes(const DecayTreeVector& DTV);

/// \return whether a decay tree has changed
const bool has_changed(const std::shared_ptr<DecayTree>& dt);

/// \return vector of trees whose data-dependent amplitude variable statuses are VariableStatus::changed
/// \param vector of trees to check in
const DecayTreeVector select_changed(const DecayTreeVector& dtv);

}

#endif
