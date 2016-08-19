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

#ifndef yap_FourMomenta_
#define yap_FourMomenta_

#include "fwd/CachedValue.h"
#include "fwd/DataPoint.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/FourVector.h"
#include "fwd/MassAxes.h"
#include "fwd/Model.h"
#include "fwd/ParticleCombination.h"
#include "fwd/StatusManager.h"

#include "StaticDataAccessor.h"

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

namespace yap {

/// \class FourMomenta
/// \brief Stores and gives access to four-momenta and invariant masses
/// \author Johannes Rauch, Daniel Greenwald
class FourMomenta : public StaticDataAccessor
{
public:

    /// Constructor
    /// \param m Owning model
    FourMomenta(Model& m);

    /// check consistency
    bool consistent() const;

    /// Fill 4-momenta
    /// \param d DataPoint to fill
    /// \param sm StatusManager to update
    virtual void calculate(DataPoint& d, StatusManager& sm) const override;

    /// \name Getters
    /// @{

    /// Access 4-momenta (const)
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return 4-momentum of
    FourVector<double> p(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const;

    /// Access invariant mass squared
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return squared mass of
    double m2(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
    { return pow(m(d, pc), 2); }

    /// Access invariant mass
    /// \param d DataPoint to get data from
    /// \param pc ParticleCombination to return mass of
    double m(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const;

    /// \return total four-momentum (const)
    /// \param d DataPoint to get data from
    const FourVector<double> totalMomentum(const DataPoint& d) const;

    /// \return vector of final-state four-momenta (const)
    /// \param d DataPoint to get data from
    const std::vector<FourVector<double> > finalStateMomenta(const DataPoint& d) const;

    /// \return masses
    std::shared_ptr<RealCachedValue> mass()
    { return M_; }

    /// \return masses (const)
    std::shared_ptr<RealCachedValue> mass() const
    { return M_; }

    /// \return momentum
    std::shared_ptr<FourVectorCachedValue> momentum()
    { return P_; }

    /// \return momentum (const)
    std::shared_ptr<FourVectorCachedValue> momentum() const
    { return P_; }

    /// @}

    /// print all masses
    std::ostream& printMasses(const DataPoint& d, std::ostream& os = std::cout) const;

    /// grant friend status to Model to call addParticleCombination and setFinalStateMomenta
    friend class Model;

protected:

    /// set final-state four-momenta
    /// \param d DataPoint to set into
    /// \param P Four-momenta to set
    /// \param sm StatusManager to be updated
    void setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, StatusManager& sm) const;

    /// looks for ISP when adding ParticleCombination's
    void addParticleCombination(std::shared_ptr<ParticleCombination> pc) override;

    /// override to do nothing, since FourMomenta doesn't rely on parents being set.
    void pruneSymmetrizationIndices() override
    {}

private:

    /// Symmetrization index of an initial state that is composed from ALL final state particles
    int TotalIndex_;

    /// Symmetrization indices of final states
    std::vector<int> FSPIndices_;

    /// four-vector of particle combinations
    std::shared_ptr<FourVectorCachedValue> P_;

    /// invariant mass of particle combinations [GeV]
    std::shared_ptr<RealCachedValue> M_;

};

/// Calculate four-momenta for final-state particles for phase-space coordinate.
/// \param initial_mass initial mass of decaying system
/// \param FSPs Vector of final state particles
/// \param axes phase-space axes
/// \param squared_masses phase-space coordinate
std::vector<FourVector<double> > calculate_four_momenta(double initial_mass, const FinalStateParticleVector& FPSs,
                                                        const MassAxes& axes, const std::vector<double>& squared_masses);

/// Calculate four-momenta for final-state particles for phase-space coordinate
/// And apply rotation into model's coordinate system
/// \param initial_mass initial mass of decaying system
/// \param M model to get final state particles and coordinate system from
/// \param axes phase-space axes
/// \param squared_masses phase-space coordinate
std::vector<FourVector<double> > calculate_four_momenta(double initial_mass, const Model& M,
                                                        const MassAxes& axes, const std::vector<double>& squared_masses);

}
#endif

