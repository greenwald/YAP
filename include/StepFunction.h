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

#ifndef yap_StepFunction_h
#define yap_StepFunction_h

#include "fwd/StepFunction.h"

#include "fwd/DataPartition.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/Parameter.h"
#include "fwd/ParticleCombination.h"
#include "fwd/StatusManager.h"

#include "MassShape.h"
#include "VariableStatus.h"

#include <complex>
#include <memory>
#include <vector>

namespace yap {

/// \class StepFunction
/// \brief A step-function mass shape
/// \author Daniel Greenwald, Paolo Di Giglio
/// \ingroup MassShapes
class StepFunction : public MassShape
{
public:

    /// \typedef LowEdgeVector
    using LowEdgeVector = std::vector<std::shared_ptr<NonnegativeRealParameter> >;

    /// individual step in a StepFunction
    class Step : public MassShape
    {
    public:
        /// constructor
        /// \param SF owner
        /// \param low_edge Low edge of step [GeV^2]
        /// \param high_edge High edge of step [GeV^2]
        Step(const StepFunction& SF, const LowEdgeVector::value_type& low_edge, const LowEdgeVector::value_type& high_edge);

        /// get raw pointer to owner throw StepFunction
        virtual DecayingParticle* owner() const override
        { return StepFunction_->owner(); }

        /// \return LowEdge_
        const LowEdgeVector::value_type& lowEdge() const
        { return LowEdge_; }

        /// \return HighEdge_
        const LowEdgeVector::value_type& highEdge() const
        { return HighEdge_; }
        
        using MassShape::calculate;

        /// Does nothing
        virtual void calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const override {}

        /// \return 1 if mass lies in [LowEdge, HighEdge), 0 otherwise
        /// \param d DataPoint
        /// \param pc shared_ptr to ParticleCombination
        virtual const std::complex<double> value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const override;

        /// does nothing
        virtual void updateCalculationStatus(StatusManager& D) const override {}

        /// Grant friend status to StepFunction
        friend StepFunction;
        
    private:
        /// Raw pointer to StepFunction to which Step belongs
        const StepFunction* StepFunction_;
        
        /// Low edge of step [GeV^2]
        LowEdgeVector::value_type LowEdge_;

        /// High edge of step [GeV^2]
        LowEdgeVector::value_type HighEdge_;
    };

    /// \typedef StepVector
    using StepVector = std::vector<std::shared_ptr<Step> >;
    
    /// Constructor;
    /// last low-edge is upper edge of last interval
    /// \param low_edges low edges of step-function intervals [GeV^2]
    StepFunction(const std::vector<double>& low_edges);

    /// Constructor;
    /// last low-edge is upper edge of last interval
    /// \param low_edges low edges of step-function intervals [GeV^2]
    template <typename ... Others>
    StepFunction(double l0, double l1, Others ... L) : StepFunction(std::vector<double>({l0, l1, L...})) {}

    /// \return LowEdges_ (const)
    const LowEdgeVector& lowEdges() const
    { return LowEdges_; }

    /// \return FreeAmplitudes_
    FreeAmplitudeVector& freeAmplitudes()
    { return FreeAmplitudes_; }

    /// \return FreeAmplitudes_ (const)
    const FreeAmplitudeVector& freeAmplitudes() const
    { return FreeAmplitudes_; }

    /// \return Steps_
    const StepVector& steps() const
    { return Steps_; }
    
    /// set a low edge (checks for monotonicity)
    const VariableStatus setLowEdge(unsigned index, double low_edge);

    /// does nothing
    virtual void updateCalculationStatus(StatusManager& D) const override {}
    
    /// \return 0; should never be called
    /// \param d DataPoint
    /// \param pc shared_ptr to ParticleCombination
    virtual const std::complex<double> value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const override
    { return 0; }

    using MassShape::calculate;
    
    /// does nothing
    virtual void calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const override {}

protected:

    // create copies of DecayTree's in owner of DecayChannel for each FreeAmplitude necessary
    virtual void addDecayChannel(std::shared_ptr<DecayChannel> c) override;

    /// call MassShape::addParticleCombination and add ParticleCombination to Steps
    virtual void addParticleCombination(const ParticleCombination& pc) override;

private:
    
    /// low edges of intervals [GeV^2]
    LowEdgeVector LowEdges_;

    /// steps
    StepVector Steps_;
    
    /// FreeAmplitudes for each interval
    FreeAmplitudeVector FreeAmplitudes_;
};

}

#endif
