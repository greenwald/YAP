#include "StepFunction.h"

#include "Attributes.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "DecayChannel.h"
#include "DecayTree.h"
#include "Exceptions.h"
#include "Filter.h"
#include "FourMomenta.h"
#include "FreeAmplitude.h"
#include "Group.h"
#include "Model.h"
#include "Parameter.h"
#include "StatusManager.h"

#include <algorithm>

#include "logging.h"

namespace yap {

//-------------------------
StepFunction::Step::Step(const StepFunction& SF, const StepFunction::LowEdgeVector::value_type& low_edge, const StepFunction::LowEdgeVector::value_type& high_edge) :
    StepFunction_(&SF),
    LowEdge_(low_edge),
    HighEdge_(high_edge)
{
    if (HighEdge_->value() <= LowEdge_->value())
        throw exceptions::Exception("high edge <= low edge", "StepFunction::Step::Step");
    addParameter(LowEdge_);
    addParameter(HighEdge_);
}

//-------------------------
const std::complex<double> StepFunction::Step::value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
{
    auto squared_mass = model()->fourMomenta()->m2(d, pc);
    if (squared_mass > LowEdge_->value() and squared_mass <= HighEdge_->value())
        return 1;
    return 0;
}

//-------------------------
StepFunction::StepFunction(const std::vector<double>& low_edges)
{
    // check side
    if (low_edges.size() < 2)
        throw exceptions::Exception("At least two low edges must be specified", "StepFunction::StepFunction");
    
    // create parameters and steps (which check for sorting of edges)
    for (size_t i = 0; i < low_edges.size(); ++i) {
        LowEdges_.push_back(std::make_shared<LowEdgeVector::value_type::element_type>(low_edges[i]));
        addParameter(LowEdges_.back());
        if (i > 0)
            Steps_.push_back(std::make_shared<Step>(*this, LowEdges_[i - 1], LowEdges_[i]));
    }
}

//-------------------------
const VariableStatus StepFunction::setLowEdge(unsigned index, double low_edge)
{
    if (index >= LowEdges_.size())
        throw exceptions::Exception("index is beyond LowEdges_ size", "StepFunction::setLowEdge");

    if (index > 0 and low_edge <= LowEdges_[index - 1]->value())
        throw exceptions::Exception("value is less than or equal to previous interval edge", "StepFunction::setLowEdge");

    if (index < LowEdges_.size() - 1 and low_edge >= LowEdges_[index + 1]->value())
        throw exceptions::Exception("value is greater than or equal to next interval edge", "StepFunction::setLowEdge");

    return LowEdges_[index]->setValue(low_edge);
}

//-------------------------
void StepFunction::addDecayChannel(std::shared_ptr<DecayChannel> c)
{
    if (!FreeAmplitudes_.empty())
        throw exceptions::Exception("DecayChannel has already been added.", "StepFunction::addDecayChannel");

    // get owners DecayTrees with DecayChannel c and MassShape this
    auto DTV = filter(ownersDecayTrees(), has_decay_channel(c), has_mass_shape(this));

    // and group by FreeAmplitude
    auto DTV_by_FA = group(DTV, [](const DecayTreeVector::value_type& A, const DecayTreeVector::value_type& B)
                           {return A->freeAmplitude() < B->freeAmplitude();});

    // loop over groups
    for (auto& dtv : DTV_by_FA) {
        // all DecayTree's in dtv have the same FreeAmplitude
        // for each Step, create copies of all DecayTree's replacing MassShape with Step and with a new copy of the FreeAmplitude
        for (const auto& step : Steps_) {
            auto fa = std::make_shared<FreeAmplitude>(*dtv[0]->freeAmplitude());
            for (const auto& dt : dtv) {
                // copy dt
                auto copy_dt = std::make_shared<DecayTree>(*dt);
                // remove this from AmplitudeComponents
                removeAmplitudeComponent(*this, *copy_dt);
                // add Step to AmpltiudeComponents
                addAmplitudeComponent(*step, *copy_dt);
                // replace FreeAmplitude
                replaceFreeAmplitude(fa, *copy_dt);
                // add new DecayTree
                ownersDecayTrees().push_back(copy_dt);
            }
            FreeAmplitudes_.push_back(fa);
        }
    }

    // new remove the pre-existing DecayTrees
    for (const auto& DT : DTV) {
        auto it = std::find(ownersDecayTrees().begin(), ownersDecayTrees().end(), DT);
        if (it != ownersDecayTrees().end())
            ownersDecayTrees().erase(it);
    }
}

//-------------------------
void StepFunction::addParticleCombination(const ParticleCombination& pc)
{
    MassShape::addParticleCombination(pc);
    for (auto& step : Steps_)
        step->addParticleCombination(pc);
}

}
