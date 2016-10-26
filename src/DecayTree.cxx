#include "DecayTree.h"

#include "AmplitudeComponent.h"
#include "Attributes.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FreeAmplitude.h"
#include "Model.h"
#include "Particle.h"
#include "ParticleCombination.h"
#include "SpinAmplitude.h"
#include "VariableStatus.h"

#include <iterator>
#include <functional>

namespace yap {

//-------------------------
DecayTree::DecayTree(std::shared_ptr<FreeAmplitude> free_amp, int two_M, const SpinProjectionVector& two_m) :
    FreeAmplitude_(free_amp),
    InitialTwoM_(two_M),
    FinalTwoM_(two_m)
{
    if (!FreeAmplitude_)
        throw exceptions::Exception("FreeAmplitude is nullptr", "DecayTree::DecayTree");

    if (!FreeAmplitude_->spinAmplitude())
        throw exceptions::Exception("FreeAmplitude's SpinAmplitude is nullptr", "DecayTree::DecayTree");

    if (!FreeAmplitude_->decayChannel())
        throw exceptions::Exception("FreeAmplitude's DecayChannel is nullptr", "DecayTree::DecayTree");
}

//-------------------------
const Model* DecayTree::model() const
{
    return (FreeAmplitude_) ? FreeAmplitude_->model() : nullptr;
}

//-------------------------
FreeAmplitudeSet free_amplitudes(const DecayTree& DT)
{
    FreeAmplitudeSet S = {DT.freeAmplitude()};
    for (auto& d_dt : DT.daughterDecayTrees()) {
        auto s = free_amplitudes(*d_dt.second);
        S.insert(s.begin(), s.end());
    }
    return S;
}

//-------------------------
FreeAmplitudeSet free_amplitudes(const DecayTreeVector& DTV)
{
    FreeAmplitudeSet S;
    for (auto& DT : DTV) {
        auto s = free_amplitudes(*DT);
        S.insert(s.begin(), s.end());
    }
    return S;
}

//-------------------------
const std::complex<double> DecayTree::dataDependentAmplitude(const DataPoint& d) const
{
    std::complex<double> A = 0;
    for (const auto& pc : FreeAmplitude_->decayChannel()->particleCombinations()) {
        A += dataDependentAmplitude(d, pc);
    }
    return A;
}

//-------------------------
const std::complex<double> DecayTree::dataDependentAmplitude(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
{
    // spin amplitude
    auto A = FreeAmplitude_->spinAmplitude()->amplitude(d, pc, InitialTwoM_, FinalTwoM_);

    // amplitude components
    for (const auto& ac : AmplitudeComponents_)
        A *= ac->value(d, pc);

    // likewise for daughters
    for (const auto& d_dt : DaughterDecayTrees_)
        A *= d_dt.second->dataDependentAmplitude(d, pc->daughters()[d_dt.first]);

    return A;
}

//-------------------------
const VariableStatus DecayTree::dataDependentAmplitudeStatus() const
{
    static variable_status S;

    // check components
    auto s = S(AmplitudeComponents_.begin(), AmplitudeComponents_.end());
    if (s == VariableStatus::changed)
        return s;
    
    // check daughters
    for (const auto& d_dt : DaughterDecayTrees_)
        if ( (s *= d_dt.second->dataDependentAmplitudeStatus()) == VariableStatus::changed)
            return s;
    
    return s;
}

//-------------------------
const std::complex<double> DecayTree::dataIndependentAmplitude() const
{
    // spin amplitude
    auto A = FreeAmplitude_->value();
    // likewise for daughters
    for (const auto& d_dt : DaughterDecayTrees_)
        A *= d_dt.second->dataIndependentAmplitude();
    return A;
}

//-------------------------
void DecayTree::setDaughterDecayTree(unsigned i, std::shared_ptr<DecayTree> dt)
{
    if (i >= FinalTwoM_.size())
        throw exceptions::Exception("index exceeds number of daughters", "DecayTree::setDaughterDecayTree");

    if (!dt->freeAmplitude())
        throw exceptions::Exception("DecayTree's free amplitude is nullptr", "DecayTree:setDaughterDecayTree");

    if (DaughterDecayTrees_.find(i) != DaughterDecayTrees_.end())
        throw exceptions::Exception("DecayTree for this daughter already set", "DecayTree::setDaughterDecayTree");

    DaughterDecayTrees_[i] = dt;
}

//-------------------------
void DecayTree::addAmplitudeComponent(const AmplitudeComponent& ac)
{
    if (!FreeAmplitude_)
        throw exceptions::Exception("FreeAmplitude is nullptr", "DecayTree::addAmplitudeComponent");

    for (const auto& pc : FreeAmplitude_->particleCombinations())
        if (!ac.validFor(*pc))
            throw exceptions::Exception("AmplitudeComponent not valid for all ParticleCombinations required by FreeAmplitude",
                                        "DecayTree::addAmplitudeComponent");

    AmplitudeComponents_.push_back(&ac);
}

//-------------------------
std::shared_ptr<DecayingParticle> decayingParticle(const DecayTree& dt)
{
    if (!dt.model())
        throw exceptions::Exception("model is nullptr", "decayingParticle(DecayTree)");

    return std::static_pointer_cast<DecayingParticle>(particle(*dt.model(), has_decay_tree(dt)));
}

//-------------------------
std::string to_string(const DecayTree& dt, std::string offset)
{
    auto s = to_string(*dt.freeAmplitude()->decayChannel(), dt.initialTwoM(), dt.finalTwoM())
        + ", L = " + std::to_string(dt.freeAmplitude()->spinAmplitude()->L())
        + ", S = " + spin_to_string(dt.freeAmplitude()->spinAmplitude()->twoS());
    offset += "     ";
    for (const auto& d_dt : dt.daughterDecayTrees())
        s += "\n" + offset + to_string(*d_dt.second, offset);
    return s;
}

//-------------------------
std::string to_string(const DecayTreeVector& dtv)
{
    return std::accumulate(dtv.begin(), dtv.end(), std::string(),
                           [](std::string& s, const DecayTreeVector::value_type& dt)
                           { return s += "\n" + to_string(*dt); }).erase(0, 1);
}

// //-------------------------
// std::string to_string(const DecayTreeVectorMap& m_dtv_map)
// {
//     return std::accumulate(m_dtv_map.begin(), m_dtv_map.end(), std::string(),
//                            [](std::string & s, const DecayTreeVectorMap::value_type & m_dtv)
//                            { return s += to_string(m_dtv.second); });
// }

//-------------------------
std::shared_ptr<FreeAmplitude> free_amplitude(const DecayTree& dt)
{
    return dt.freeAmplitude();
}

//-------------------------
unsigned depth(const DecayTree& DT)
{
    unsigned d = 0;
    for (const auto& i_dt : DT.daughterDecayTrees())
        d = std::max(d, depth(*i_dt.second));
    return d + 1;
}

//-------------------------
const bool has_changed(const std::shared_ptr<DecayTree>& dt)
{
    return dt->dataDependentAmplitudeStatus() == VariableStatus::changed;
}

//-------------------------
const DecayTreeVector select_changed(const DecayTreeVector& dtv)
{
    DecayTreeVector C;
    C.reserve(dtv.size());
    std::copy_if(dtv.begin(), dtv.end(), std::back_inserter(C), &has_changed);
    return C;
}

}
