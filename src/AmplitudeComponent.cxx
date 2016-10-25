#include "AmplitudeComponent.h"

#include "Attributes.h"
#include "Parameter.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
const bool StaticAmplitudeComponent::validFor(const ParticleCombination& pc) const
{
    return symmetrizationIndices().find(pc.shared_from_this()) != symmetrizationIndices().end();
}

//-------------------------
const bool RecalculableAmplitudeComponent::validFor(const ParticleCombination& pc) const
{
    return symmetrizationIndices().find(pc.shared_from_this()) != symmetrizationIndices().end();
}

//-------------------------
const VariableStatus RecalculableAmplitudeComponent::status() const
{
    static variable_status S;
    return S(static_cast<const RecalculableDataAccessor&>(*this));
}

}
