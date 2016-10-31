#include "FreeAmplitude.h"

#include "DecayChannel.h"
#include "container_utils.h"
#include "Exceptions.h"
#include "Model.h"
#include "Particle.h"
#include "ParticleCombination.h"
#include "Spin.h"
#include "SpinAmplitude.h"

namespace yap {

//-------------------------
FreeAmplitude::FreeAmplitude(std::shared_ptr<DecayChannel> dc, std::shared_ptr<SpinAmplitude> sa,
                             std::complex<double> a) :
    ComplexParameter(a),
    DecayChannel_(dc),
    SpinAmplitude_(sa)
{
    if (!DecayChannel_)
        throw exceptions::Exception("DecayChannel is nullptr", "FreeAmplitude::FreeAmplitude");

    if (!SpinAmplitude_)
        throw exceptions::Exception("SpinAmpliutde is nullptr", "FreeAmplitude::FreeAmplitude");
}

//-------------------------
const Model* FreeAmplitude::model() const
{
    return (DecayChannel_) ? DecayChannel_->model() : nullptr;
}

//-------------------------
const ParticleCombinationSet& FreeAmplitude::particleCombinations() const
{
    if (!DecayChannel_)
        throw exceptions::Exception("DecayChannel_ is nullptr", "FreeAmplitude::particleCombinations");
    
    return DecayChannel_->particleCombinations();
}

//-------------------------
std::string to_string(const FreeAmplitude& fa)
{
    return to_string(*fa.decayChannel())
        + ", L = " + std::to_string(fa.spinAmplitude()->L())
        + ", S = " + spin_to_string(fa.spinAmplitude()->twoS())
        + (fa.variableStatus() == VariableStatus::fixed ? " [" + to_string(fa.variableStatus())+ "]" : "");
}

}
