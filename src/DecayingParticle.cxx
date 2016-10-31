#include "DecayingParticle.h"

#include "DecayChannel.h"
#include "DecayTree.h"
#include "FreeAmplitude.h"

#include <memory>

namespace yap {

//-------------------------
const is_of_type<DecayingParticle> is_decaying_particle{};

//-------------------------
std::shared_ptr<DecayChannel> DecayingParticle::addWeakDecay(const ParticleVector& daughters)
{
    return addDecayChannel(std::make_shared<DecayChannel>(daughters), false);
}

//-------------------------
std::shared_ptr<DecayChannel> DecayingParticle::addStrongDecay(const ParticleVector& daughters)
{
    return addDecayChannel(std::make_shared<DecayChannel>(daughters), true);
}

}