#include "DecayingParticle.h"

#include "DecayChannel.h"
#include "DecayTree.h"
#include "FreeAmplitude.h"
#include "VariableStatus.h"

#include "logging.h"

namespace yap {

//-------------------------
std::shared_ptr<DecayChannel> DecayingParticle::addChannel(const ParticleVector& daughters)
{
    FLOG(INFO) << "in";
    return addChannel(std::make_shared<DecayChannel>(daughters));
}

//-------------------------
void DecayingParticle::fixSolitaryFreeAmplitudes()
{
    // // loop over entries in map of (spin projection) -> (decay tree vector)
    // for (auto& m_dtv : decayTrees())
    //     // if only available decay tree
    //     if (m_dtv.second.size() == 1)
    //         m_dtv.second[0]->freeAmplitude()->variableStatus() = VariableStatus::fixed;
    // for (auto& c : channels()) {
    //     for (auto& d : c->daughters()) {
    //         if (d and is_decaying_particle(*d))
    //             std::static_pointer_cast<DecayingParticle>(d)->fixSolitaryFreeAmplitudes();
    //     }
    // }
}

}
