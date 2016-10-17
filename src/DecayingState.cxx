#include "DecayingState.h"

#include "BlattWeisskopf.h"
#include "Parameter.h"

namespace yap {

//-------------------------
DecayingState::DecayingState(const std::string& name, const QuantumNumbers& q, double radial_size) :
    Particle(name, q),
    RadialSize_(std::make_shared<PositiveRealParameter>(radial_size))
{
}

//-------------------------
void DecayingState::addParticleCombination(BlattWeisskopf& bw, const ParticleCombination& pc) const
{
    bw.addParticleCombination(pc);
}

//-------------------------
void DecayingState::registerWithModel(BlattWeisskopf& bw) const
{
    bw.registerWithModel();
}

}
