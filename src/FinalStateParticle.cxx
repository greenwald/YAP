#include "FinalStateParticle.h"

#include "Exceptions.h"
#include "logging.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"

namespace yap {

//-------------------------
const is_of_type<FinalStateParticle> is_final_state_particle{};

//-------------------------
FinalStateParticle::FinalStateParticle(const std::string& name, const QuantumNumbers& q, double m) :
    Particle(name, q),
    Model_(nullptr),
    Mass_(m)
{
    if (Mass_ < 0)
        throw exceptions::Exception("Mass is negative", "FinalStateParticle::FinalStateParticle");
}

//-------------------------
FinalStateParticle::FinalStateParticle(const ParticleTableEntry& pte)
    : FinalStateParticle(pte.name(), pte.quantumNumbers(), pte.mass())
{
}

//-------------------------
void FinalStateParticle::addParticleCombination(const ParticleCombination& pc)
{
    // pc must be final state particle
    if (!is_final_state_particle_combination(pc))
        throw exceptions::Exception("pc is not final state particle", "FinalStateParticle::addParticleCombination");

    // if pc not already in particleCombinations vector, add it
    if (particleCombinations().find(pc.shared_from_this()) == particleCombinations().end())
        Particle::addParticleCombination(pc);
}

}

