#include "Particle.h"

#include "logging.h"
#include "Model.h"
#include "Parameter.h"

namespace yap {

//-------------------------
Particle::Particle(const QuantumNumbers& q, std::string name) :
    std::enable_shared_from_this<Particle>(),
    QuantumNumbers_(q),
    Name_(name)
{}

//-------------------------
bool Particle::consistent() const
{
    bool C = true;

    if (ParticleCombinations_.empty()) {
        FLOG(ERROR) << "ParticleCombinations_ is empty";
        C &= false;
    }

    return C;
}

//-------------------------
void Particle::addParticleCombination(const std::shared_ptr<ParticleCombination>& pc)
{
    ParticleCombinations_.insert(pc);
}

//-------------------------
void Particle::pruneParticleCombinations()
{
    prune_particle_combinations(ParticleCombinations_);
}

//-------------------------
const SpinVector spins(const ParticleVector& v)
{
    SpinVector s;
    s.reserve(v.size());
    std::transform(v.begin(), v.end(), std::back_inserter(s),
    [](const ParticleVector::value_type & p) {return p->quantumNumbers().twoJ();});
    return s;
}

//-------------------------
const bool decays_to_full_final_state(const Particle& p)
{
    return std::any_of(p.particleCombinations().begin(), p.particleCombinations().end(),
                       [&p](const std::shared_ptr<ParticleCombination>& pc)
                       {return pc->indices().size() == p.model()->finalStateParticles().size();});
}

}
