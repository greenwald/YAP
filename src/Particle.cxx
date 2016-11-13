#include "Particle.h"

#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleFactory.h"

namespace yap {

//-------------------------
Particle::Particle(const ParticleTableEntry& pte)
    : Particle(pte.name(), pte.quantumNumbers())
{
}

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
void Particle::addParticleCombination(const ParticleCombination& pc)
{
    ParticleCombinations_.insert(pc.shared_from_this());
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
                       [&p](const ParticleCombinationSet::value_type& pc)
                       {return pc->indices().size() == p.model()->finalStateParticles().size();});
}

//-------------------------
std::string to_string(const ParticleVector& p, const SpinProjectionVector& two_m)
{
    if (p.size() != two_m.size())
        throw exceptions::Exception("vector size mismatch", "two_string(ParticleVector, SpinProjectionVector)");
    std::string s;
    for (size_t i = 0; i < p.size(); ++i)
        s += ", " + p[i]->name() + " [m = " + spin_to_string(two_m[i]) + "]";
    return s.erase(0, 2);
}

}
