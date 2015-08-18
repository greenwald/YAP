#include "FinalStateParticle.h"

#include "logging.h"

namespace yap {


//-------------------------
FinalStateParticle::FinalStateParticle(const QuantumNumbers& q, double mass, std::string name, int pdg, std::vector<ParticleIndex>& indices)
    : Particle(q, mass, name),
      PDGCode_(pdg)
{
    for (ParticleIndex i : indices) {
        this->addSymmetrizationIndex(yap::ParticleCombination::uniqueSharedPtr(i));
    }
}

//-------------------------
bool FinalStateParticle::consistent() const
{
    bool consistent = true;

    consistent &= Particle::consistent();

    if (SymmetrizationIndices_.empty()) {
        LOG(ERROR) << "FinalStateParticle::consistent() - SymmetrizationIndices_ are empty!";
        return false;
    }

    for (auto i : SymmetrizationIndices_) {
        if (i->indices().size() != 1) {
            LOG(ERROR) << "FinalStateParticle::consistent() - SymmetrizationIndices_ don't have size 1!";
            return false;
        }
        if (i->daughters().size() != 0) {
            LOG(ERROR) << "FinalStateParticle::consistent() - SymmetrizationIndices_ have daughters!";
            return false;
        }
    }

    return consistent;
}

//-------------------------
void FinalStateParticle::addSymmetrizationIndex(std::shared_ptr<ParticleCombination> c)
{
    // check if not yet there
    if (std::find(SymmetrizationIndices_.begin(), SymmetrizationIndices_.end(), c) == SymmetrizationIndices_.end()) {
        SymmetrizationIndices_.push_back(c);
    } else {
        LOG(WARNING) << "FinalStateParticle::addSymmetrizationIndex() - Index already existing!";
    }
}

}
