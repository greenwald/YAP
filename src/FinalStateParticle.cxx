#include "FinalStateParticle.h"

#include "logging.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
FinalStateParticle::FinalStateParticle(const QuantumNumbers& q, double m, std::string name)
    : Particle(q, m, name)
{
    // final state particles have fixed mass
    mass()->setVariableStatus(kFixed);
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
void FinalStateParticle::addSymmetrizationIndex(std::shared_ptr<const ParticleCombination> c)
{
    // check if not yet there
    if (std::find(SymmetrizationIndices_.begin(), SymmetrizationIndices_.end(), c) == SymmetrizationIndices_.end()) {
        SymmetrizationIndices_.push_back(c);
    } else {
        LOG(WARNING) << "FinalStateParticle::addSymmetrizationIndex() - Index already existing!";
    }
}

//-------------------------
void FinalStateParticle::setSymmetrizationIndexParents()
{
    ParticleCombinationVector PCs = SymmetrizationIndices_;

    // check if already set
    if (PCs[0]->parent() != nullptr)
        return;

    SymmetrizationIndices_.clear();

    for (auto& PC : PCs) {
        for (auto& pc : ParticleCombination::particleCombinationSet()) {
            if (ParticleCombination::equivDown(PC, pc)) {
                SymmetrizationIndices_.push_back(pc);
            }
        }
    }

}

}

