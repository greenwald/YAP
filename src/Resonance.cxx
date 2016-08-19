#include "Resonance.h"

#include "DecayTree.h"
#include "Exceptions.h"
#include "logging.h"
#include "MassShape.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
Resonance::Resonance(const QuantumNumbers& q, std::string name, double radialSize, std::shared_ptr<MassShape> massShape) :
    DecayingParticle(q, name, radialSize),
    MassShape_(massShape)
{
    if (!MassShape_)
        throw exceptions::Exception("MassShape unset", "Resonance::Resonance");

    MassShape_->setResonance(this);
}

//-------------------------
std::shared_ptr<DecayChannel> Resonance::addChannel(std::shared_ptr<DecayChannel> c)
{
    auto dc = DecayingParticle::addChannel(c);
    MassShape_->addDecayChannel(dc);
    return dc;
}

//-------------------------
void Resonance::checkDecayChannel(const DecayChannel& c) const
{
    DecayingParticle::checkDecayChannel(c);
    MassShape_->checkDecayChannel(c);
}

//-------------------------
bool Resonance::consistent() const
{
    bool C = DecayingParticle::consistent();

    if (!MassShape_->consistent()) {
        FLOG(ERROR) << "MassShape not consistent";
        C &= false;
    }
    if (MassShape_->resonance() != this) {
        FLOG(ERROR) << "MassShape's resonance is not this";
        C &= false;
    }

    return C;
}

//-------------------------
void Resonance::registerWithModel()
{
    DecayingParticle::registerWithModel();
    MassShape_->registerWithModel();
}

//-------------------------
void Resonance::addParticleCombination(const std::shared_ptr<ParticleCombination>& c)
{
    DecayingParticle::addParticleCombination(c);
    MassShape_->addParticleCombination(c);
}

//-------------------------
void Resonance::modifyDecayTree(DecayTree& dt) const
{
    DecayingParticle::modifyDecayTree(dt);
    dt.addRecalculableDataAccessor(*MassShape_);
}

//-------------------------
const bool is_resonance(const std::shared_ptr<Particle>& p)
{
    return std::dynamic_pointer_cast<Resonance>(p) != nullptr;
}

}
