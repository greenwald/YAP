#include "Resonance.h"

#include "DecayTree.h"
#include "Exceptions.h"
#include "logging.h"
#include "MassShape.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
const is_of_type<Resonance> is_resonance{};

//-------------------------
Resonance::Resonance(const std::string& name, const QuantumNumbers& q, double radialSize, std::shared_ptr<MassShape> massShape) :
    DecayingParticle(name, q, radialSize),
    MassShape_(massShape)
{
    if (!MassShape_)
        throw exceptions::Exception("MassShape unset", "Resonance::Resonance");

    MassShape_->setResonance(this);
}

//-------------------------
std::shared_ptr<DecayChannel> Resonance::addDecayChannel(std::shared_ptr<DecayChannel> c, bool conserve_parity)
{
    MassShape_->checkDecayChannel(*c);
    auto dc = DecayingParticle::addDecayChannel(c, conserve_parity);
    MassShape_->addDecayChannel(c);
    return dc;
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
void Resonance::addParticleCombination(const ParticleCombination& pc)
{
    DecayingParticle::addParticleCombination(pc);
    MassShape_->addParticleCombination(pc);
}

//-------------------------
void Resonance::modifyDecayTree(DecayTree& dt)
{
    DecayingParticle::modifyDecayTree(dt);
    addAmplitudeComponent(dt, *MassShape_);
}

}
