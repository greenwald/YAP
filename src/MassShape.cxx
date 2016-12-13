#include "MassShape.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "logging.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "VariableStatus.h"

namespace yap {

//-------------------------
MassShape::MassShape() :
    RecalculableAmplitudeComponent(equal_by_orderless_content),
    Owner_(nullptr)
{}

//-------------------------
void MassShape::calculate(DataPartition& D) const
{
    // loop over (ParticleCombination --> symmetrization index) map
    for (const auto& pc_si : symmetrizationIndices())
        calculate(D, pc_si.first, pc_si.second);
}

//-------------------------
bool MassShape::consistent() const
{
    bool C = DataAccessor::consistent();

    // check if owner is set
    if (!Owner_) {
        FLOG(ERROR) << "Owner isn't set";
        C &= false;
    }

    return C;
}

//-------------------------
void MassShape::setOwner(DecayingParticle* dp)
{
    if (Owner_)
        throw exceptions::Exception("MassShape already has owner", "MassShape::setOwner");

    Owner_ = dp;
}

//-------------------------
DecayTreeVector& MassShape::ownersDecayTrees()
{
    if (!Owner_)
        throw exceptions::Exception("MassShape has no owner", "MassShape::ownersDecayTrees");
    return Owner_->DecayTrees_;
}

//-------------------------
void MassShape::addAmplitudeComponent(const AmplitudeComponent& ac, DecayTree& dt) const
{
    if (!Owner_)
        throw exceptions::Exception("MassShape has no owner", "MassShape::addAmplitudeComponent");
    Owner_->addAmplitudeComponent(ac, dt);
}

//-------------------------
void MassShape::removeAmplitudeComponent(const AmplitudeComponent& ac, DecayTree& dt) const
{
    if (!Owner_)
        throw exceptions::Exception("MassShape has no owner", "MassShape::removeAmplitudeComponent");
    Owner_->removeAmplitudeComponent(ac, dt);
}

//-------------------------
void MassShape::replaceFreeAmplitude(std::shared_ptr<FreeAmplitude> ac, DecayTree& dt) const
{
    if (!Owner_)
        throw exceptions::Exception("MassShape has no owner", "MassShape::replaceFreeAmplitude");
    Owner_->replaceFreeAmplitude(ac, dt);
}

//-------------------------
const Model* MassShape::model() const
{
    return Owner_->model();
}

}
