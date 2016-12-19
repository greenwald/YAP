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
    if (!owner()) {
        FLOG(ERROR) << "Owner isn't set";
        C &= false;
    }

    return C;
}

//-------------------------
void MassShape::setOwner(DecayingParticle* dp)
{
    if (owner())
        throw exceptions::Exception("MassShape already has owner", "MassShape::setOwner");

    Owner_ = dp;
}

//-------------------------
DecayTreeVector& MassShape::ownersDecayTrees()
{
    if (!owner())
        throw exceptions::Exception("MassShape has no owner", "MassShape::ownersDecayTrees");
    return owner()->DecayTrees_;
}

//-------------------------
void MassShape::addAmplitudeComponent(const AmplitudeComponent& ac, DecayTree& dt) const
{
    if (!owner())
        throw exceptions::Exception("MassShape has no owner", "MassShape::addAmplitudeComponent");
    owner()->addAmplitudeComponent(ac, dt);
}

//-------------------------
void MassShape::removeAmplitudeComponent(const AmplitudeComponent& ac, DecayTree& dt) const
{
    if (!owner())
        throw exceptions::Exception("MassShape has no owner", "MassShape::removeAmplitudeComponent");
    owner()->removeAmplitudeComponent(ac, dt);
}

//-------------------------
void MassShape::replaceFreeAmplitude(std::shared_ptr<FreeAmplitude> ac, DecayTree& dt) const
{
    if (!owner())
        throw exceptions::Exception("MassShape has no owner", "MassShape::replaceFreeAmplitude");
    owner()->replaceFreeAmplitude(ac, dt);
}

//-------------------------
const Model* MassShape::model() const
{
    return owner()->model();
}

}
