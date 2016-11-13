#include "MassShape.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "DecayingState.h"
#include "Exceptions.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "VariableStatus.h"
#include "logging.h"

namespace yap {

//-------------------------
MassShape::MassShape() :
    RecalculableAmplitudeComponent(equal_by_orderless_content),
    T_(ComplexCachedValue::create(*this))
{}

//-------------------------
void MassShape::calculate(DataPartition& D) const
{
    // loop over (ParticleCombination --> symmetrization index) map
    for (const auto& pc_si : symmetrizationIndices()) {

        // recalculate & cache, if necessary
        if (D.status(*T(), pc_si.second) == CalculationStatus::uncalculated) {

            calculateT(D, pc_si.first, pc_si.second);
            D.status(*T(), pc_si.second) = CalculationStatus::calculated;
        }

    }
}

//-------------------------
void MassShape::updateCalculationStatus(StatusManager& D) const
{
    if (status() == VariableStatus::changed)
        D.set(*T(), CalculationStatus::uncalculated);
}

//-------------------------
const std::complex<double> MassShape::value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
{
    return T_->value(d, symmetrizationIndex(pc));
}

//-------------------------
bool MassShape::consistent() const
{
    bool C = DataAccessor::consistent();

    // check if owning DecayingState is set
    if (!DecayingState_) {
        FLOG(ERROR) << "Owning DecayingState isn't set";
        C &= false;
    }

    // check owning DecayingState's mass shape
    if (DecayingState_ and DecayingState_->massShape().get() != this) {
        FLOG(ERROR) << "Owning DecayingState's mass shape is not this MassShape";
        C &= false;
    }
    
    return C;
}

//-------------------------
void MassShape::setDecayingState(DecayingState* r)
{
    if (DecayingState_)
        throw exceptions::Exception("MassShape already has owning DecayingState", "MassShape::setDecayingState");

    DecayingState_ = r;
}

//-------------------------
const Model* MassShape::model() const
{
    return DecayingState_->model();
}

}
