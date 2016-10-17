#include "BlattWeisskopf.h"

#include "CalculationStatus.h"
#include "DataPartition.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "logging.h"
#include "MeasuredBreakupMomenta.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
const double squared_barrier_factor(unsigned l, double z)
{
    switch (l) {
    case 0:
        return 1;
    case 1:
        return 2 * z / (z + 1);
    case 2:
        return 13 * pow(z, l) / (z * (z + 3) + 9);
    case 3:
        return 277 * pow(z, l) / (z * (z * (z + 6) + 45) + 225);
    case 4:
        return 12746 * pow(z, l) / (z * (z * (z * (z + 10) + 135) + 1575) + 11025);
    case 5:
        return 998881 * pow(z, l) / (z * (z * (z * (z * (z + 15) + 315) + 6300) + 99225) + 893025);
    case 6:
        return 118394977 * pow(z, l) / (z * (z * (z * (z * (z * (z + 21) + 630) + 18900) + 496125) + 9823275) + 108056025);
    case 7:
        return 19727003738LL * pow(z, l) / (z * (z * (z * (z * (z * (z * (z + 28) + 1134) + 47250) + 1819125) + 58939650) + 1404728325L) + 18261468225LL);
    default:
        /// \todo: speed this up. And check how well it approximates.
        // approximated:
        double num = 0;
        double den = 0;
         for (int n = l; n >= 0; --n) {
             double coef = exp(0.5 * l * n);
             num += coef;
             den += coef / pow(z, n);
         }
         return num / den;
    }
}

//-------------------------
BlattWeisskopf::BlattWeisskopf(unsigned L, const DecayingState& ds) :
    RecalculableAmplitudeComponent(equal_down_by_orderless_content),
    DecayingState_(&ds),
    L_(L)
{
    if (!model())
        throw exceptions::Exception("Model unset", "BlattWeisskopf::BlattWeisskopf");

    if (L_ > 0) {
        addParameter(DecayingState_->radialSize());
        BarrierFactor_ = RealCachedValue::create(*this);
    }

    // if L == 0, values are all always 1, no storage in DataPoint necessary

}

//-------------------------
const std::complex<double> BlattWeisskopf::value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
{
    return (L_ == 0) ? 1 : BarrierFactor_->value(d, symmetrizationIndex(pc));
}

//-------------------------
void BlattWeisskopf::calculate(DataPartition& D) const
{
    // spin 0 always has factors 1, no calculations necessary
    if (L_ == 0)
        return;

    // loop over (ParticleCombination --> symmetrization index) map
    for (const auto& pc_symIndex : symmetrizationIndices()) {

        // check if barrier factor is uncalculated
        if (D.status(*BarrierFactor_, pc_symIndex.second) == CalculationStatus::uncalculated) {

            DEBUG("calculate BlattWeisskopf");

            double r2 = pow(DecayingState_->radialSize()->value(), 2);

            // calculate on all data points in D
            for (auto& d : D) {
                // measured breakup momentum
                double q2_meas = measured_breakup_momenta::q2(d, pc_symIndex.first, *model());
                BarrierFactor_->setValue(sqrt(squared_barrier_factor(L(), q2_meas * r2)), d, pc_symIndex.second, D);
            }

            // update status
            D.status(*BarrierFactor_, pc_symIndex.second) = CalculationStatus::calculated;
        }
    }
}

//-------------------------
void BlattWeisskopf::updateCalculationStatus(StatusManager& D) const
{
    if (status() == VariableStatus::changed)
        D.set(*BarrierFactor_, CalculationStatus::uncalculated);
}

//-------------------------
const Model* BlattWeisskopf::model() const
{
    return DecayingState_->model();
}

//-------------------------
void BlattWeisskopf::addParticleCombination(const ParticleCombination& pc)
{
    if (pc.daughters().size() != 2)
        throw exceptions::NotTwoBodyParticleCombination("cannot calculate Blatt-Weisskopf barrier factor for "
                                                        + std::to_string(pc.daughters().size()) + "-body decay",
                                                        "BlattWeisskopf::addParticleCombination");

    return RecalculableDataAccessor::addParticleCombination(pc);
}

}

