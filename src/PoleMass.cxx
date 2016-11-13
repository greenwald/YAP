#include "PoleMass.h"

#include "CalculationStatus.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleFactory.h"

namespace yap {

//-------------------------
PoleMass::PoleMass(std::complex<double> mass) :
    MassShape(),
    Mass_(std::make_shared<ComplexParameter>(mass))
{
    addParameter(Mass_);
}

//-------------------------
PoleMass::PoleMass(const ParticleTableEntry& pte) :
    PoleMass(std::complex<double>(pte.mass(), get_nth_entry(pte, 0, "PoleMass::PoleMass") / 2))
{
}

//-------------------------
void PoleMass::calculateT(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    auto M2 = pow(Mass_->value(), 2);

    // T := 1 / (M^2 - m^2)
    for (auto& d : D)
        T()->setValue(1. / (M2 - model()->fourMomenta()->m2(d, pc)), d, si, D);
}

//-------------------------
bool PoleMass::consistent() const
{
    bool C = MassShape::consistent();

    if (real(Mass_->value()) <= 0) {
        FLOG(ERROR) << "real(mass) <= 0";
        C &= false;
    }
    if (imag(Mass_->value()) <= 0) {
        FLOG(ERROR) << "imag(mass) <= 0";
        C &= false;
    }

    return C;
}

}




