#include "PoleMass.h"

#include "CalculationStatus.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleFactory.h"
#include "Resonance.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
PoleMass::PoleMass(std::complex<double> mass) :
    MassShape(),
    Mass_(std::make_shared<ComplexParameter>(mass))
{
    addParameter(Mass_);
}

//-------------------------
void PoleMass::setParameters(const ParticleTableEntry& entry)
{
    // copy current value
    auto m = Mass_->value();

    if (real(m) < 0)
        m.real(entry.Mass);
    
    if (imag(m) < 0 and !entry.MassShapeParameters.empty())
        m.imag(entry.MassShapeParameters[0] / 2.);

    *Mass_ = m;
}

//-------------------------
void PoleMass::calculateT(DataPartition& D, const std::shared_ptr<ParticleCombination>& pc, unsigned si) const
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




