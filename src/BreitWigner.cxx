#include "BreitWigner.h"

#include "CachedValue.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleFactory.h"

namespace yap {

//-------------------------
BreitWigner::BreitWigner(double m, double w) :
    MassShapeWithNominalMass(m),
    Width_(std::make_shared<RealParameter>(w))
{
    addParameter(Width_);
}

//-------------------------
BreitWigner::BreitWigner(const ParticleTableEntry& pte) :
    BreitWigner(pte.mass(), get_nth_entry(pte, 0, "BreitWigner::BreitWigner"))
{
}

//-------------------------
void BreitWigner::calculateT(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    // common factor := M^2 - i * M * Gamma
    auto M2_iMG = pow(mass()->value(), 2) - 1_i * mass()->value() * Width_->value();

    // T := 1 / (M^2 - m^2 - i * M * Gamma)
    for (auto& d : D)
        T()->setValue(1. / (M2_iMG - model()->fourMomenta()->m2(d, pc)), d, si, D);
}

//-------------------------
bool BreitWigner::consistent() const
{
    bool C = MassShapeWithNominalMass::consistent();

    if (Width_->value() <= 0) {
        FLOG(ERROR) << "width <= 0";
        C &= false;
    }

    return C;
}

}




