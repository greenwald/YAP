#include "RelativisticBreitWigner.h"

#include "BlattWeisskopf.h"
#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "DecayChannel.h"
#include "DecayingParticle.h"
#include "FourMomenta.h"
#include "logging.h"
#include "MathUtilities.h"
#include "MeasuredBreakupMomenta.h"
#include "Model.h"
#include "Parameter.h"

namespace yap {

//-------------------------
void RelativisticBreitWigner::addDecayChannel(std::shared_ptr<DecayChannel> c)
{
    auto it = owner()->blattWeisskopfs().find(owner()->quantumNumbers().twoJ() / 2);
    if (it == owner()->blattWeisskopfs().end())
        throw exceptions::Exception("Could not find Blatt-Weisskopf barrier factor in Owner_",
                                    "RelativisticBreitWigner::addDecayChannel");
    BlattWeisskopf_ = it->second;
    // add parameters of BlattWeisskopf to this object
    for (const auto& p : BlattWeisskopf_->parameters())
        addParameter(p);
    addParameter(owner()->radialSize());
    addParameter(mass());
}

//-------------------------
void RelativisticBreitWigner::checkDecayChannel(const DecayChannel& c) const
{
    // check channel is only to two daughters
    if (c.daughters().size() != 2)
        throw exceptions::Exception("Wrong number of daughters (" + std::to_string(c.daughters().size()) + " != 2",
                                    "RelativisticBreitWigner::checkDecayChannel");

    // check channel is to spin-0 particles
    auto sv = spins(c.daughters());
    if (std::any_of(sv.begin(), sv.end(), std::bind(std::not_equal_to<unsigned>(), std::placeholders::_1, 0)))
        throw exceptions::Exception("Non-spin-0 daughter", "RelativisticBreitWigner::checkDecayChannel");
}

//-------------------------
void RelativisticBreitWigner::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    // if no calculation necessary, exit
    if (D.status(*T(), si) != CalculationStatus::uncalculated)
        return;

    /////////////////////////
    // common factors:

    // radial size squared
    double r2 = pow(owner()->radialSize()->value(), 2);

    // squared resonance mass
    double m2_R = pow(mass()->value(), 2);

    // i * mass^2 * nominal width
    auto im2w_R = 1_i * m2_R * width()->value();

    // nominal mass * nominal width
    auto MG = mass()->value() * width()->value();

    // J + 1/2
    unsigned twoLp1 = BlattWeisskopf_->L() + 1;

    /////////////////////////

    // calculate BlattWeisskopf
    BlattWeisskopf_->calculate(D);

    // T := 1 / (M^2 - m^2 - i * M * Gamma)
    for (auto& d : D) {

        // retrieve masses
        double m_ab = model()->fourMomenta()->m(d, pc);
        double m2_ab = m_ab * m_ab;
        double m_a  = model()->fourMomenta()->m(d, pc->daughters()[0]);
        double m_b  = model()->fourMomenta()->m(d, pc->daughters()[1]);

        // calculate nominal breakup momentum
        double q2_nomi = measured_breakup_momenta::q2(m2_R, m_a, m_b);

        // retrieve measured breakup momentum
        double q2_meas = measured_breakup_momenta::q2(m2_ab, m_a, m_b);

        double Q = sqrt(q2_meas / q2_nomi);

        auto imw = im2w_R / m_ab * pow(Q, twoLp1) * pow(BlattWeisskopf_->value(d, pc), 2) / squared_barrier_factor(BlattWeisskopf_->L(), q2_nomi * r2);

        T()->setValue(MG / (m2_R - m2_ab - imw), d, si, D);

    }

    D.status(*T(), si) = CalculationStatus::calculated;
}

}




