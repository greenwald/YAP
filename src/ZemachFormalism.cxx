#include "ZemachFormalism.h"

#include "CachedValue.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "logging.h"
#include "LorentzTransformation.h"
#include "Model.h"
#include "Spin.h"

namespace yap {

//-------------------------
bool equal_zemach(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B)
{
    //check if either empty
    if (!A or !B)
        return false;

    if (A->indices().size() > 3 or B->indices().size() > 3)
        throw exceptions::Exception("Zemach formalism cannot be used with 4 or more particles",
                                    "equal_zemach::operator()");

    // compare shared_ptr addresses
    if (A == B)
        return true;

    // check if both size < 3
    if (A->indices().size() < 3 and B->indices().size() < 3)
        return true;

    // now check if sizes same
    if (A->indices().size() != B->indices().size())
        return false;

    // find resonance and spectator
    auto rA = A->daughters()[0];
    auto sA = A->daughters()[1];
    if (sA->indices().size() == 2)
        std::swap(sA, rA);
    if (rA->indices().size() != 2 and sA->indices().size() != 1)
        throw exceptions::Exception("could not find resonance and spectator in A",
                                    "ParticleCombination::EqualZemach::operator()");
    auto rB = B->daughters()[0];
    auto sB = B->daughters()[1];
    if (sB->indices().size() == 2)
        std::swap(sB, rB);
    if (rB->indices().size() != 2 and sB->indices().size() != 1)
        throw exceptions::Exception("could not find resonance and spectator in B",
                                    "ParticleCombination::EqualZemach::operator()");

    return equal_by_orderless_content(rA, rB) and equal_by_orderless_content(sA, sB);
}

//-------------------------
ZemachSpinAmplitude::ZemachSpinAmplitude(Model& m, unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s) :
    SpinAmplitude(m, two_J, two_j, l, two_s, ParticleCombinationEqualTo(equal_zemach))
{
    if (finalTwoJ().size() != 2)
        throw exceptions::Exception("Wrong number of daughter spins specified (" + std::to_string(finalTwoJ().size()) + " != 2)",
                                    "ZemachSpinAmplitude::ZemachSpinAmplitude");

    // check j1j2S triangle
    if (!triangle(finalTwoJ()[0], finalTwoJ()[1], twoS()))
        throw exceptions::AngularMomentumNotConserved("ZemachSpinAmplitude::ZemachSpinAmplitude");

    if (is_odd(initialTwoJ()) or is_odd(twoS()) or
            std::any_of(finalTwoJ().begin(), finalTwoJ().end(), std::function<bool(int)>(is_odd)))
        throw exceptions::Exception("only supporting integer spins", "ZemachSpinAmplitude::ZemachSpinAmplitude");

    if (initialTwoJ() != 0) {
        if (!all_zero(finalTwoJ()))
            throw exceptions::Exception("Zemach not valid for " + spin_to_string(initialTwoJ())
                                        + " -> " + to_string(finalTwoJ()),
                                        "ZemachSpinAmplitude::ZemachSpinAmplitude");
        if (2 * L() != initialTwoJ() and twoS() != 0)
            throw exceptions::Exception("Zemach not valid for " + spin_to_string(initialTwoJ())
                                        + " -> " + to_string(finalTwoJ())
                                        + " with L = " + std::to_string(L()) + " and S = " + spin_to_string(twoS()),
                                        "ZemachSpinAmplitude::ZemachSpinAmplitude");
    } else {

        if (std::all_of(finalTwoJ().begin(), finalTwoJ().end(), [](const SpinVector::value_type & j) {return j != 0;}))
        throw exceptions::Exception("one daughter must be a scalar", "ZemachSpinAmplitude::ZemachSpinAmplitude");

        // get resonance spin
        auto two_J_r = *std::max_element(finalTwoJ().begin(), finalTwoJ().end());

        if (two_J_r != 2 * L())
            throw exceptions::Exception("orbital angular momentum != resonance spin (" + std::to_string(L()) + " != " + spin_to_string(two_J_r) + ")",
                                        "ZemachSpinAmplitude::ZemachSpinAmplitude");

        if (two_J_r != twoS())
            throw exceptions::Exception("total spin != resonance spin (" + std::to_string(L()) + " != " + spin_to_string(twoS()) + ")",
                                        "ZemachSpinAmplitude::ZemachSpinAmplitude");

    }

    if (twoS() > 2 * 2)
        throw exceptions::Exception("currently only supporting up to spin 2", "ZemachSpinAmplitude::ZemachSpinAmplitude");

    addAmplitude(0, {0, 0});
    // dependencies?
}

//-------------------------
void ZemachSpinAmplitude::addParticleCombination(const ParticleCombination& pc)
{
    if (pc.indices().size() > 2) SpinAmplitude::addParticleCombination(pc);
}

//-------------------------
const std::complex<double> ZemachSpinAmplitude::amplitude(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc,
                                                          int two_M, const SpinProjectionVector& two_m) const
{
    return (twoS() == 0 or pc->indices().size() < 3) ? 1. : SpinAmplitude::amplitude(d, pc, two_M, two_m);
}

//-------------------------
const std::complex<double> ZemachSpinAmplitude::calc(int two_M, const SpinProjectionVector& two_m,
        const DataPoint& d, const StatusManager& sm,
        const std::shared_ptr<const ParticleCombination>& pc) const
{
    if (pc->indices().size() < 3)
        return 1.;

    if (pc->indices().size() > 3)
        throw exceptions::Exception("Zemach not valid for more than 3 particles",
                                    "ZemachSpinAmplitude::calc");

    if (twoS() == 0)
        return 1.;

    // find resonance and spectator ParticleCombinations
    auto pcR = pc->daughters()[0];
    auto pcS = pc->daughters()[1];
    if (pcR->indices().size() == 1)
        std::swap(pcR, pcS);

    // get resonance four-momentum
    auto R4 = model()->fourMomenta()->p(d, pcR);
    // get spectator four-momentum
    auto p4 = model()->fourMomenta()->p(d, pcS);
    // get four-momentum of R's daughter
    auto q4 = model()->fourMomenta()->p(d, pcR->daughters()[0]);

    // boost p and q into R rest frame, and get spacial components
    auto B = lorentzTransformation<double>(-R4);
    auto p = vect(B * p4);
    auto q = vect(B * q4);

    if (twoS() == 2)
        return -2. * (p * q);

    // else twoS() == 4
    return 4 * (pow(p * q, 2) - norm(p) * norm(q) / 3.);
}

}
