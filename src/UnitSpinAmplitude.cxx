#include "UnitSpinAmplitude.h"

#include "ClebschGordan.h"
#include "Exceptions.h"
#include "ParticleCombination.h"
#include "Spin.h"

namespace yap {

//-------------------------
UnitSpinAmplitude::UnitSpinAmplitude(Model& m, unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s) :
    SpinAmplitude(m, two_J, two_j, l, two_s, equal_always)
{
    // check if all spins are 0 for non-res decays
    if (two_j.size() > 2 and (two_J != 0 or l != 0 or two_s != 0 or !all_zero(two_j)))
        throw exceptions::Exception("More than 2 daughters, but not all spins are 0", "UnitSpinAmplitude::UnitSpinAmplitude");
    
    // for 2 daughters, loop over all possibilities
    for (auto two_M : projections(two_J)) {
        for (const auto& two_m : projections(two_j)) {
            if (two_m.size() > 2) {
                addAmplitude(two_M, two_m, true);
            } else {
                try {
                    if (ClebschGordan::nonzeroCoupling(two_j[0], two_m[0], two_j[1], two_m[1], l, two_s, two_J))
                        addAmplitude(two_M, two_m, true);
                } catch (const exceptions::InconsistentSpinProjection&) { /* ignore */ }
            }
        }
    }
}

}
