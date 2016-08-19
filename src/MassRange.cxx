#include "MassRange.h"

#include "container_utils.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "MassAxes.h"
#include "Parameter.h"
#include "ParticleCombination.h"

#include <algorithm>
#include <functional>

namespace yap {

//-------------------------
const FinalStateParticleVector::value_type& fsp_with_index(const FinalStateParticleVector& FSPs, unsigned i)
{
    for (const auto& fsp : FSPs)
        for (const auto& pc : fsp->particleCombinations())
            if (pc->indices()[0] == i)
                return fsp;
    throw exceptions::Exception("fsp with index " + std::to_string(i) + " not found", "fsp_with_index");
}

//-------------------------
const MassRange mass_range(double isp_mass, const std::shared_ptr<ParticleCombination>& pc, const FinalStateParticleVector& FSPs)
{
    if (isp_mass < 0)
        throw exceptions::Exception("isp_mass is negative", "mass_range");

    if (std::any_of(pc->indices().begin(), pc->indices().end(), [&](unsigned i){return i >= FSPs.size();}))
        throw exceptions::Exception("pc contains index beyond size of FSP", "mass_range");

    MassRange m = {0, isp_mass};

    for (size_t i = 0; i < FSPs.size(); ++i) {
        if (std::find(pc->indices().begin(), pc->indices().end(), i) != pc->indices().end())
            // add mass to low end
            m[0] += FSPs[i]->mass();
        else
            // subtract mass from high end
            m[1] -= FSPs[i]->mass();
    }

    return m;
}

//-------------------------
const std::vector<MassRange> mass_range(double isp_mass, const MassAxes& A, const FinalStateParticleVector& FSPs)
{
    std::vector<MassRange> R;
    R.reserve(A.size());
    std::transform(A.begin(), A.end(), std::back_inserter(R), [&](const MassAxes::value_type & a) {return mass_range(isp_mass, a, FSPs);});
    return R;
}

//-------------------------
const MassRange squared(MassRange mr)
{
    std::transform(mr.begin(), mr.end(), mr.begin(), std::bind(pow, std::placeholders::_1, 2));
    return mr;
}

//-------------------------
const std::vector<MassRange> squared(std::vector<MassRange> mr)
{
    std::transform(mr.begin(), mr.end(), mr.begin(), [](const MassRange & r) {return squared(r);});
    return mr;
}

}
