#include "ModelIntegral.h"

#include "Attributes.h"
#include "DecayingParticle.h"
#include "DecayTreeVectorIntegral.h"
#include "Exceptions.h"
#include "Filter.h"
#include "Model.h"
#include "Parameter.h"

#include <memory>

namespace yap {

//-------------------------
ModelIntegral::ModelIntegral(const Model& model)
{
    // for each initial state particle
    for (const auto& isp_mix : model.initialStateParticles())
        // for each spin projection
        for (const auto& m_b : isp_mix.second)
            // create new DecayTreeVectorIntegral
            Integrals_.emplace(m_b.second, DecayTreeVectorIntegral(filter(isp_mix.first->decayTrees(), m_equals(m_b.first))));
}

//-------------------------
const RealIntegralElement integral(const ModelIntegral& MI)
{
    return std::accumulate(MI.integrals().begin(), MI.integrals().end(), RealIntegralElement(),
                           [](RealIntegralElement& i, const IntegralMap::value_type& b_I)
                           { return i += b_I.first->value() * integral(b_I.second); });
}

}
