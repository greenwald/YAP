#include "Attributes.h"

#include "BlattWeisskopf.h"
#include "DecayChannel.h"
#include "DecayTree.h"
#include "DecayingParticle.h"
#include "Filter.h"
#include "FinalStateParticle.h"
#include "FreeAmplitude.h"
#include "MassShape.h"
#include "MassShapeWithNominalMass.h"
#include "Model.h"
#include "Particle.h"
#include "RecalculableDataAccessor.h"
#include "Resonance.h"
#include "SpinAmplitude.h"
#include "container_utils.h"

#include <numeric>

namespace yap {

//-------------------------
const unsigned orbital_angular_momentum::operator()(const SpinAmplitude& sa) const
{
    return sa.L();
}

//-------------------------
const unsigned orbital_angular_momentum::operator()(const FreeAmplitude& fa) const
{
    return operator()(fa.spinAmplitude());
}

//-------------------------
const unsigned orbital_angular_momentum::operator()(const BlattWeisskopf& bw) const
{
    return bw.L();
}

//-------------------------
const unsigned orbital_angular_momentum::operator()(const DecayTree& dt) const
{
    return operator()(dt.freeAmplitude());
}

//-------------------------
const unsigned spin_angular_momentum::operator()(const SpinAmplitude& sa) const
{
    return sa.twoS();
}

//-------------------------
const unsigned spin_angular_momentum::operator()(const FreeAmplitude& fa) const
{
    return operator()(fa.spinAmplitude());
}

//-------------------------
const unsigned spin_angular_momentum::operator()(const Particle& p) const
{
    return p.quantumNumbers().twoJ();
}

//-------------------------
const unsigned spin_angular_momentum::operator()(const DecayTree& dt) const
{
    return operator()(dt.freeAmplitude());
}

//-------------------------
const int spin_projection::operator()(const DecayTree& dt) const
{
    return dt.initialTwoM();
}

//-------------------------
const int charge::operator()(const QuantumNumbers& q) const
{
    return q.Q();
}

//-------------------------
const int charge::operator()(const Particle& p) const
{
    return operator()(p.quantumNumbers());
}

//-------------------------
const int charge::operator()(const DecayChannel& dc) const
{
    return std::accumulate(dc.daughters().begin(), dc.daughters().end(), 0,
                           [&](int q, const ParticleVector::value_type& p)
                           { return q += (*this)(p); });
}

//-------------------------
const bool to::operator()(const DecayChannel& dc) const
{
    return std::all_of(objects().begin(), objects().end(),
                       [&](const Particle* const p)
                       { return std::count(objects().begin(), objects().end(), p)
                               <= std::count_if(dc.daughters().begin(), dc.daughters().end(),
                                                [&](const std::shared_ptr<Particle>& d)
                                                { return d.get() == p; });});
}

//-------------------------
const bool to::operator()(const FreeAmplitude& fa) const
{
    return operator()(fa.decayChannel());
}

//-------------------------
const bool to::operator()(const DecayTree& dt) const
{
    return operator()(dt.freeAmplitude());
}

//-------------------------
const bool to::operator()(const Particle& p) const
{
    if (!is_decaying_particle(p))
        return false;
    for (const auto& dc : static_cast<const DecayingParticle*>(&p)->channels())
        if (operator()(dc))
            return true;
    return false;
}

//-------------------------
const bool exactly_to::operator()(const DecayChannel& dc) const
{
    return std::all_of(objects().begin(), objects().end(),
                       [&](const Particle* const p)
                       { return std::count(objects().begin(), objects().end(), p)
                               == std::count_if(dc.daughters().begin(), dc.daughters().end(),
                                                [&](const ParticleVector::value_type& d)
                                                { return d.get() == p; });});
}

//-------------------------
const bool is_fixed::operator()(const ParameterBase& p) const
{
    return p.variableStatus() == VariableStatus::fixed;
}

//-------------------------
const bool is_fixed::operator()(const DecayTree& dt) const
{
    for (const auto& fa : free_amplitudes(dt))
        if (!operator()(fa))
            return false;
    return true;
}

//-------------------------
const bool is_not_fixed::operator()(const ParameterBase& p) const
{
    return p.variableStatus() != VariableStatus::fixed;
}

//-------------------------
const bool is_not_fixed::operator()(const DecayTree& dt) const
{
    for (const auto& fa : free_amplitudes(dt))
        if (operator()(fa))
            return true;
    return false;
}

//-------------------------
const VariableStatus variable_status::operator()(const ParameterBase& p) const
{
    return p.variableStatus();
}

//-------------------------
const VariableStatus variable_status::operator()(const RecalculableDataAccessor& rda) const
{
    return operator()(rda.parameters());
}

//-------------------------
const VariableStatus variable_status::operator()(const AmplitudeComponent& a) const
{
    return a.status();
}

//-------------------------
const bool has_free_amplitude::operator()(const DecayTree& dt) const
{
    for (const auto& fa : free_amplitudes(dt))
        if (std::find(objects().begin(), objects().end(), fa.get()) != objects().end())
            return true;
    return false;
}

//-------------------------
const bool has_free_amplitude::operator()(const Particle& p) const
{
    if (!is_decaying_particle(p))
        return false;
    for (const auto& dt : dynamic_cast<const DecayingParticle&>(p).decayTrees())
        if(dt and std::find(objects().begin(), objects().end(), dt->freeAmplitude().get()) != objects().end())
            return true;
    return false;
}

//-------------------------
const bool has_decay_tree::operator()(const Particle& p) const
{
    if (!is_decaying_particle(p))
        return false;
    for (const auto& dt : dynamic_cast<const DecayingParticle&>(p).decayTrees())
        if (dt and std::find(objects().begin(), objects().end(), dt.get()) != objects().end())
            return true;
    return false;
}

//-------------------------
const bool has_decay_tree::operator()(const DecayTree& dt) const
{
    for (const auto& i_dt : dt.daughterDecayTrees()) {
        if (std::find(objects().begin(), objects().end(), i_dt.second.get()) != objects().end())
            return true;
        if (operator()(i_dt.second))
            return true;
    }
    return false;
}

//-------------------------
const bool has_decay_channel::operator()(const Particle& p) const
{
    if (!is_decaying_particle(p))
        return false;
    for (const auto& dc : dynamic_cast<const DecayingParticle&>(p).channels())
        if (std::find(objects().begin(), objects().end(), dc.get()) != objects().end())
            return true;
    return false;
}

//-------------------------
const bool has_decay_channel::operator()(const FreeAmplitude& fa) const
{
    return std::find(objects().begin(), objects().end(), fa.decayChannel().get()) != objects().end();
}

//-------------------------
const bool has_decay_channel::operator()(const DecayTree& dt) const
{
    for (const auto& fa : free_amplitudes(dt))
        if (operator()(fa))
            return true;
    return false;
}

//-------------------------
const bool has_particle_combination::operator()(const ParticleCombinationSet& S) const
{
    for (const auto& pc : S)
        if (std::find(objects().begin(), objects().end(), pc.get()) != objects().end())
            return true;
    return false;
}

//-------------------------
const bool has_particle_combination::operator()(const DecayChannel& dc) const
{
    return operator()(dc.particleCombinations());
}

//-------------------------
const bool has_particle_combination::operator()(const FreeAmplitude& fa) const
{
    return operator()(fa.decayChannel());
}

//-------------------------
const bool has_particle_combination::operator()(const DecayTree& dt) const
{
    return operator()(dt.freeAmplitude());
}

//-------------------------
const bool has_particle_combination::operator()(const Particle& p) const
{
    return operator()(p.particleCombinations());
}

//-------------------------
std::shared_ptr<const DecayingState> parent_state::operator()(const DecayTree& dt) const
{
    if (!dt.model())
        throw exceptions::Exception("model is nullptr", "parent_state::operator()(DecayTree)");
    return std::static_pointer_cast<DecayingState>(particle(*dt.model(), has_decay_tree(dt)));
}

//-------------------------
std::shared_ptr<const DecayingState> parent_state::operator()(const FreeAmplitude& fa) const
{
    if (!fa.model())
        throw exceptions::Exception("model is nullptr", "parent_state::operator()(FreeAmplitude)");
    return std::static_pointer_cast<DecayingState>(particle(*fa.model(), has_free_amplitude(fa)));
}

//-------------------------
std::shared_ptr<const DecayingState> parent_state::operator()(const DecayChannel& dc) const
{
    if (!dc.model())
        throw exceptions::Exception("model is nullptr", "parent_state::operator()(DecayChannel)");
    return std::static_pointer_cast<DecayingState>(particle(*dc.model(), has_decay_channel(dc)));
}

//-------------------------
std::shared_ptr<const DecayingState> parent_state::operator()(const BlattWeisskopf& bw) const
{
    if (!bw.decayingState())
        return nullptr;
    return std::static_pointer_cast<const DecayingState>(bw.decayingState()->shared_from_this());
}

//-------------------------
std::shared_ptr<const DecayingState> parent_state::operator()(const MassShape& m) const
{
    if (!m.resonance())
        return nullptr;
    return std::static_pointer_cast<DecayingState>(m.resonance()->shared_from_this());
}

//-------------------------
const bool has_a_mass::operator()(const MassShape& m) const
{
    return dynamic_cast<const MassShapeWithNominalMass*>(&m) != nullptr;
}

//-------------------------
const bool has_a_mass::operator()(const Particle& p) const
{
    return dynamic_cast<const FinalStateParticle*>(&p) != nullptr
        or (dynamic_cast<const Resonance*>(&p)
            and operator()(dynamic_cast<const Resonance&>(p).massShape()));
}

//-------------------------
const RealParameter& mass_parameter::operator()(const MassShape& m) const
{
    if (!dynamic_cast<const MassShapeWithNominalMass*>(&m))
        throw exceptions::Exception("Mass shape does not inherit from MassShapeWithNominalMass", "mass_parameter::operator()(MassShape)");
    auto mp = dynamic_cast<const MassShapeWithNominalMass&>(m).mass();
    if (!mp)
        throw exceptions::Exception("Mass parameter is nullptr", "mass_parameter::operator()(MassShape)");
    return *mp;
}

//-------------------------
const RealParameter& mass_parameter::operator()(const Particle& p) const
{
    if (!dynamic_cast<const Resonance*>(&p))
        throw exceptions::Exception("Particle is not Resonance", "mass_parameter::operator()(Particle)");
    return operator()(dynamic_cast<const Resonance&>(p).massShape());
}
    
//-------------------------
const bool from::operator()(const Particle& p) const
{
    if (!p.model())
        throw exceptions::Exception("model is nullptr", "from::operator()(Particle)");
    auto parents = particles(*p.model(), to(p));
    for (const auto& P : parents)
        if (std::dynamic_pointer_cast<DecayingParticle>(P)
            and contains(dynamic_cast<DecayingParticle*>(P.get())))
            return true;
    return false;
}

//-------------------------
const bool from::contains(const DecayingState* const p) const
{        
    return std::find(objects().begin(), objects().end(), p) != objects().end();
}

    
    
}
