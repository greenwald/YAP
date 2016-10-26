#include "DecayingState.h"

#include "Attributes.h"
#include "BlattWeisskopf.h"
#include "container_utils.h"
#include "DecayChannel.h"
#include "DecayTree.h"
#include "Filter.h"
#include "FreeAmplitude.h"
#include "Parameter.h"
#include "SpinAmplitude.h"

#include "logging.h"

namespace yap {

//-------------------------
DecayingState::DecayingState(const std::string& name, const QuantumNumbers& q, double radial_size) :
    Particle(name, q),
    RadialSize_(std::make_shared<PositiveRealParameter>(radial_size))
{
}

//-------------------------
bool DecayingState::consistent() const
{
    bool C = Particle::consistent();

    if (DecayChannels_.empty()) {
        FLOG(ERROR) << "no channels specified.";
        return false; // further checks require at least one channel
    }

    // check no channel is empty
    if (std::any_of(DecayChannels_.begin(), DecayChannels_.end(), std::logical_not<DecayChannelVector::value_type>())) {
        FLOG(ERROR) << "DecayChannel vector contains nullptr";
        C &= false;
    }

    // check consistency of all channels
    std::for_each(DecayChannels_.begin(), DecayChannels_.end(), [&](const std::shared_ptr<DecayChannel>& dc) {if (dc) C &= dc->consistent();});

    return C;
}

//-------------------------
void DecayingState::checkDecayChannel(const DecayChannel& dc) const
{
    static charge Q;
    if (Q(dc) != Q(*this))
        throw exceptions::Exception("Incompatible charge (" + std::to_string(Q(dc)) + " != " + std::to_string(Q(*this)) + ")",
                                    "DecayingState::checkDecayChannel");
}

//-------------------------
std::shared_ptr<DecayChannel> DecayingState::addChannel(std::shared_ptr<DecayChannel> c)
{
    FLOG(INFO) << "in";
    if (!c)
        throw exceptions::Exception("DecayChannel empty", "DecayingState::addChannel");

    if (c->particleCombinations().empty())
        throw exceptions::Exception(std::string("DecayChannel has no ParticleCombinations - ") + to_string(*c),
                                    "DecayingState::addChannel");

    // check ISP
    if (!DecayChannels_.empty() and c->model() != model())
        throw exceptions::Exception("Model mismatch", "DecayingState::addChannel");

    // check if valid for DecayingState
    checkDecayChannel(*c);

    DecayChannels_.push_back(c);

    // if spin amplitudes haven't been added by hand, add all possible
    if (DecayChannels_.back()->spinAmplitudes().empty())
        DecayChannels_.back()->addAllPossibleSpinAmplitudes(quantumNumbers().twoJ());

    // create necessary BlattWeisskopf objects
    for (const auto& sa : DecayChannels_.back()->spinAmplitudes()) {
        // if BW is not already stored for L, add it
        if (BlattWeisskopfs_.find(sa->L()) == BlattWeisskopfs_.end())
            BlattWeisskopfs_.emplace(sa->L(), std::make_shared<BlattWeisskopf>(sa->L(), *this));
    }

    // add particle combinations
    for (auto pc : DecayChannels_.back()->particleCombinations())
        addParticleCombination(*pc);

    /// create decay trees for channel:

    /// loop over spin amplitudes of channel
    for (auto& sa : DecayChannels_.back()->spinAmplitudes()) {

        // create new FreeAmplitude
        auto fa = std::make_shared<FreeAmplitude>(DecayChannels_.back(), sa);

        // loop over possible parent spin projections
        for (const auto& two_M : sa->twoM()) {

            // loop over possible daughter spin projections
            for (const auto& two_m : sa->twoM(two_M)) {

                // loop over particle combinations of this decaying particle
                for (const auto& pc : particleCombinations()) {

                    // initialize vector of possible decay trees from
                    // initial DecayTree created above
                    DecayTreeVector DTV(1, std::make_shared<DecayTree>(fa, two_M, two_m));
                    modifyDecayTree(*DTV[0]);

                    // loop over daughters in channel
                    for (size_t d = 0; d < DecayChannels_.back()->daughters().size(); ++d) {

                        if (!DecayChannels_.back()->daughters()[d])
                            throw exceptions::Exception("daughter is nullptr", "DecayingState::addChannel");

                        // if daughter d is a DecayingState
                        if (is_decaying_state(*DecayChannels_.back()->daughters()[d])) {

                            // if daughter doesn't have tree for spin projection
                            if (std::none_of(DecayTrees_.begin(), DecayTrees_.end(), m_equals(two_m[d])))
                                DTV.clear();
                            else {
                                // create temp tree vector to store new copies into
                                DecayTreeVector DTV_temp;
                                
                                // loop over decay trees of daughter with appropriate spin projection and particle combination
                                for (const auto& dt : filter(std::dynamic_pointer_cast<DecayingState>(DecayChannels_.back()->daughters()[d])->decayTrees(),
                                                             m_equals(two_m[d]), has_particle_combination(pc->daughters()[d])))
                                    for (const auto& DT : DTV) {
                                        // add copy of DT to DTV_temp
                                        DTV_temp.push_back(std::make_shared<DecayTree>(*DT));
                                        // add decay tree to it
                                        DTV_temp.back()->setDaughterDecayTree(d, dt);
                                    }
                                // replace DTV with DTV_temp
                                DTV = DTV_temp;
                            }
                        }
                        // else if daughter d is a FinalStateParticle
                        else if (is_final_state_particle(*DecayChannels_.back()->daughters()[d])) {
                            // if it doesn't have the right ParticleCombination, clear DTV
                            if (!has_particle_combination(pc->daughters()[d])(DecayChannels_.back()->daughters()[d]))
                                DTV.clear();
                        } else
                            throw exceptions::Exception("Daughter is neither DecayingState nor FinalStateParticle", "DecayingState::addChannel");
                        
                        // if DTV is now empty, break
                        if (DTV.empty())
                            break;
                        
                    } // ends loop over daughters

                    // if decay trees were created, add them into DecayTrees_
                    // if they aren't already present in it
                    for (const auto& DT : DTV)
                        if (std::none_of(DecayTrees_.begin(), DecayTrees_.end(), [&DT](const DecayTreeVector::value_type& dt){return *dt == *DT;}))
                            DecayTrees_.push_back(DT);

                } // ends loop over particle combinations of this decaying particle
            } // ends loop over spin projections of daughters
        } // ends loop over spin projection of parent
    } // ends loop over spin amplitude

    return DecayChannels_.back();
}

//-------------------------
const Model* DecayingState::model() const
{
    return DecayChannels_.empty() ? nullptr : DecayChannels_[0]->model();
}

//-------------------------
void DecayingState::registerWithModel()
{
    for (auto& l_bw : BlattWeisskopfs_)
        if (l_bw.second)
            l_bw.second->registerWithModel();

    for (auto& c : DecayChannels_)
        c->registerWithModel();
}

//-------------------------
void DecayingState::addParticleCombination(const ParticleCombination& pc)
{
    Particle::addParticleCombination(pc);

    // add also to all BlattWeiskopf barrier factors
    if (pc.daughters().size() == 2)
        for (auto& l_bw : BlattWeisskopfs_)
            if (l_bw.second)
                l_bw.second->addParticleCombination(pc);

    // add to DecayChannels,
    // if DecayChannel contains particle combination with same content (without checking parent)
    // this is for the setting of ParticleCombination's with parents
    for (auto& dc : DecayChannels_) {
        if (std::any_of(dc->particleCombinations().begin(), dc->particleCombinations().end(), std::bind(&equal_down, pc.shared_from_this(), std::placeholders::_1)))
            dc->addParticleCombination(pc);
    }
}

//-------------------------
void DecayingState::pruneParticleCombinations()
{
    Particle::pruneParticleCombinations();

    for (auto& dc : DecayChannels_)
        dc->pruneParticleCombinations();
}

//-------------------------
void DecayingState::modifyDecayTree(DecayTree& dt) const
{
    if (!dt.freeAmplitude())
        throw exceptions::Exception("DecayTree has nullptr free amplitude", "DecayingState::modifyDecayTree");

    if (!dt.freeAmplitude()->spinAmplitude())
        throw exceptions::Exception("FreeAmplitude's SpinAmplitude is nullptr", "DecayingState::modifyDecayTree");

    // find BlattWeisskopf object
    if (dt.freeAmplitude()->spinAmplitude()->L() > 0) {
        auto bw = BlattWeisskopfs_.find(dt.freeAmplitude()->spinAmplitude()->L());
        if (bw == BlattWeisskopfs_.end())
            throw exceptions::Exception("No Blatt-Weisskopf factor found for L = "
                                        + std::to_string(dt.freeAmplitude()->spinAmplitude()->L()),
                                        "DecayingState::modifyDecayTree");
        
        if (!bw->second)
            throw exceptions::Exception("BlattWeisskopf is nullptr", "DecayingState::modifyDecayTree");

        // Add BlattWeisskopf object
        dt.addAmplitudeComponent(*bw->second);
    }
}

//-------------------------
ParticleSet particles(DecayingState& ds)
{
    ParticleSet S = {ds.shared_from_this()};
    for (const auto& dc : ds.channels()) {
        auto s = particles(*dc);
        S.insert(s.begin(), s.end());
    }
    return S;
}

//-------------------------
FreeAmplitudeSet free_amplitudes(const DecayingState& ds)
{
    return free_amplitudes(ds.decayTrees());
}

//-------------------------
// helper function
void get_paddings(const DecayingState& ds, size_t& name_padding, size_t& spinamp_padding)
{
    name_padding = std::max(name_padding, ds.name().length());

    // loop over channels
    for (auto& c : ds.channels()) {

        spinamp_padding = std::max(spinamp_padding, to_string(c->spinAmplitudes()).length());

        // loop over daughters of channel
        for (const auto& d : c->daughters()) {
            
            name_padding = std::max(name_padding, d->name().length());

            // recursively call
            if (d and is_decaying_state(*d))
                get_paddings(dynamic_cast<const DecayingState&>(*d), name_padding, spinamp_padding);

        }
    }
}

//-------------------------
std::string pad_right(const std::string& s, size_t len, char c = ' ')
{ return s + std::string(s.length() < len ? len - s.length() : 0, ' '); }

//-------------------------
std::string to_decay_string(const DecayingState& ds, unsigned level)
{
    // get maximum length of relevant strings
    static size_t name_padding = 0;
    static size_t spinamp_padding = 0;
    
    if (level == 0)
        get_paddings(ds, name_padding, spinamp_padding);

    std::string s;
    
    for (size_t i = 0; i < ds.channels().size(); ++i) {
        if (i > 0)
            s += "\n" + pad_right("", level * (3 * name_padding + 8 + spinamp_padding));

        s += pad_right(ds.name(), name_padding) + " ->";

        for (const auto& d : ds.channels()[i]->daughters())
            s += " " + pad_right(d->name(), name_padding);

        s += pad_right(to_string(ds.channels()[i]->spinAmplitudes()), spinamp_padding);
        
        for (const auto& d : ds.channels()[i]->daughters())
            if (d and is_decaying_state(*d))
                s += ", " + to_decay_string(static_cast<const DecayingState&>(*d), level + 1);
    }
    
    return s;
}

}
