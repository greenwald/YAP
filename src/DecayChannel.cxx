#include "DecayChannel.h"

#include "container_utils.h"
#include "Attributes.h"
#include "DecayingState.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "Model.h"
#include "ParticleCombination.h"
#include "ParticleCombinationCache.h"
#include "Spin.h"
#include "SpinAmplitude.h"
#include "SpinAmplitudeCache.h"

#include <algorithm>
#include <iterator>

namespace yap {

//-------------------------
DecayChannel::DecayChannel(const ParticleVector& daughters) :
    Daughters_(daughters)
{
    // check daughter size
    if (Daughters_.empty())
        throw exceptions::Exception("No daughters", "DecayChannel::DecayChannel");
    if (Daughters_.size() == 1)
        throw exceptions::Exception("Only one daughter", "DecayChannel::DecayChannel");

    // check no Daughters_ are empty
    if (std::any_of(Daughters_.begin(), Daughters_.end(), std::logical_not<ParticleVector::value_type>()))
        throw exceptions::Exception("Empty daughter", "DecayChannel::DecayChannel");

    // if more than two daughters, check all are spin 0
    if (Daughters_.size() > 2 and
        std::any_of(Daughters_.begin(), Daughters_.end(), [](const ParticleVector::value_type & d) {return d->quantumNumbers().twoJ() != 0;}))
        throw exceptions::Exception("Attempted to create nonresonant decay with spinfull daughters---use DecayingParticle without mass shape instead.", "DecayChannel::DecayChannel");

    // check that (first daughter's) Model is not nullptr
    if (model() == nullptr)
        throw exceptions::Exception(std::string("Model unset in ") + to_string(*Daughters_[0]),
                                    "DecayChannel::DecayChannel");

    // check that all daughters have same Model (trivially checks first daughter against itself)
    if (std::any_of(Daughters_.begin(), Daughters_.end(), [&](const ParticleVector::value_type & d) {return model() != d->model();}))
        throw exceptions::Exception("Model mismatch", "DecayChannel::DecayChannel");

    //////////////////////////////////////////////////
    // create ParticleCombination's for parent
    //

    // create vector of unique entries of daughters
    auto u_daughters = Daughters_;
    std::sort(u_daughters.begin(), u_daughters.end());
    u_daughters.erase(std::unique(u_daughters.begin(), u_daughters.end()), u_daughters.end());

    // collect combinations of ParticleCombination's of daughters
    // map: daughter -> vector of combinations of particle combinations
    std::map<ParticleVector::value_type, std::vector<ParticleCombinationVector> > d_CV_map;
    std::vector<ParticleCombinationVector> PCs;
    for (auto d : u_daughters) {
        // get daughter's list of particle combinations
        const auto& all_d_pcs = d->particleCombinations();

        // create set to copy them into
        ParticleCombinationVector d_pcs;
        
        // copy them if they aren't nullptr, aren't empty,
        // and aren't already contained in v, using equalDown to compare regardless of parent.
        std::copy_if(all_d_pcs.begin(), all_d_pcs.end(), std::back_inserter(d_pcs),
                     [&](const ParticleCombinationVector::value_type & pc) {
                         return pc and !pc->indices().empty() and
                             std::none_of(d_pcs.begin(), d_pcs.end(), std::bind(equal_down, std::placeholders::_1, pc));
                     });
        
        if (d_pcs.empty())
            throw exceptions::Exception("No ParticleCombinations for daughter " + to_string(*d)
                                        + " in DecayChannel " + to_string(*this),
                                        "DecayChannel::DecayChannel");

        // create combinations of ParticleCombinations to accomodate occurance of daughter in final state
        auto cv = combinations(d_pcs, std::count(Daughters_.begin(), Daughters_.end(), d));
        if (cv.empty())
            throw exceptions::Exception("Could not form combinations of ParticleCombinations for daughter " + to_string(*d),
                                        "DecayChannel::DecayChannel");

        // place combinations into map:
        d_CV_map.emplace(d, cv);
    }

    // create vector of iterators to combinations for each unique daughter
    // initialized to begin()
    std::vector<std::vector<ParticleCombinationVector>::const_iterator> Its;
    Its.reserve(u_daughters.size());
    std::transform(u_daughters.begin(), u_daughters.end(), std::back_inserter(Its),
    [&](const ParticleVector::value_type & ud) {return d_CV_map[ud].begin();});

    // use 'odometer" style increment of iterators:
    // repeat until last iterator hits its end()
    while (Its.back() != d_CV_map[u_daughters.back()].end()) {

        // fill vector of ParticleCombination's to use to create parent ParticleCombination from
        ParticleCombinationVector pcv(Daughters_.size());
        // loop through iterators to combinations
        for  (size_t d = 0; d < Its.size(); ++d) {
            auto prev_it = Daughters_.begin();
            // loop over elements in particular combination
            for (const auto& pc : *Its[d]) {
                // find daughter to set value for, starting from one beyond previously found daughter
                auto it = std::find(prev_it, Daughters_.end(), u_daughters[d]);
                if (it == Daughters_.end())
                    throw exceptions::Exception("failed to create particle combinations", "DecayChannel::DecayChannel");
                // set daughter's ParticleCombination
                pcv[it - Daughters_.begin()] = pc;
                // increment prev_it to one beyond last found daughter
                prev_it = it + 1;
            }
        }

        // check vector for possible overlaps
        if (disjoint(pcv))
            // create new ParticleCombination for parent;
            // ParticleCombinationCache::composite copies all elts of
            // pcv setting the parents of all copies to the newly
            // created ParticleCombination
            addParticleCombination(*const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->particleCombinationCache().composite(pcv));

        // increment the "odometer":
        // increment first iterator
        ++Its[0];
        // loop through iterators: if current is at end(), reset it and increment next one
        for (size_t i = 0; (i < Its.size() - 1) and Its[i] == d_CV_map[u_daughters[i]].end(); ++i) {
            Its[i] = d_CV_map[u_daughters[i]].begin();
            ++Its[i + 1];
        }
        // the "odometer" has ticked over completely if the last iterator is now at its end()
    }

    // check that at least one combination was created above
    if (particleCombinations().empty())
        throw exceptions::Exception("ParticleCombinations_ is empty", "DecayChannel::DecayChannel");

    // if more than two daughters, add nonresonant spin amplitude
    if (Daughters_.size() > 2)
        // add SpinAmplitude retrieved from cache
        addSpinAmplitude(const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->spinAmplitudeCache()->spinAmplitude(0, spins(Daughters_), 0, 0));
}

//-------------------------
void DecayChannel::addParticleCombination(const ParticleCombination& pc)
{
    // check number of daughters in pc
    if (pc.daughters().size() != Daughters_.size())
        throw exceptions::Exception("ParticleCombination has wrong number of daughters ("
                                    + std::to_string(pc.daughters().size()) + " != "
                                    + std::to_string(Daughters_.size()) + ")",
                                    "DecayChannel::addParticleCombination");

    ParticleCombinations_.insert(pc.shared_from_this());

    // add pc's daughters to daughter particles;
    // pc's daughters have their parents set correctly.
    for (size_t i = 0; i < pc.daughters().size(); ++i)
        Daughters_[i]->addParticleCombination(*pc.daughters()[i]);

    // add to SpinAmplitude's
    for (auto& sa : SpinAmplitudes_)
        sa->addParticleCombination(pc);

}

//-------------------------
void DecayChannel::registerWithModel()
{
    for (auto& d : Daughters_)
        d->registerWithModel();

    for (auto& sa : SpinAmplitudes_)
        sa->registerWithModel();
}

//-------------------------
void DecayChannel::pruneParticleCombinations()
{
    prune_particle_combinations(ParticleCombinations_);

    for (auto& p : Daughters_)
        p->pruneParticleCombinations();
}

//-------------------------
void DecayChannel::addAllPossibleSpinAmplitudes(unsigned two_J, int p)
{
    auto two_j = spins(Daughters_);
    p *= parity(Daughters_);

    // create spin amplitudes
    // loop over possible S: |j1-j2| <= S <= (j1+j2)
    for (unsigned two_S = std::abs<int>(two_j[0] - two_j[1]); two_S <= two_j[0] + two_j[1]; two_S += 2)
        // loop over possible L: |J-s| <= L <= (J+s)
        for (unsigned L = std::abs<int>(two_J - two_S) / 2; L <= (two_J + two_S) / 2; ++L)
            // check parity conservation (also fulfilled if parity = 0)
            if (p * pow_negative_one(L) >= 0)
                // add SpinAmplitude retrieved from cache
                addSpinAmplitude(const_cast<Model*>(static_cast<const DecayChannel*>(this)->model())->spinAmplitudeCache()->spinAmplitude(two_J, two_j, L, two_S));
}

//-------------------------
const Model* DecayChannel::model() const
{
    return Daughters_[0]->model();
}

//-------------------------
void DecayChannel::addSpinAmplitude(std::shared_ptr<SpinAmplitude> sa)
{
    // if spin amplitude already held, do nothing
    if (std::find(SpinAmplitudes_.begin(), SpinAmplitudes_.end(), sa) != SpinAmplitudes_.end())
        return;

    // check number of daughters
    if (sa->finalTwoJ().size() != Daughters_.size())
        throw exceptions::Exception("Number of daughters doesn't match", "DecayChannel::addSpinAmplitude");

    // check against daughter quantum numbers
    for (size_t i = 0; i < Daughters_.size(); ++i)
        if (Daughters_[i]->quantumNumbers().twoJ() != sa->finalTwoJ()[i])
            throw exceptions::Exception("Spins don't match daughter's", "DecayChannel::addSpinAmplitude");

    if (!SpinAmplitudes_.empty() and SpinAmplitudes_[0]->initialTwoJ() != sa->initialTwoJ())
        throw exceptions::Exception("Spins don't match previously added", "DecayChannel::addSpinAmplitude");

    // add to SpinAmplitudes_
    SpinAmplitudes_.push_back(sa);

    // add this' ParticleCombination's to it
    for (auto& pc : particleCombinations())
        SpinAmplitudes_.back()->addParticleCombination(*pc);
}

//-------------------------
bool DecayChannel::consistent() const
{
    bool C = true;

    // check no daughters is empty
    if (std::any_of(Daughters_.begin(), Daughters_.end(), [](std::shared_ptr<Particle> d) {return !d;})) {
        FLOG(ERROR) << "null pointer in daughters vector.";
        C &= false;
    }
    // check daughters
    std::for_each(Daughters_.begin(), Daughters_.end(), [&](std::shared_ptr<Particle> d) {if (d) C &= d->consistent();});

    // loop over SpinAmplitude's
    for (const auto& sa : SpinAmplitudes_) {
        // check SpinAmplitude
        if (!sa) {
            FLOG(ERROR) << "A SpinAmplitude is empty";
            C &= false;
        } else
            C &= sa->consistent();
    }

    return C;
}

//-------------------------
ParticleSet particles(const DecayChannel& dc)
{
    ParticleSet S;
    for (const auto& d : dc.daughters()) {
        if (is_final_state_particle(*d))
            S.insert(d);
        else if (is_decaying_state(*d)) {
            auto s = particles(static_cast<const DecayingState&>(*d));
            S.insert(s.begin(), s.end());
        } else
            throw exceptions::Exception("neither final-state particle nor decaying state", "find_particles");
    }
    return S;
}

//-------------------------
std::string to_string(const DecayChannel& dc)
{
    return dc.daughters().empty() ? "[nothing]"
        : particle(*dc.model(), has_decay_channel(&dc))->name() + " --> " + to_string(dc.daughters());
}

//-------------------------
std::string to_string(const DecayChannel& dc, int two_M, const SpinProjectionVector& two_m)
{
    return dc.daughters().empty() ? "[nothing]"
        : particle(*dc.model(), has_decay_channel(&dc))->name() + " [m = " + spin_to_string(two_M)
        + "] --> " + to_string(dc.daughters(), two_m);
}

}
