#include "ParticleCombination.h"

#include "container_utils.h"
#include "DecayingState.h"
#include "Exceptions.h"
#include "logging.h"
#include "Model.h"

#include <algorithm>
#include <set>

namespace yap {

//-------------------------
void ParticleCombination::addDaughter(ParticleCombination& daughter)
{
    if (daughter.indices().empty())
        throw exceptions::Exception("daughter contains no indices", "ParticleCombination::addDaughter");

    /// Check that new daughter does not share content with other daughters?
    if (overlap(daughter.indices(), Indices_))
        throw exceptions::Exception("daughter overlaps with other daughters", "ParticleCombination::addDaughter");

    /// Check that new daughter doesn't already have parent set
    if (daughter.parent())
        throw exceptions::Exception("daughter's parent is already set", "ParticleCombination::addDaughter");

    // set daughter's parent to shared_from_this
    daughter.Parent_ = shared_from_this();

    // add daughter to vector
    Daughters_.push_back(daughter.shared_from_this());

    // copy daughter's indices into Indices_
    Indices_.insert(Indices_.end(), Daughters_.back()->indices().begin(), Daughters_.back()->indices().end());
}

//-------------------------
bool ParticleCombination::consistent() const
{
    bool C = true;

    // should have no daughters or 2 or more daughters:
    if (Daughters_.size() == 1) {
        FLOG(ERROR) << "has only one daughter.";
        C &= false;
    }

    if (Indices_.empty()) {
        FLOG(ERROR) << "has no indices.";
        C &= false;
    }

    // if has daughters
    if (!Daughters_.empty()) {
        // Check Indices_ doesn't have duplicates
        // create unique_copy of Indices_ (as set)
        std::set<unsigned> U(Indices_.begin(), Indices_.end());
        // check unique_copy-object's size == this object's size
        if (U.size() != Indices_.size()) {
            FLOG(ERROR) << "index vector contains duplicate entries (" << U.size() << " != " << Indices_.size() << ").";
            C &= false;
        }

        // count number of daughters with parent not set to this
        auto n = std::count_if(Daughters_.begin(), Daughters_.end(),
                               [&](const ParticleCombinationVector::value_type & d) {return d->parent().get() != this;});
        if (n != 0) {
            FLOG(ERROR) << n << " daughters' parent not set to this ParticleCombination.";
            C &= false;
        }

        // check consistency of daughters
        std::for_each(Daughters_.begin(), Daughters_.end(),
                      [&](const ParticleCombinationVector::value_type & d) {C &= d->consistent();});
    }
    // if Daugthers_ empty, should have one and only index (as FSP)
    else if (Indices_.size() != 1) {
        FLOG(ERROR) << "contains wrong number of indices for final-state particle (" << Indices_.size() << " != 1)";
        C &= false;
    }

    return C;
}

//-------------------------
bool disjoint(const ParticleCombinationVector& pcv)
{
    // check for any entry in pcv overlapping with another
    for (size_t i = 0; i < pcv.size(); ++i)
        for (size_t j = i + 1; j < pcv.size(); ++j)
            if (overlap(pcv[i]->indices(), pcv[j]->indices()))
                return false;
    return true;
}

//-------------------------
const bool is_initial_state_particle_combination(const ParticleCombination& pc, const Model& m)
{
    for (auto& isp_am : m.initialStates())
        if (isp_am.first->particleCombinations().find(std::const_pointer_cast<ParticleCombination>(pc.shared_from_this()))
                != isp_am.first->particleCombinations().end())
            return true;
    return false;
}

//-------------------------
void prune_particle_combinations(ParticleCombinationSet& PCs)
{
    // get ISP PCs
    size_t ispNIndices(0);
    for (auto& pc : PCs)
        ispNIndices = std::max(ispNIndices, origin(*pc).indices().size());

    for (auto it = PCs.begin(); it != PCs.end(); ) {
        if (origin(**it).indices().size() < ispNIndices)
            it = PCs.erase(it);
        else
            ++it;
    }
}

//-------------------------
std::string indices_string(const ParticleCombination& pc, std::string before, std::string after)
{
    if (pc.indices().empty())
        return "(empty)";
    std::string s = before;
    for (auto i : pc.indices())
        s += std::to_string(i);
    s += after;
    return s;
}

//-------------------------
std::string to_string(const ParticleCombination& pc)
{
    auto s = indices_string(pc);

    if (pc.daughters().empty())
        return s;

    s += " -> ";

    for (auto& d : pc.daughters()) {
        s += "(";
        for (auto i : d->indices())
            s += std::to_string(i);
        s += ") + ";
    }
    s.erase(s.size() - 3, 3);
    for (auto& d : pc.daughters())
        if (!is_final_state_particle_combination(*d))
            s += "; " + to_string(*d);
    return s;
}

//-------------------------
std::string to_string_with_parent(const ParticleCombination& pc)
{
    if (!pc.parent())
        return to_string(pc);

    return to_string(pc) + " in " + to_string(origin(pc));
}

//-------------------------
bool equal_by_ordered_content(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B)
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    // Check indices
    if (A->indices() != B->indices())
        return false;

    // a match!
    return true;
}

//-------------------------
bool equal_down(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B)
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    // check ordered content
    if (A->indices() != B->indices())
        return false;

    // Check daughters
    if (A->daughters().size() != B->daughters().size())
        return false;

    for (unsigned i = 0; i < A->daughters().size(); ++i)
        if (!equal_down(A->daughters()[i], B->daughters()[i]))
            return false;

    // a match!
    return true;
}

//-------------------------
bool equal_up(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B)
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!equal_by_ordered_content(A, B))
        return false;

    // if no more parents, return true
    if (!A->parent() and !B->parent())
        return true;

    // else continue up
    return equal_up(A->parent(), B->parent());
}

//-------------------------
bool equal_up_and_down(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B)
{
    if (!equal_down(A, B))
        return false;

    // if parents, return true
    if (!A->parent() and !B->parent())
        return true;

    // else check up
    return equal_up(A->parent(), B->parent());
}

//-------------------------
bool equal_by_orderless_content(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B)
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    // check size of indices vectors
    if (A->indices().size() != B->indices().size())
        return false;

    // check contents of indices vectors
    // (creating a set will sort entries for easy comparison,
    // since order doesn't matter)
    std::set<unsigned> a(A->indices().begin(), A->indices().end());
    std::set<unsigned> b(B->indices().begin(), B->indices().end());

    return std::equal(a.begin(), a.end(), b.begin());
}

//-------------------------
bool equal_down_by_orderless_content(const std::shared_ptr<const ParticleCombination>& A, const std::shared_ptr<const ParticleCombination>& B)
{
    // check if either empty
    if (!A or !B)
        return false;

    // compare shared_ptr addresses
    if (A == B)
        return true;

    if (!equal_by_orderless_content(A, B))
        return false;

    // Check daughters
    if (A->daughters().size() != B->daughters().size()) {
        return false;
    }

    unsigned matches(0);
    for (unsigned i = 0; i < A->daughters().size(); ++i)
        for (unsigned j = 0; j < B->daughters().size(); ++j)
            if (equal_by_orderless_content(A->daughters()[i], B->daughters()[j])) {
                ++matches;
                break;
            }

    if (matches == A->daughters().size())
        // a match!
        return true;

    return false;
}

}
