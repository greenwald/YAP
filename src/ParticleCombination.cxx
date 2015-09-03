#include "ParticleCombination.h"

#include "logging.h"

#include <algorithm>
#include <assert.h>

namespace yap {

//-------------------------
ParticleCombination::ParticleCombination() :
    Parent_(nullptr)
{
}

//-------------------------
ParticleCombination::ParticleCombination(ParticleIndex index) :
    Parent_(nullptr),
    Indices_(1, index)
{
}

//-------------------------
ParticleCombination::ParticleCombination(std::vector<std::shared_ptr<ParticleCombination> > c) :
    Parent_(nullptr)
{
    for (std::shared_ptr<ParticleCombination>& d : c)
        addDaughter(d);
}

//-------------------------
bool ParticleCombination::addDaughter(std::shared_ptr<ParticleCombination> daughter)
{
    if (daughter->indices().empty()) {
        LOG(ERROR) << "ParticleCombination::addDaughter - daughter contains no indices.";
        return false;
    }

    /// \todo Check that new daughter does not share content with other daughters?



    /*if (daughter->Parent_ == nullptr) {
        // daughter has no parent yet. Set this as parent and add
        daughter->Parent_ = this;
    } else if (daughter->Parent_ == this) {
        // fine
    } else {
        LOG(ERROR) << "this should not happen!";

        // daughter has already different parent -> make copy, set this as parent and get unique shared_ptr
        std::shared_ptr<ParticleCombination> copy(new ParticleCombination(*daughter));
        copy->Parent_ = this;
        std::shared_ptr<ParticleCombination> uniqueCopy = uniqueSharedPtr(copy);
        // need to swap so that parent of argument daughter is now set

        // \todo does not work since it would have to operate on the original shared_ptr object
        daughter.swap(uniqueCopy);
    }

    assert(daughter->Parent_ == this);*/

    // add daughter to vector
    Daughters_.push_back(daughter);

    // copy daughter's indices into Indices_
    Indices_.insert(Indices_.end(), daughter->indices().begin(), daughter->indices().end());

    return true;
}

//-------------------------
bool ParticleCombination::consistent() const
{
    // should have no daughters or 2 or more daughters:
    if (Daughters_.size() == 1) {
        LOG(ERROR) << "ParticleCombination::consistent() - has only one daughter.";
        return false;
    }

    // if empty, return inconsistent
    if (Indices_.empty()) {
        LOG(ERROR) << "ParticleCombination::consistent() - has no indices.";
        return false;
    }

    // if Daugthers_ empty, should have one and only index
    if (Daughters_.empty()) {
        if (Indices_.size() != 1) {
            LOG(ERROR) << "ParticleCombination::consistent() - contains wrong number of indices for final-state particle (" << Indices_.size() << " != 1)";
            return false;
        } else
            // don't need to check indices below
            return true;
    }

    // Check indices & and then daughters
    bool result = true;

    // create unique_copy of Indices_ (as set)
    std::set<ParticleIndex> U(Indices_.begin(), Indices_.end());
    // check unique_copy-object's size == this object's size
    if (U.size() != Indices_.size()) {
        LOG(ERROR) << "ParticleCombination::consistent - index vector contains duplicate entries (" << U.size() << " != " << Indices_.size() << ").";
        result = false;
    }

    // check daughters
    /*for (std::shared_ptr<ParticleCombination> d : Daughters_) {
        if (d->parent() != this) {
            LOG(ERROR) << "ParticleCombination::consistent - daughter's parent is not this ParticleCombination.";
            result = false;
        }
        result &= d->consistent();
    }*/

    return result;
}

//-------------------------
ParticleCombination::operator std::string() const
{
    std::string result = "(";
    for (ParticleIndex i : Indices_)
        result += std::to_string(static_cast<unsigned>(i));
    result += ")";

    if (Daughters_.empty() || (Daughters_.size() == 2 and Indices_.size() == 2))
        return result;

    result += " -> ";

    for (auto& d : Daughters_) {
        if (d->daughters().empty() || (d->daughters().size() == 2 and d->indices().size() == 2))
            result += std::string(*d);
        else
            result += "[" + std::string(*d) + "]";

    }


    return result;
}

//-------------------------
bool ParticleCombination::sharesIndices(std::shared_ptr<ParticleCombination> B)
{
    for (ParticleIndex a : Indices_)
        for (ParticleIndex b : B->indices())
            if (a == b)
                return true;

    return false;
}

//-------------------------
bool operator==(const ParticleCombination& A, const ParticleCombination& B)
{
    // Check indices
    if (A.Indices_.size() != B.Indices_.size())
        return false;
    if (!std::equal(A.Indices_.begin(), A.Indices_.end(), B.Indices_.begin()))
        return false;

    // check parents
    if (A.parent() != B.parent())
        return false;

    // Check daughters
    if (A.Daughters_.size() != B.Daughters_.size())
        return false;
    for (unsigned i = 0; i < A.Daughters_.size(); ++i)
        if (*(A.Daughters_[i]) != *(B.Daughters_[i]))
            return false;

    // a match!
    return true;
}

/////////////////////////
// Static stuff:

std::set<std::shared_ptr<ParticleCombination> > ParticleCombination::ParticleCombinationSet_;

//-------------------------
std::shared_ptr<ParticleCombination> ParticleCombination::uniqueSharedPtr(std::shared_ptr<ParticleCombination> pc)
{
    for (std::shared_ptr<ParticleCombination> d : ParticleCombinationSet_)
        if ((*pc) == (*d))
            return d;
    ParticleCombinationSet_.insert(pc);
    return pc;
}

//-------------------------
std::shared_ptr<ParticleCombination> ParticleCombination::uniqueSharedPtr(ParticleIndex i)
{
    return uniqueSharedPtr(std::make_shared<ParticleCombination>(i));
}

//-------------------------
std::shared_ptr<ParticleCombination> ParticleCombination::uniqueSharedPtr(std::vector<std::shared_ptr<ParticleCombination> > c)
{
    return uniqueSharedPtr(std::make_shared<yap::ParticleCombination>(c));
}

//-------------------------
void ParticleCombination::printParticleCombinationSet()
{
  std::cout << "ParticleCombination set:\n";
  for (auto pc : ParticleCombinationSet_) {
    std::cout << "  " << std::string(*pc) << "\n";
  }
  std::cout << std::endl;
}

}
