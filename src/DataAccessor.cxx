#include "DataAccessor.h"

#include "CachedValue.h"
#include "FourMomenta.h"
#include "HelicityAngles.h"
#include "logging.h"
#include "MeasuredBreakupMomenta.h"
#include "Model.h"
#include "RequiresHelicityAngles.h"
#include "RequiresMeasuredBreakupMomenta.h"

namespace yap {

//-------------------------
DataAccessor::DataAccessor(const ParticleCombinationEqualTo& equal) :
    Equal_(equal),
    NIndices_(0),
    Size_(0),
    Index_(-1)
{
}

//-------------------------
void DataAccessor::printParticleCombinations() const
{
    for (auto& kv : SymmetrizationIndices_)
        LOG(INFO) << kv.second << " : " << to_string_with_parent(*(kv.first));
}

//-------------------------
bool DataAccessor::consistent() const
{
    bool C = true;

    // check CachedValues_
    for (auto& c : CachedValues_)
        if (c->owner() != this) {
            FLOG(ERROR) << "CachedValue's owner != this";
            C &= false;
        }

    return C;
}

//-------------------------
void DataAccessor::addParticleCombination(const ParticleCombination& c)
{
    // if already indexed, do nothing
    if (SymmetrizationIndices_.find(c.shared_from_this()) != SymmetrizationIndices_.end())
        return;

    // search for match using Equal
    auto it = std::find_if(SymmetrizationIndices_.begin(), SymmetrizationIndices_.end(),
                           [&](const ParticleCombinationMap<unsigned>::value_type & kv)
                           {return Equal_(c.shared_from_this(), kv.first);});

    // if found, use found index
    if (it != SymmetrizationIndices_.end()) {
        SymmetrizationIndices_.emplace(c.shared_from_this(), it->second);
        return;
    }

    // else assign to one higher than current highest index
    SymmetrizationIndices_.emplace(c.shared_from_this(), static_cast<unsigned>(NIndices_));
    // and increase current highest index
    ++NIndices_;
}

//-------------------------
void DataAccessor::pruneSymmetrizationIndices()
{
    if (!model())
        throw exceptions::Exception("Model not set", "DataAccessor::pruneSymmetrizationIndices");

    // remove entries that don't trace back to an ISP
    for (auto it = SymmetrizationIndices_.begin(); it != SymmetrizationIndices_.end(); ) {
        if (is_from_initial_state_particle_combination(*it->first, *model()))
            ++it;
        else
            it = SymmetrizationIndices_.erase(it);
    }

    //
    // fix indices now for holes
    //

    // collect used indices
    std::set<unsigned> used;
    for (const auto& kv : SymmetrizationIndices_)
        used.insert(kv.second);

    // repair
    unsigned index = 0;
    while (index < used.size()) {

        // if index is not used
        if (used.find(index) == used.end()) {
            // clear used
            used.clear();
            // reduce all indices greater than index by 1
            // and rebuild used
            for (auto& kv : SymmetrizationIndices_) {
                if (kv.second > index)
                    kv.second -= 1;
                used.insert(kv.second);
            }
        }

        //if index is (now) used, increment by 1
        if (used.find(index) != used.end())
            index += 1;

    }

    // reset NIndices_
    NIndices_ = 0;
    for (const auto& kv : SymmetrizationIndices_)
        NIndices_ = std::max(kv.second + 1, NIndices_);
}

//-------------------------
void DataAccessor::registerWithModel()
{
    if (!model())
        throw exceptions::Exception("Model unset", "DataAccessor::registerWithModel");

    if (model()->locked())
        throw exceptions::Exception("Model is locked and cannot be modified", "DataAccessor::registerWithModel");

    for (auto pc_i : symmetrizationIndices())
        const_cast<Model*>(model())->fourMomenta()->addParticleCombination(*pc_i.first);

    // if HelicityAngles is required
    if (dynamic_cast<RequiresHelicityAngles*>(this) and dynamic_cast<RequiresHelicityAngles&>(*this).requiresHelicityAngles()) {
        const_cast<Model*>(model())->requireHelicityAngles();
        for (auto pc_i : symmetrizationIndices())
            const_cast<Model*>(model())->helicityAngles()->addParticleCombination(*pc_i.first);
    }

    // if MeasuredBreakupMomenta is required
    if (dynamic_cast<RequiresMeasuredBreakupMomenta*>(this) and dynamic_cast<RequiresMeasuredBreakupMomenta*>(this)->requiresMeasuredBreakupMomenta()) {
        const_cast<Model*>(model())->requireMeasuredBreakupMomenta();
        for (auto pc_i : symmetrizationIndices())
            const_cast<Model*>(model())->measuredBreakupMomenta()->addParticleCombination(*pc_i.first);
    }

    // if stores nothing, do nothing
    if (size() == 0)
        return;

    // insert this into Model's DataAccessors_
    const_cast<Model*>(model())->DataAccessors_.insert(this);
}

//-------------------------
void DataAccessor::addCachedValue(std::shared_ptr<CachedValue> c)
{
    if (not c)
        throw exceptions::Exception("CachedValue is NULL", "DataAccessor::addCachedValue");

    // add CachedValue
    if (CachedValues_.insert(c).second) {
        // if insertion was successful

        // set its index
        c->setIndex(CachedValues_.size() - 1);

        // set its position
        c->setPosition(size());

        // increase data size to accommodate CachedValue
        increaseSize(c->size());
    }
}

//-------------------------
void remove_expired(DataAccessorSet& S)
{
    for (auto it = S.begin(); it != S.end(); )
        if (!*it) it = S.erase(it);
        else ++it;
}


}

