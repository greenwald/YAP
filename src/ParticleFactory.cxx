#include "ParticleFactory.h"

#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "MassShape.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <utility>

namespace yap {

//-------------------------
ParticleTableEntry::ParticleTableEntry(int pdg, const std::string& name, QuantumNumbers q, double mass, std::vector<double> parameters) :
    PDG_(pdg),
    Name_(name),
    QuantumNumbers_(q),
    Mass_(mass),
    MassShapeParameters_(parameters)
{
    if (Name_.empty())
        throw exceptions::Exception("Name is empty", "ParticleTableEntry::ParticleTableEntry");
    if (Mass_ < 0)
        throw exceptions::Exception("Mass is negative", "ParticleTableEntry::ParticleTableEntry");
}

//-------------------------
const double get_nth_entry(const ParticleTableEntry& pte, size_t n, const std::string& where)
{
    if (pte.massShapeParameters().size() <= n)
        throw exceptions::Exception("pte.massShapeParameters() does not contain an " + std::to_string(n) + "'th element", where);
    return pte.massShapeParameters()[n];
}

//-------------------------
ParticleFactory& ParticleFactory::operator+=(const ParticleFactory& rhs)
{
    std::transform(rhs.ParticleTable_.begin(), rhs.ParticleTable_.end(), inserter(*this),
                   std::bind(&ParticleTableMap::value_type::second, std::placeholders::_1));
    return *this;
}

//-------------------------
const ParticleTableEntry& ParticleFactory::operator[](int PDG) const
{
    if (ParticleTable_.count(PDG) == 0)
        throw exceptions::Exception("No particle table entry for PDG " + std::to_string(PDG), "ParticleFactory::operator[]");
    return ParticleTable_.at(PDG);
}

//-------------------------
const ParticleTableEntry& ParticleFactory::operator[](const std::string& name) const
{
    auto it = std::find_if(ParticleTable_.begin(), ParticleTable_.end(),
                           [&](const std::map<int, ParticleTableEntry>::value_type & p) {return p.second.name() == name;});
    if (it == ParticleTable_.end())
        throw exceptions::Exception("particle with name \"" + name + "\" not found", "ParticleFactory::pdg");
    return it->second;
}


//-------------------------
std::pair<ParticleTableMap::iterator, bool> ParticleFactory::insert(const ParticleTableEntry& entry)
{
    auto it_b = ParticleTable_.insert(ParticleTableMap::value_type(entry.pdg(), entry));

    // if insertion failed because key value entry.PDG was already contained
    if (!it_b.second and it_b.first != ParticleTable_.end()) {
        LOG(WARNING) << "PDG code " << entry.pdg() << " already exists. Overwriting entry.";
        it_b.first->second = entry;
        it_b.second = true;
    }

    return it_b;
}

//-------------------------
ParticleTableMap::iterator ParticleFactory::insert(ParticleTableMap::iterator hint, const ParticleTableEntry& entry)
{
    if (ParticleTable_.count(entry.pdg()) == 0)
        return ParticleTable_.insert(hint, ParticleTableMap::value_type(entry.pdg(), entry));

    LOG(WARNING) << "PDG code " << entry.pdg() << " already exists. Overwriting entry.";
    auto it = ParticleTable_.find(entry.pdg());
    it->second = entry;
    return it;
}

}
