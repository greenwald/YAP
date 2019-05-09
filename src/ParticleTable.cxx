#include "ParticleTable.h"

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
ParticleTableEntry cp_conjugate(ParticleTableEntry t, const std::string& prefix)
{
    t.quantumNumbers().setQ(-t.quantumNumbers().Q());
    auto cp_name = t.name();
    switch (cp_name.back()) {
    case '+' :
        cp_name.back() = '-';
        break;
    case '-' :
        cp_name.back() = '+';
        break;
    default:
        cp_name = prefix + cp_name;
    }
    t.setName(cp_name);
    t.pdg(-t.pdg());
    return t;
}

//-------------------------
ParticleTableEntry cp_conjugate(const ParticleTableMap::value_type& t, const std::string& prefix)
{
    return cp_conjugate(std::get<1>(t), prefix);
}

//-------------------------
double get_nth_element(const ParticleTableEntry& pde, size_t n, const std::string& where)
{
    if (n >= pde.massShapeParameters().size())
        throw exceptions::Exception("MassShapeParameters_ has no " + std::to_string(n) + "'th element", where);
    return pde.massShapeParameters()[n];
}

//-------------------------
ParticleTable& ParticleTable::operator+=(const ParticleTable& rhs)
{
    std::transform(rhs.ParticleTableMap_.begin(), rhs.ParticleTableMap_.end(), inserter(*this),
                   std::bind(&ParticleTableMap::value_type::second, std::placeholders::_1));
    return *this;
}

//-------------------------
const ParticleTableEntry& ParticleTable::operator[](int PDG) const
{
    if (ParticleTableMap_.count(PDG) == 0)
        throw exceptions::Exception("No particle table entry for PDG " + std::to_string(PDG), "ParticleTable::operator[]");
    return ParticleTableMap_.at(PDG);
}

//-------------------------
const ParticleTableEntry& ParticleTable::operator[](const std::string& name) const
{
    auto it = std::find_if(ParticleTableMap_.begin(), ParticleTableMap_.end(),
    [&](const std::map<int, ParticleTableEntry>::value_type & p) {return p.second.name() == name;});
    if (it == ParticleTableMap_.end())
        throw exceptions::Exception("particle with name \"" + name + "\" not found", "ParticleTable::operator[]");
    return it->second;
}

//-------------------------
std::pair<ParticleTableMap::iterator, bool> ParticleTable::insert(const ParticleTableEntry& entry)
{
    auto it_b = ParticleTableMap_.insert(ParticleTableMap::value_type(entry.pdg(), entry));

    // if insertion failed because key value entry.PDG was already contained
    if (!it_b.second and it_b.first != ParticleTableMap_.end()) {
        LOG(WARNING) << "PDG code " << entry.pdg() << " already exists. Overwriting entry.";
        it_b.first->second = entry;
        it_b.second = true;
    }

    return it_b;
}

//-------------------------
ParticleTableMap::iterator ParticleTable::insert(ParticleTableMap::iterator hint, const ParticleTableEntry& entry)
{
    if (ParticleTableMap_.count(entry.pdg()) == 0)
        return ParticleTableMap_.insert(hint, ParticleTableMap::value_type(entry.pdg(), entry));

    LOG(WARNING) << "PDG code " << entry.pdg() << " already exists. Overwriting entry.";
    auto it = ParticleTableMap_.find(entry.pdg());
    it->second = entry;
    return it;
}

}
