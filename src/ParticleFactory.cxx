#include "ParticleFactory.h"

#include "DecayingParticle.h"
#include "FinalStateParticle.h"
#include "InitialStateParticle.h"
#include "logging.h"
#include "make_unique.h"
#include "Resonance.h"
#include "SpinAmplitude.h"

#include <fstream>
#include <iostream>

namespace yap {

//-------------------------
ParticleTableEntry::ParticleTableEntry(int pdg, std::string name, QuantumNumbers q, double mass, std::vector<double> parameters) :
    QuantumNumbers(q),
    PDG_(pdg),
    Name_(name),
    Mass_(mass),
    MassShapeParameters_(parameters)
{
}

//-------------------------
bool ParticleTableEntry::consistent() const
{
    bool result = QuantumNumbers::consistent();

    if (Name_.empty()) {
        LOG(ERROR) << "ParticleTableEntry::consistent : No name specified.";
        result = false;
    }

    return result;
}

//-------------------------
ParticleFactory::ParticleFactory(const std::string pdlFile)
{
    readPDT(pdlFile);
}

//-------------------------
std::unique_ptr<FinalStateParticle> ParticleFactory::createFinalStateParticle(int PDG)
{
    const auto& p = particleTableEntry(PDG);
    DEBUG("make FinalStateParticle " << p.Name_ << " with quantum numbers " << p);
    return std::make_unique<FinalStateParticle>(p, p.Mass_, p.Name_);
}

//-------------------------
std::unique_ptr<InitialStateParticle> ParticleFactory::createInitialStateParticle(int PDG, double radialSize)
{
    const auto& p = particleTableEntry(PDG);

    if (p.twoJ() != 0)
        LOG(ERROR) << "InitialStateParticle has spin != 0. ";

    DEBUG("make InitialStateParticle " << p.Name_ << " with quantum numbers " << p);
    return std::make_unique<InitialStateParticle>(p, p.Mass_, p.Name_, radialSize);
}

//-------------------------
std::unique_ptr<Resonance> ParticleFactory::createResonance(int PDG, double radialSize, std::unique_ptr<MassShape>&& massShape)
{
    const auto& p = particleTableEntry(PDG);
    DEBUG("make Resonance " << p.Name_ << " with quantum numbers " << p);
    massShape->setParameters(p);
    return std::make_unique<Resonance>(p, p.Mass_, p.Name_, radialSize, std::move(massShape));
}

//-------------------------
const ParticleTableEntry& ParticleFactory::particleTableEntry(int PDG) const
{
    if (particleTable_.count(PDG) == 0) {
        LOG(FATAL) << "ParticleFactory::particleTableEntry : No particle table entry for PDG " << PDG;
    }
    return particleTable_.at(PDG);
}

//-------------------------
bool ParticleFactory::addParticleTableEntry(ParticleTableEntry entry)
{
    if (!entry.consistent()) {
        LOG(ERROR) << "ParticleFactory::addParticleTableEntry : entry with PDG code " << entry.PDG_ << " inconsistent";
        return false;
    }

    if (particleTable_.count(entry.PDG_) != 0) {
        LOG(WARNING) << "ParticleFactory::addParticleTableEntry : PDG code " << entry.PDG_ << " already exists. Overwriting entry.";
    }

    particleTable_[entry.PDG_] = entry;
    // } else
    //     particleTable_.insert(std::make_pair<int, ParticleTableEntry>(entry.PDG_, entry));

    return true;
}

//-------------------------
int ParticleFactory::pdgCode(std::string name) const
{
    auto it = std::find_if(particleTable_.begin(), particleTable_.end(),
    [&](const std::map<int, ParticleTableEntry>::value_type & p) {return p.second.Name_ == name;});
    if (it == particleTable_.end()) {
        LOG(ERROR) << "ParticleFactory::pdgCode - particle with name \"" << name << "\" not found.";
        return 0;
    }
    return it->first;
}

//-------------------------
void ParticleFactory::readPDT(const std::string fname)
{

    /**
     * This function was taken from EvtGen and modified
     *
     * // Copyright Information: See EvtGen/COPYRIGHT
     * //      Copyright (C) 1998      Caltech, UCSB
     *
     */

    std::ifstream indec;

    indec.open(fname.c_str());

    char cmnd[100];
    char xxxx[100];

    char pname[100];
    int  stdhepid;
    double mass;
    double pwidth;
    double pmaxwidth;
    int    chg3;
    int    spin2;
    double ctau;
    int    lundkc;
    //EvtId i;

    if (!indec) {
        LOG(ERROR) << "Could not open:" << fname.c_str() << "EvtPDL";
        return;
    }

    do {

        char ch, ch1;

        // ignoring commented lines
        do {
            indec.get(ch);
            if (ch == '\n') {
                indec.get(ch);
            }
            if (ch != '*') {
                indec.putback(ch);
            } else {
                while (indec.get(ch1), ch1 != '\n');
            }
        } while (ch == '*');

        indec >> cmnd;

        if (strcmp(cmnd, "end")) {
            if (!strcmp(cmnd, "add")) {
                indec >> xxxx;
                indec >> xxxx;
                indec >> pname;
                indec >> stdhepid;
                indec >> mass;
                indec >> pwidth;
                indec >> pmaxwidth;
                indec >> chg3;
                indec >> spin2;
                indec >> ctau;
                indec >> lundkc;

                // note: isospin & parity are missing from .pdl format
                addParticleTableEntry(ParticleTableEntry(stdhepid, pname,
                                      QuantumNumbers(0, spin2, 0, std::round(1. * chg3 / 3)),
                                      mass, {pwidth}));
            }

            // if find a set read information and discard it
            if (!strcmp(cmnd, "set")) {
                indec >> xxxx;
                indec >> xxxx;
                indec >> xxxx;
                indec >> xxxx;
            }
        }

    } while (strcmp(cmnd, "end"));
}

}
