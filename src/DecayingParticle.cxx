#include "DecayingParticle.h"

#include "HelicitySpinAmplitude.h"
#include "FinalStateParticle.h"
#include "logging.h"
#include "QuantumNumbers.h"

#include <iomanip>
#include <memory>

namespace yap {

//-------------------------
DecayingParticle::DecayingParticle(InitialStateParticle* isp, const QuantumNumbers& q, double mass, std::string name, double radialSize) :
    AmplitudeComponent(),
    Particle(q, mass, name),
    AmplitudeComponentDataAccessor(isp),
    RadialSize_(radialSize)
{}

//-------------------------
Amp DecayingParticle::calcAmplitude(DataPoint& d, std::shared_ptr<ParticleCombination> pc)
{
    // \todo check
    Amp a = Complex_0;

    for (auto& c : channels()) {
        a += c->freeAmplitude() * c->amplitude(d, pc);
    }

    DEBUG("DecayingParticle: amplitude = " << a);

    return a;
}

//-------------------------
bool DecayingParticle::consistent() const
{
    bool result = true;

    result &= DataAccessor::consistent();
    result &= Particle::consistent();

    if (RadialSize_ <= 0.) {
        LOG(ERROR) << "DecayingParticle::consistent() - Radial size not positive.";
        result = false;
    }

    if (Channels_.empty()) {
        LOG(ERROR) << "DecayingParticle::consistent() - no channels specified.";
        return false; // further checks require at least one channel
    }

    for (auto& c : Channels_) {
        if (!c) {
            LOG(ERROR) << "DecayingParticle::consistent() - DecayChannels contains null pointer.";
            result = false;
        } else if (this != c->parent()) {
            LOG(ERROR) << "DecayingParticle::consistent() - DecayChannels does not point back to this DecayingParticle.";
            result = false; // channel consistency check requires correct pointer
        }
        result &= c->consistent();
    }

    // check if all channels lead to same final state particles
    std::vector<std::shared_ptr<FinalStateParticle> > fsps0 = finalStateParticles(0);
    for (unsigned i = 1; i < nChannels(); ++i)
        if (finalStateParticles(i) != fsps0) {
            LOG(ERROR) << "DecayingParticle::consistent() - final state of channel " << i << " does not match.";
            result = false;
        }

    return result;
}

//-------------------------
void DecayingParticle::addChannel(DecayChannel* c)
{
    Channels_.push_back(std::unique_ptr<yap::DecayChannel>(c));
    Channels_.back()->setParent(this);

    if (c->particleCombinations().empty())
        LOG(ERROR) << "DecayingParticle::addChannel(c) - c->particleCombinations().empty()";

    for (std::shared_ptr<ParticleCombination> pc : c->particleCombinations()) {
        addSymmetrizationIndex(ParticleCombination::uniqueSharedPtr(pc));
    }
}

//-------------------------
void DecayingParticle::addChannels(std::vector<std::shared_ptr<Particle> > A, std::vector<std::shared_ptr<Particle> > B, unsigned maxTwoL)
{
    // consistency check
    for (auto& vec : {A, B}) {
        std::string name = vec[0]->name();
        std::set<int> helicities;
        for (auto& part : vec) {
            if (part->name() != name) {
                LOG(ERROR) << "DecayingParticle::addChannels: particles in vector are not of the same type!";
                return;
            }
            int hel = part->quantumNumbers().twoHelicity();
            if (helicities.find(hel) != helicities.end()) {
                LOG(ERROR) << "DecayingParticle::addChannels: particles in vector don't have different helicities!";
                return;
            }
            helicities.insert(hel);
        }
    }

    // loop over possible l
    for (unsigned twoL = 0; twoL < maxTwoL; ++twoL) {
        if (!SpinAmplitude::angularMomentumConserved(quantumNumbers(), A[0]->quantumNumbers(), B[0]->quantumNumbers(), twoL))
            continue;

        // loop over helicity combinations
        for (auto& partA : A) {
            for (auto& partB : B) {
                double cg = HelicitySpinAmplitude::calculateClebschGordanCoefficient(quantumNumbers(), partA->quantumNumbers(), partB->quantumNumbers(), twoL);
                if (cg == 0.)
                    continue;

                addChannel(new DecayChannel(partA, partB,
                                            std::make_shared<HelicitySpinAmplitude>(initialStateParticle(), quantumNumbers(),
                                                    partA->quantumNumbers(), partB->quantumNumbers(),
                                                    twoL, cg)));

                //std::cout << "add channel " << std::string(*Channels_.back()) << "\n";
            }
        }
    }

}


//-------------------------
std::vector< std::shared_ptr<FinalStateParticle> > DecayingParticle::finalStateParticles(unsigned i) const
{
    if (!Channels_[i])
        return std::vector<std::shared_ptr<FinalStateParticle>>();
    return Channels_[i]->finalStateParticles();
}

//-------------------------
void DecayingParticle::setParentHelicities(int twoParentHelicity)
{
    std::cout << "DecayingParticle::setParentHelicities(" << twoParentHelicity << ") of " << name() << "\n";

    /*bool found(false);
    for (unsigned i = 0; i < Channels_.size(); ++i) {
      std::shared_ptr<HelicitySpinAmplitude> amp = std::dynamic_pointer_cast<HelicitySpinAmplitude>(channel(i)->spinAmplitude());
      if (amp->twoParentHelicity() == twoParentHelicity) {
          found = true;
          break;
      }
    }*/



    /*for (unsigned i = 0; i < Channels_.size(); ++i) {

          std::shared_ptr<HelicitySpinAmplitude> amp = std::dynamic_pointer_cast<HelicitySpinAmplitude>(channel(i)->spinAmplitude());

          std::array<int, 2> helicities({0,0});
          if (amp) {
              helicities = amp->helicities();

              if (not amp->twoParentHelicitySet()) {
                  amp->setTwoParentHelicity(twoParentHelicity);
                  std::cout << "  setTwoParentHelicity of " << std::string(*channel(i)) << "\n";
              }
              else if (amp->twoParentHelicity() == twoParentHelicity) {
                // fine
                std::cout << "  twoParentHelicity already set of " << std::string(*channel(i)) << "\n";
              }
              else {
                std::shared_ptr<HelicitySpinAmplitude> newAmp = std::make_shared<HelicitySpinAmplitude>(initialStateParticle(), quantumNumbers(),
                       channel(i)->daughters()[0]->quantumNumbers(), channel(i)->daughters()[1]->quantumNumbers(),
                       amp->twoL(), amp->clebschGordanCoefficient());
                newAmp->setTwoParentHelicity(twoParentHelicity);
                std::unique_ptr<yap::DecayChannel> newChannel(new DecayChannel(channel(i)->daughters()[0], channel(i)->daughters()[1], newAmp));
                Channels_[i].swap(newChannel);
                std::cout << "  setTwoParentHelicity of " << std::string(*newAmp) << " to " << twoParentHelicity << "/2 and add new channel.\n";
              }
          }
          else {
              std::cout << "  meh\n";
          }


          unsigned iDaugh(0);

          for (std::shared_ptr<Particle> d : channel(i)->daughters()) {
              if (std::dynamic_pointer_cast<DecayingParticle>(d))
                  std::static_pointer_cast<DecayingParticle>(d)->setParentHelicities(helicities[iDaugh++]);
          }
    }*/
}

//-------------------------
void DecayingParticle::optimizeSpinAmplitudeSharing()
{
    /// \todo Will not work if we have more than one initial state particle
    /// make ampSet member of initialStateParticle ???
    static std::set<std::shared_ptr<SpinAmplitude> > ampSet;
    static SharedSpinAmplitudeComparator comp;

    for (unsigned i = 0; i < nChannels(); ++i) {

        std::shared_ptr<SpinAmplitude> amp = channel(i)->spinAmplitude();

        bool found(false);
        for (auto sa : ampSet) {
            if (comp(sa, amp)) {
                LOG(INFO) << "Sharing spin amplitude " << std::string(*sa) << " and " << std::string(*amp) << " (" << amp.get() << ")";
                channel(i)->spinAmplitude() = sa;
                found = true;
                break;
            }
        }
        if (!found)
            ampSet.insert(amp);

        for (std::shared_ptr<Particle> d : channel(i)->daughters()) {
            if (std::dynamic_pointer_cast<DecayingParticle>(d))
                std::static_pointer_cast<DecayingParticle>(d)->optimizeSpinAmplitudeSharing();
        }
    }
}

//-------------------------
void DecayingParticle::printDecayChainLevel(int level) const
{
    // get maximum length of particle names
    static size_t padding = 0;
    static size_t paddingSpinAmp = 0;
    if (padding == 0 || level == -1) {
        padding = std::max(padding, name().length());
        for (auto& c : Channels_) {
            for (std::shared_ptr<Particle> d : c->daughters()) {
                padding = std::max(padding, d->name().length());
                if (std::dynamic_pointer_cast<DecayingParticle>(d))
                    std::static_pointer_cast<DecayingParticle>(d)->printDecayChainLevel(-1);
                paddingSpinAmp = std::max(paddingSpinAmp, std::string(*c->spinAmplitude()).length());
            }
        }
        if (level == -1)
            return;
    }

    for (unsigned int i = 0; i < nChannels(); ++i) {
        if (i > 0)
            std::cout << "\n" << std::setw(level * (padding * 3 + 7 + paddingSpinAmp)) << "";

        std::cout << std::left << std::setw(padding) << this->name() << " ->";
        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            std::cout << " " << std::setw(padding) << d->name();
        std::cout << std::left << std::setw(paddingSpinAmp)
                  << std::string(*(channel(i)->spinAmplitude()));

        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
                std::cout << ",  ";
                std::static_pointer_cast<DecayingParticle>(d)->printDecayChainLevel(level + 1);
            }
    }

    if (level == 0)
        std::cout << "\n";
}

//-------------------------
void DecayingParticle::printSpinAmplitudes(int level)
{
    static std::set<std::shared_ptr<SpinAmplitude> > ampSet;

    // get maximum length of particle names
    static size_t padding = 0;
    static size_t paddingSpinAmp = 0;
    if (padding == 0 || level == -1) {
        ampSet.clear();
        padding = std::max(padding, name().length());
        for (auto& c : Channels_) {
            for (std::shared_ptr<Particle> d : c->daughters()) {
                padding = std::max(padding, d->name().length());
                if (std::dynamic_pointer_cast<DecayingParticle>(d))
                    std::static_pointer_cast<DecayingParticle>(d)->printSpinAmplitudes(-1);
                paddingSpinAmp = std::max(paddingSpinAmp, std::string(*c->spinAmplitude()).length());
            }
        }
        if (level == -1)
            return;
    }

    for (unsigned int i = 0; i < nChannels(); ++i) {
        if (i > 0)
            std::cout << "\n" << std::setw(level * (padding * 3 + 7 + paddingSpinAmp)) << "";

        std::cout << std::left << std::setw(padding) << this->name() << " ->";
        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            std::cout << " " << std::setw(padding) << d->name();

        std::shared_ptr<SpinAmplitude> amp = channel(i)->spinAmplitude();
        if (ampSet.find(amp) == ampSet.end()) {
            ampSet.insert(amp);
            std::cout << std::left << std::setw(paddingSpinAmp) << std::string(*amp);
        } else {
            std::cout << std::left << std::setw(paddingSpinAmp - 1) << "sharedSpinAmp";
        }


        for (std::shared_ptr<Particle> d : channel(i)->daughters())
            if (std::dynamic_pointer_cast<DecayingParticle>(d)) {
                std::cout << ",  ";
                std::static_pointer_cast<DecayingParticle>(d)->printSpinAmplitudes(level + 1);
            }
    }

    if (level == 0)
        std::cout << "\n";
}

//-------------------------
void DecayingParticle::setSymmetrizationIndexParents()
{
    //std::cout << "DecayingParticle::setSymmetrizationIndexParents()\n";

    // clean up PCs without parents
    std::vector<std::shared_ptr<ParticleCombination> > PCsParents = particleCombinations();
    auto it = PCsParents.begin();
    while (it != PCsParents.end()) {
        if (!(*it)->parent()) {
            it = PCsParents.erase(it);
        } else
            ++it;
    }
    clearSymmetrizationIndices();

    for (auto& pc : PCsParents)
        addSymmetrizationIndex(pc);


    /*std::cout << "  Particle combinations in DecayingParticle " <<  name() << "\n";
    for (auto& chPC : particleCombinations()) {
      std::cout << "    " << std::string(*chPC);
      if (chPC->parent())
        std::cout << " from decay " << std::string(*chPC->parent());
      std::cout << "\n";
    }*/


    // next level
    for (auto& ch : channels())
        ch->setSymmetrizationIndexParents();


}


}
