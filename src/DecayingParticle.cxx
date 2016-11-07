#include "DecayingParticle.h"

#include "DecayChannel.h"
#include "Model.h"
#include "Particle.h"

#include <memory>

namespace yap {

//-------------------------
const is_of_type<DecayingParticle> is_decaying_particle{};

//-------------------------
void DecayingParticle::addAllPossibleSpinAmplitudes(DecayChannel& dc, bool conserve_parity) const
{
    auto two_j = spins(dc.daughters());
    auto p = (conserve_parity) ? quantumNumbers().P() * parity(dc.daughters()) : 0;

    // create spin amplitudes
    // loop over possible S: |j1-j2| <= S <= (j1+j2)
    for (unsigned two_S = std::abs<int>(two_j[0] - two_j[1]); two_S <= two_j[0] + two_j[1]; two_S += 2)
        // loop over possible L: |J-s| <= L <= (J+s)
        for (unsigned L = std::abs<int>(quantumNumbers().twoJ() - two_S) / 2; L <= (quantumNumbers().twoJ() + two_S) / 2; ++L)
            // check parity conservation (also fulfilled if parity = 0)
            if (p * pow_negative_one(L) >= 0)
                // add SpinAmplitude retrieved from cache
                dc.addSpinAmplitude(const_cast<Model*>(dc.model())->spinAmplitudeCache()->spinAmplitude(quantumNumbers().twoJ(), two_j, L, two_S));
}

//-------------------------
std::shared_ptr<DecayChannel> DecayingParticle::addDecay(const ParticleVector& daughters, bool conserve_parity)
{
    auto dc = std::make_shared<DecayChannel>(daughters);
    addAllPossibleSpinAmplitudes(*dc, conserve_parity);
    addDecayChannel(dc);
    return dc;
}

}
