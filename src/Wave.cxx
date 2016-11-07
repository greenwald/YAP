#include "Wave.h"

#include "DecayChannel.h"
#include "Model.h"
#include "Spin.h"

#include <cstdlib>

namespace yap {

//-------------------------
Wave::Wave(const std::string& name, const QuantumNumbers& q, unsigned l, unsigned two_s,
           double radial_size, const ParticleVector& daughters) :
    DecayingState(name, q, radial_size)
{
    if (daughters.size() != 2)
        throw exceptions::Exception("Only two-body waves allowed", "Wave::Wave");

    if (!triangle(quantumNumbers().twoJ(), 2 * l, two_s))
        throw exceptions::Exception("JLS triangle violated", "Wave::Wave");

    auto two_j = spins(daughters);
    if (!triangle(two_s, two_j[0], two_j[1]))
        throw exceptions::Exception("j1j2S triangle violated", "Wave::Wave");

    // create DecayChannel
    auto dc = std::make_shared<DecayChannel>(daughters);
    if (!dc->model())
        throw exceptions::Exception("DecayChannel's model is nullptr", "Wave::Wave");

    // create spin amplitude
    dc->addSpinAmplitude(const_cast<Model*>(dc->model())->spinAmplitudeCache()->spinAmplitude(quantumNumbers().twoJ(), two_j, l, two_s));

    addDecayChannel(dc);
}

//-------------------------
std::shared_ptr<Wave> Wave::create(const std::string& name, const QuantumNumbers& q, unsigned l, double radial_size, const ParticleVector& daughters)
{
    if (daughters.size() != 2)
        throw exceptions::Exception("Only two-body waves allowed", "Wave::create");

    auto two_j = spins(daughters);
    
    // if only one possible s from daughter spins
    if (std::abs(static_cast<int>(two_j[0]) - static_cast<int>(two_j[1])) == two_j[0] + two_j[1])
        return create(name, q, l, two_j[0] + two_j[1], radial_size, daughters);

    // if only one possible s from J and L
    if (std::abs(static_cast<int>(q.twoJ()) - 2 * static_cast<int>(l)) == q.twoJ() + 2 * l)
        return create(name, q, l, q.twoJ() + 2 * l, radial_size, daughters);

    // else throw
    throw exceptions::Exception("S is ambiguous", "Wave::create");
}

//-------------------------
std::shared_ptr<Wave> Wave::create(const std::string& name, const QuantumNumbers& q, double radial_size, const ParticleVector& daughters)
{
    if (daughters.size() != 2)
        throw exceptions::Exception("Only two-body waves allowed", "Wave::create");

    auto two_j = spins(daughters);
    
    // check for unambiguous s
    if (std::abs(static_cast<int>(two_j[0]) - static_cast<int>(two_j[1])) != two_j[0] + two_j[1])
        throw exceptions::Exception("S is ambiguous", "Wave::create");

    auto two_s = two_j[0] + two_j[1];
    
    // check for unambiguous l
    if (std::abs(static_cast<int>(q.twoJ()) - static_cast<int>(two_s)) != q.twoJ() + two_s)
        throw exceptions::Exception("L is ambiguous", "Wave::create");

    // else return
    return create(name, q, q.twoJ() + two_s, two_s, radial_size, daughters);
}

//-------------------------
void Wave::addDecayChannel(std::shared_ptr<DecayChannel> dc)
{
    if (!channels().empty())
        throw exceptions::Exception("Wave already contains DecayChannel. Only one allowed.", "Wave::addDecayChannel");

    if (dc->spinAmplitudes().size() != 1)
        throw exceptions::Exception("DecayChannel has more than one SpinAmplitude.", "Wave::addDecayChannel");

    DecayingState::addDecayChannel(dc);
}

}
