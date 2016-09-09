/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/// \file

#ifndef yap__Filters_h
#define yap__Filters_h

#include "fwd/Filters.h"

#include "fwd/DecayChannel.h"
#include "fwd/DecayingParticle.h"
#include "fwd/DecayTree.h"
#include "fwd/FinalStateParticle.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/MassShape.h"
#include "fwd/MassShapeWithNominalMass.h"
#include "fwd/Particle.h"
#include "fwd/Resonance.h"
#include "fwd/SpinAmplitude.h"

#include "Exceptions.h"

#include <algorithm>
#include <memory>
#include <set>

namespace yap {

/// \defgroup filters
/// @{

/// Throws if container's size is not 1
/// \return lone element in container
template <typename container>
typename container::value_type lone_elt(container& C)
{
    if (C.size() != 1)
        throw yap::exceptions::Exception("Container size not 1 (" + std::to_string(C.size()) + ")", "lone_elt");
    return *C.begin();
}

/// Throws if container's size is not 1
/// \return lone element in container
template <typename container>
typename container::value_type lone_elt(container&& C)
{
    if (C.size() != 1)
        throw yap::exceptions::Exception("Container size not 1 (" + std::to_string(C.size()) + ")", "lone_elt");
    return *C.begin();
}

/// filter through only the members of a container of shared_ptr's
/// that evaluate to true with all of the predicates given
/// \todo generalize beyond std::set
template <typename T, typename Last, typename ... UnaryPredicates>
std::set<std::shared_ptr<T> > filter(const std::set<std::shared_ptr<T> >& S, Last p, UnaryPredicates ... P)
{
    auto s = filter<T, UnaryPredicates...>(S, P...);
    for (auto it = s.begin(); it != s.end(); ) {
        if (p(**it))
            ++it;
        else
            it = s.erase(it);
    }
    return s;
}

// does nothing, needed to terminate above variadic template
template <typename T>
const std::set<std::shared_ptr<T> >& filter(const std::set<std::shared_ptr<T> >& S)
{ return S; }

/// \class filter_decay_tree
/// \brief helper base class for functor classes that can filter DecayTree's
struct filter_decay_tree {
    /// DecayTree functor
    virtual const bool operator()(const DecayTree& dt) const = 0;
};

/// \class filter_free_amplitude
/// \brief helper base class for functor classes that can filter FreeAmplitude's
struct filter_free_amplitude : public filter_decay_tree {
    using filter_decay_tree::operator();

    /// DecayTree functor
    virtual const bool operator()(const DecayTree& dt) const override;

    /// FreeAmplitude functor
    virtual const bool operator()(const FreeAmplitude& fa) const = 0;
};

/// \class filter_decay_channel
/// \brief helper base class for functor classes that can filter DecayChannel's
struct filter_decay_channel : public filter_free_amplitude {
    using filter_free_amplitude::operator();

    /// FreeAmplitude functor
    virtual const bool operator()(const FreeAmplitude& fa) const override;

    /// DecayChannel functor
    virtual const bool operator()(const DecayChannel& dc) const = 0;
};

/// \class filter_spin_amplitude
/// \brief helper base class for functor classes that can filter SpinAmplitude's
struct filter_spin_amplitude : public filter_free_amplitude {
    using filter_free_amplitude::operator();

    /// FreeAmplitude functor
    virtual const bool operator()(const FreeAmplitude& fa) const override;

    /// SpinAmplitude functor
    virtual const bool operator()(const SpinAmplitude& dc) const = 0;
};

/// \class filter_particle
/// \brief helper base class for filtering Particle's
struct filter_particle {
    /// Particle functor
    virtual const bool operator()(const Particle& p) const = 0;
};

/// \class filter_decaying_particle
/// \brief helper base class for filtering DecayingParticle's
struct filter_decaying_particle : public filter_particle {
    using filter_particle::operator();

    /// Particle functor
    virtual const bool operator()(const Particle& p) const override;

    /// DecayingParticle functor
    virtual const bool operator()(const DecayingParticle& dp) const = 0;
};

/// \class filter_resonance
/// \brief helper base class for filtering Resonance's
struct filter_resonance : public filter_decaying_particle {
    using filter_decaying_particle::operator();

    /// DecayingParticle functor
    virtual const bool operator()(const DecayingParticle& dp) const override;

    /// Resonance functor
    virtual const bool operator()(const Resonance& dp) const = 0;
};

/// \class filter_mass_shape
/// \brief helper base class for filtering MassShape's
struct filter_mass_shape : public filter_resonance {
    using filter_resonance::operator();

    /// Resonance functor
    virtual const bool operator()(const Resonance& dp) const override;

    /// MassShape functor
    virtual const bool operator()(const MassShape& dp) const = 0;

};

/// calls functor with shared_ptr to object;
/// \param f functor to call
/// \param ptr shared_ptr to object to call on
template <class F, typename T>
const bool by_ptr(const F& f, const std::shared_ptr<T>& ptr)
{ return ptr and f(*ptr); }

/// calls functor with raw pointer to object
/// \param f functor to call
/// \param ptr raw pointer to object to call on
template <class F, typename T>
const bool by_ptr(const F& f, const T* ptr)
{ return ptr and f(*ptr); }

/// \class has_decay_channel
/// \brief Functor class for filtering by decay channel
class has_decay_channel : public filter_decay_channel, public filter_decaying_particle
{
public:
    /// constructor;
    /// takes a pointer so that checking against nullptr is possible
    has_decay_channel(const DecayChannel* dc)
        : filter_decay_channel(), filter_decaying_particle(), DecayChannel_(dc) {}

    using filter_decay_channel::operator();
    using filter_decaying_particle::operator();

    /// FreeAmplitude functor
    const bool operator()(const FreeAmplitude& fa) const override;

    /// DecayChannel functor
    const bool operator()(const DecayChannel& dc) const override
    { return operator()(&dc); }

    /// shared_ptr<DecayChannel> functor
    const bool operator()(const std::shared_ptr<DecayChannel> dc) const
    { return operator()(dc.get()); }

    /// DecayChannel* functor
    const bool operator()(const DecayChannel* dc) const
    { return dc == DecayChannel_; }

    /// DecayingParticle functor
    const bool operator()(const DecayingParticle& dp) const override;

private:
    /// DecayChannel pointer to check equality to
    const DecayChannel* DecayChannel_;
};

/// \class has_spin_amplitude
/// \brief Functor class for filtering by spin amplitude
class has_spin_amplitude : public filter_spin_amplitude
{
public:
    /// Constructor;
    /// takes a pointer so that checking against nullptr is possible
    has_spin_amplitude(const SpinAmplitude* sa)
        : filter_spin_amplitude(), SpinAmplitude_(sa) {}

    using filter_spin_amplitude::operator();

    /// FreeAmplitude functor
    const bool operator()(const FreeAmplitude& fa) const override;

    /// SpinAmplitude functor
    const bool operator()(const SpinAmplitude& sa) const override
    { return operator()(&sa); }

    /// shared_ptr<SpinAmplitude> functor
    const bool operator()(const std::shared_ptr<SpinAmplitude>& sa) const
    { return operator()(sa.get()); }

    /// SpinAmplitude* functor
    const bool operator()(const SpinAmplitude* sa) const
    { return sa == SpinAmplitude_; }

private:
    /// SpinAmplitude pointer to check equality to
    const SpinAmplitude* SpinAmplitude_;
};

/// \class to
/// \brief Functor class for filtering by particle content
class to : public filter_decay_channel, public filter_decaying_particle
{
public:

    /// constructor
    /// \param daughters Daughter particles to check for
    to(const ParticleVector& daughters)
        : filter_decay_channel(), filter_decaying_particle(), Daughters_(daughters)
    {
        if (std::any_of(Daughters_.begin(), Daughters_.end(), std::logical_not<ParticleVector::value_type>()))
            throw exceptions::Exception("nullptr daughter provided", "to::to");
    }

    /// constructor (variadic)
    template <typename ... Others>
        to(std::shared_ptr<Particle> A, Others ... others) : to(ParticleVector{A, others...}) {}

    using filter_decay_channel::operator();
    using filter_decaying_particle::operator();

    /// DecayChannel functor
    virtual const bool operator()(const DecayChannel& dc) const override;

    /// DecayingParticle functor
    virtual const bool operator()(const DecayingParticle& dp) const override;

protected:

    /// Particle content to check equality to
    ParticleVector Daughters_;

};

/// \class exactly_to
/// \brief Functor class for filtering by exact particle content
class exactly_to : public to
{
public:
    /// constructor
    /// \param daughters Daughter particles to check for
    exactly_to(const ParticleVector& daughters) : to(daughters) {}

    /// constructor (variadic)
    template <typename ... Others>
    exactly_to(const ParticleVector::value_type& A, Others ... others)
        : to(A, others...) {}

    using to::operator();

    /// DecayChannel functor
    virtual const bool operator()(const DecayChannel& dc) const override;
};

/// Functor object for filtering by final state
class final_state : public filter_decay_channel, public filter_decaying_particle
{
    /// constructor
    /// \param daughters Daughter particles to check for
    final_state(const FinalStateParticleVector& fspv)
        : filter_decay_channel(), filter_decaying_particle(), FSP_(fspv)
    {
        if (std::any_of(FSP_.begin(), FSP_.end(), std::logical_not<FinalStateParticleVector::value_type>()))
            throw exceptions::Exception("nullptr daughter provided", "final_state::final_state");
    }

    /// constructor (variadic)
    template <typename ... Others>
    to(std::shared_ptr<FinalStateParticle> A, Others ... others) : final_state(FinalStateParticleVector{A, others...}) {}

    using filter_decay_channel::operator();
    using filter_decaying_particle::operator();

    /// DecayChannel functor
    virtual const bool operator()(const DecayChannel& dc) const override;

    /// DecayingParticle functor
    virtual const bool operator()(const DecayingParticle& dp) const override;

protected:

    /// Particle content to check equality to
    FinalStateParticleVector FSP_;

};

/// \class l_equals
/// \brief Functor object for filtering by orbital angular momentum
class l_equals : public filter_spin_amplitude
{
public:

    /// constructor
    l_equals(unsigned l) : filter_spin_amplitude(), L_(l) {}

    using filter_spin_amplitude::operator();

    /// SpinAmplitude functor
    virtual const bool operator()(const SpinAmplitude& sa) const override;

private:

    /// orbital angular momentum to check equality to
    unsigned L_;

};

/// \class m_equals
/// \brief Functor object for filtering by spin projection
class m_equals : filter_free_amplitude
{
public:

    /// constructor
    m_equals(int two_m) : filter_free_amplitude(), TwoM_(two_m) {}

    using filter_free_amplitude::operator();

    /// FreeAmplitude functor
    virtual const bool operator()(const FreeAmplitude& fa) const override;

private:

    /// spin projection to check equality to
    int TwoM_;

};

/// \class is_named
/// \brief Functor class for filtering objects that have a #name() function
class is_named
{
public:
    /// constructor
    is_named(std::string name) : Name_(name) {}

    /// functor
    template <typename T>
    const bool operator()(const T& t) const
    { return t.name() == Name_; }

private:
    /// name to check equality to
    std::string Name_;
};

/// \class has_decay_tree
/// \brief Functor class for checking if DecayingParticle has specified DecayTree
class has_decay_tree : public filter_decaying_particle, public filter_decay_tree
{
public:
    /// constructor
    has_decay_tree(const DecayTree* dt)
        : filter_decaying_particle(), filter_decay_tree(), DecayTree_(dt) {}

    using filter_decaying_particle::operator();
    using filter_decay_tree::operator();

    /// DecayTree functor
    const bool operator()(const DecayTree& dt) const
    { return operator()(&dt); }

    /// shared_ptr<DecayTree> functor
    const bool operator()(const std::shared_ptr<DecayTree>& dt) const
    { return operator()(dt.get()); }

    /// DecayTree* functor
    const bool operator()(const DecayTree* dt) const
    { return dt == DecayTree_; }

    /// DecayingParticle functor
    const bool operator()(const DecayingParticle& dp) const;

private:
    /// DecayTree pointer to search for
    const DecayTree* DecayTree_;
};

/// \class has_free_amplitude
/// \brief Functor class for checking whether objects have a particular free amplitude
class has_free_amplitude : public filter_decaying_particle, public filter_free_amplitude
{
public:
    /// constructor
    has_free_amplitude(const FreeAmplitude* fa)
        : filter_decaying_particle(), filter_free_amplitude(), FreeAmplitude_(fa) {}

    using filter_decaying_particle::operator();
    using filter_free_amplitude::operator();

    /// DecayTree functor
    const bool operator()(const DecayTree& dt) const override;

    /// DecayingParticle functor
    const bool operator()(const DecayingParticle& dp) const override;

    /// FreeAmplitude functor
    const bool operator()(const FreeAmplitude& fa) const override
    { return operator()(&fa); }

    /// shared_ptr<FreeAmplitude> functor
    const bool operator()(const std::shared_ptr<FreeAmplitude>& fa) const
    { return operator()(fa.get()); }

    /// FreeAmplitude* functor
    const bool operator()(const FreeAmplitude* fa) const
    { return fa == FreeAmplitude_; }

private:
    /// FreeAmplitude to look for
    const FreeAmplitude* FreeAmplitude_;
};

/// \class from
/// \brief Functor object for filtering by parent particle
class from : public filter_decay_channel, public filter_particle
{
public:
    /// Constructor
    from(const DecayingParticle& dp)
        : filter_decay_channel(), filter_particle(), DecayingParticle_(&dp) {}

    /// Constructor
    from(const std::shared_ptr<DecayingParticle>& dp)
        : filter_decay_channel(), filter_particle(), DecayingParticle_(dp.get()) {}

    using filter_decay_channel::operator();
    using filter_particle::operator();

    /// DecayTree functor
    const bool operator()(const DecayTree& dt) const override
    { return has_decay_tree(&dt)(*DecayingParticle_); }

    /// FreeAmplitude functor
    const bool operator()(const FreeAmplitude& fa) const override
    { return has_free_amplitude(&fa)(*DecayingParticle_); }

    /// DecayChannel functor
    const bool operator()(const DecayChannel& dc) const override
    { return has_decay_channel(&dc)(*DecayingParticle_); }

    /// Particle functor
    const bool operator()(const Particle& p) const override;

private:
    /// DecayingParticle to check for
    const DecayingParticle* DecayingParticle_;
};


/// \class has_mass
/// \brief Functors return true if particle is a final state particle
/// or if a mass shape can be found that inherits from
/// MassShapeWithNominalMass
class has_mass : public filter_mass_shape
{
public:
    /// constructor
    has_mass() : filter_mass_shape() {}

    using filter_mass_shape::operator();

    /// Particle functor
    virtual const bool operator()(const Particle& p) const override;

    /// FinalStateParticle functor
    const bool operator()(const FinalStateParticle& p) const
    { return true; }

    /// MassShape functor
    virtual const bool operator()(const MassShape& m) const override;

    /// MassShapeWithNominalMass functor
    const bool operator()(const MassShapeWithNominalMass& m) const
    { return true; }

};

/// @}
// end of defgroup filters

}

#endif
