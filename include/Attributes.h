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

#ifndef yap__Attributes_h
#define yap__Attributes_h

#include "fwd/Attributes.h"

#include "fwd/BlattWeisskopf.h"
#include "fwd/DecayChannel.h"
#include "fwd/DecayTree.h"
#include "fwd/DecayingState.h"
#include "fwd/FreeAmplitude.h"
#include "fwd/MassShape.h"
#include "fwd/Parameter.h"
#include "fwd/Particle.h"
#include "fwd/SpinAmplitude.h"

#include "AttributeUtilities.h"
#include "Exceptions.h"
#include "Particle.h"

#include <functional>
#include <memory>
#include <vector>

namespace yap {

/// Functor class to return orbital angular momentum of argument
struct orbital_angular_momentum : public attribute_of<const unsigned, SpinAmplitude, FreeAmplitude, BlattWeisskopf, DecayTree>
{
    /// \note functors inherited
    using attribute_of::operator();

    /// SpinAmplitude& functor
    virtual const unsigned operator()(const SpinAmplitude& sa) const override;

    /// FreeAmplitude& functor
    virtual const unsigned operator()(const FreeAmplitude& fa) const override;

    /// BlattWeisskopf functor
    virtual const unsigned operator()(const BlattWeisskopf& bw) const override;

    /// DecayTree functor
    /// checks top level FreeAmplitude of DecayTree
    virtual const unsigned operator()(const DecayTree& dt) const override;
 };

/// Functor class to return (twice the) spin angular momentum of argument
struct spin_angular_momentum : public attribute_of<const unsigned, SpinAmplitude, FreeAmplitude, DecayTree, Particle>
{
    /// \note functors inherited
    using attribute_of::operator();

    /// SpinAmplitude& functor
    virtual const unsigned operator()(const SpinAmplitude& sa) const override;

    /// FreeAmplitude& functor
    virtual const unsigned operator()(const FreeAmplitude& fa) const override;

    /// DecayTree functor
    /// checks top level FreeAmplitude of DecayTree
    virtual const unsigned operator()(const DecayTree& dt) const override;

    /// BlattWeisskopf functor
    virtual const unsigned operator()(const Particle& bw) const override;
};

/// Functor class to return (twice the) spin angular momentum of argument
struct spin_projection : public attribute_of<const int, DecayTree>
{
    /// \note functors inherited
    using attribute_of::operator();

    /// DecayTree& functor
    virtual const int operator()(const DecayTree& dt) const override;
};

/// Functor class to check whether argument is fixed
/// \ingroup Attributes
struct is_fixed : public attribute_of<const bool, ParameterBase, DecayTree>
{
    /// \note functors inherited
    using attribute_of::operator();

    /// ParameterBase& functor
    virtual const bool operator()(const ParameterBase& p) const override;

    /// DecayTree& functor
    virtual const bool operator()(const DecayTree& dt) const override;
};

/// Functor class to check whether argument is not fixed
/// \ingroup Attributes
struct is_not_fixed : attribute_of<const bool, ParameterBase, DecayTree>
{
    /// \note functors inherited
    using attribute_of::operator();

    /// ParameterBase& functor
    virtual const bool operator()(const ParameterBase& p) const override;

    /// DecayTree& functor
    virtual const bool operator()(const DecayTree& dt) const override;
};

/// An is_not_fixed object
extern const is_not_fixed is_free;

/// Functor class to check if state is decayed to by argument
/// \ingroup Attributes
struct to : public has_pointed_to_object<Particle, DecayChannel, FreeAmplitude, DecayTree, Particle>
{
    /// \note constructors inherited
    using has_pointed_to_object::has_pointed_to_object;

    /// \note functors inherited
    using has_pointed_to_object::operator();

    /// DecayChannel& functor
    virtual const bool operator()(const DecayChannel& dc) const override;
    
    /// FreeAmplitude& functor
    virtual const bool operator()(const FreeAmplitude& fa) const override;
    
    /// DecayTree& functor
    virtual const bool operator()(const DecayTree& dt) const override;

    /// Particle& functor
    /// returns true if true for any DecayChannel (if DecayingState)
    virtual const bool operator()(const Particle& p) const override;
};

/// Functor class to check if state is decayed exactly to by argument
/// \ingroup Attributes
struct exactly_to : public to
{
    /// \note constructors inherited from to
    using to::to;

    /// \note functors inherited
    using to::operator();

    /// DecayChannel& functor
    virtual const bool operator()(const DecayChannel& dc) const override;
};

/// Functor class to check whether argument has a particular FreeAmplitude
/// \ingroup Attributes
struct has_free_amplitude : public has_pointed_to_object<FreeAmplitude, DecayTree, Particle>
{
    /// \note constructors inherited
    using has_pointed_to_object::has_pointed_to_object;

    /// \note functors inherited
    using has_pointed_to_object::operator();

    /// DecayTree functor
    /// checks if DecayTree contains FreeAmplitude in FreeAmplitudes_
    virtual const bool operator()(const DecayTree& dt) const override;

    /// Particle functor
    /// checks all DecayTrees in (Decaying)Particle for a match (to the top-most FreeAmplitude of the tree)
    virtual const bool operator()(const Particle& p) const override;
};

/// Functor class to check whether argument has a particular DecayTree
/// \ingroup Attributes
struct has_decay_tree : public has_pointed_to_object<DecayTree, Particle, DecayTree>
{
    /// \note constructors inherited
    using has_pointed_to_object::has_pointed_to_object;

    /// \note functors inherited
    using has_pointed_to_object::operator();

    /// Particle& functor
    virtual const bool operator()(const Particle& p) const override;

    /// DecayTree& functor
    /// checks daughter decay trees of argument
    virtual const bool operator()(const DecayTree& dt) const override;
};

/// Functor class to check whehter argument has a particular DecayChannel
/// \ingroup Attributes
struct has_decay_channel : public has_pointed_to_object<DecayChannel, Particle, FreeAmplitude, DecayTree>
{
    /// \note constructors inherited
    using has_pointed_to_object::has_pointed_to_object;

    /// \note functors inherited
    using has_pointed_to_object::operator();

    /// Particle& functor
    virtual const bool operator()(const Particle& p) const override;

    /// FreeAmplitude& functor
    virtual const bool operator()(const FreeAmplitude& fa) const override;

    /// DecayTree& functor
    /// checks all FreeAmplitudes in DecayTree
    virtual const bool operator()(const DecayTree& dt) const override;
};

/// struct to access parent particle of an object
/// \ingroup Attributes
struct parent_state : public attribute_of<std::shared_ptr<const DecayingState>,
                                         DecayTree, FreeAmplitude, DecayChannel, BlattWeisskopf, MassShape>
{
    /// \note functors inherited
    using attribute_of::operator();

    /// DecayTree& functor
    virtual std::shared_ptr<const DecayingState> operator()(const DecayTree& dt) const override;

    /// FreeAmplitude& functor
    virtual std::shared_ptr<const DecayingState> operator()(const FreeAmplitude& fa) const override;

    /// DecayChannel& functor
    virtual std::shared_ptr<const DecayingState> operator()(const DecayChannel& dc) const override;

    /// BlattWeisskopf& functor
    virtual std::shared_ptr<const DecayingState> operator()(const BlattWeisskopf& bw) const override;

    /// MassShape& functor
    virtual std::shared_ptr<const DecayingState> operator()(const MassShape& m) const override;
};

/// functor to get name of return value of attribute
/// \ingroup Attributes
template <typename A>
class name_of : public with_return_type<std::string>
{
public:

    /// Constructor
    name_of() = default;
    
    /// Constructor
    name_of(const A& a) : Attr_(a) {}
    
    /// functor
    template <typename U>
    std::string operator()(const U& u) const
    { static auto get_name = std::mem_fn(&Particle::name); return get_name(Attr_(u)); }

private:
    /// default constructed attribute object
    A Attr_;
};

/// functor to return whether an argument has a MassShape
/// \ingroup Attributes
struct has_a_mass_shape : public attribute_of<const bool, Particle>
{
    /// \note functors inherited
    using attribute_of::operator();

    /// Particle& functor
    virtual const bool operator()(const Particle& p) const override;
};

/// functor to return whether an argument has a mass
/// \ingroup Attributes
struct has_a_mass : public attribute_of<const bool, MassShape, Particle>
{
    /// \note functor inherited
    using attribute_of::operator();
    
    /// MassShape& functor
    virtual const bool operator()(const MassShape& m) const override;
    
    /// Particle& functor
    virtual const bool operator()(const Particle& p) const override;
};

/// functor to return mass parameter
/// \ingroup Attributes
struct mass_parameter : public attribute_of<const RealParameter&, MassShape, Particle>
{
    /// \note functors inherited
    using attribute_of::operator();
    
    /// MassShape& functor
    virtual const RealParameter& operator()(const MassShape& m) const override;
    
    /// Particle& functor
    virtual const RealParameter& operator()(const Particle& p) const override;
};

/// functor to return whether object comes from a particular parent
/// \ingroup Attributes
struct from : public has_pointed_to_object<DecayingState, Particle>
{
    /// \note constructors inherited
    using has_pointed_to_object::has_pointed_to_object;
    
    /// \note functors inherited
    using has_pointed_to_object::operator();

    /// Particle& functor
    const bool operator()(const Particle& p) const;

    /// Generic functor
    template <typename U>
    const bool operator()(const U& u) const
    { static parent_state ParentOf_; return contains(ParentOf_(u).get()); }

    /// checks if 
    const bool contains(const DecayingState* const p) const;
};

}

#endif
