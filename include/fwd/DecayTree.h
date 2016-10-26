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

#ifndef yap_DecayTreeFwd_h
#define yap_DecayTreeFwd_h

#include "fwd/FreeAmplitude.h"

#include <map>
#include <memory>
#include <set>
#include <vector>

namespace yap {

class DecayTree;

/// \typedef DecayTreeVector
using DecayTreeVector = std::vector<std::shared_ptr<DecayTree> >;

/// \typedef DecayTreeSet
using DecayTreeSet = std::set<std::shared_ptr<DecayTree> >;

/// \typedef DecayTreeVectorMap
/// Map of spin projection (int) to DecayTreeVector
/* using DecayTreeVectorMap = std::map<int, DecayTreeVector>; */

/// \return the highest free amplitude of a decay tree
std::shared_ptr<FreeAmplitude> free_amplitude(const DecayTree& dt);

}

#endif
