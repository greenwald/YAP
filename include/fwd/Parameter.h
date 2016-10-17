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

#ifndef yap_ParameterFwd_h
#define yap_ParameterFwd_h

#include <memory>
#include <set>
#include <vector>

namespace yap {

class ParameterBase;
template <typename T> class Parameter;
class ComplexParameter;
class RealParameter;
class PositiveRealParameter;
class ComplexComponentParameter;
class RealComponentParameter;
class ImaginaryComponentParameter;

/// \typedef ParameterVector
/// \ingroup Parameters
/// Use when order of parameters must be kept constant
using ParameterVector = std::vector<std::shared_ptr<ParameterBase> >;

/// \typedef ParameterSet
/// \ingroup Parameters
/// Use when enforcement of no duplicates is required, and order does not matter
using ParameterSet = std::set<std::shared_ptr<ParameterBase> >;

/// \typedef ComplexParameterVector
/// \ingroup Parameters
using ComplexParameterVector = std::vector<std::shared_ptr<ComplexParameter> >;

/// \typedef RealParameterVector
/// \ingroup Parameters
using RealParameterVector = std::vector<std::shared_ptr<RealParameter> >;

}

#endif
