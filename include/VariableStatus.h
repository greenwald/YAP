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

#ifndef yap_VariableStatus_h
#define yap_VariableStatus_h

namespace yap {

/// \enum VariableStatus
enum VariableStatus {
    kChanged   = -1,        ///< Parameter is free and has been changed
    kFixed     = 0,         ///< Parameter is fixed
    kUnchanged = +1,        ///< Parameter is free but has not been changed
};

inline std::ostream& operator<<(std::ostream& str, VariableStatus s)
{
    switch (s) {
        case kChanged:
            return str << "kChanged";
        case kFixed:
            return str << "kFixed";
        case kUnchanged:
            return str << "kUnchanged";
        default:
            return str << (int) s;
    }
}

}

#endif
