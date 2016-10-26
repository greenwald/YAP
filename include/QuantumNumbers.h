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

#ifndef yap_QuantumNumbers_h
#define yap_QuantumNumbers_h

#include "Exceptions.h"
#include "Spin.h"

#include <string>

namespace yap {

/// \class QuantumNumbers
/// \brief Quantum numbers of a Particle
/// \author Johannes Rauch, Daniel Greenwald
class QuantumNumbers
{
public:

    /// \name Constructors
    /// @{

    /// JQPICG constructor
    /// \todo restore constexpr if using c++14
    QuantumNumbers(unsigned twoJ, int Q, int P, unsigned twoI, int C, int G)
        : TwoJ_(twoJ), P_(signum(P)), C_(signum(C)), TwoI_(twoI), G_(signum(G)), Q_(Q)
    {
        if (Q_ != 0 and C_ != 0)
            throw exceptions::Exception("charged particle has nonzero charge parity", "QuantumNumbers::QuantumNumbers");
    }

    /// JQPI constructor
    constexpr QuantumNumbers(unsigned twoJ, int Q, int P = 0, unsigned twoI = 0)
        : TwoJ_(twoJ), P_(signum(P)), C_(0), TwoI_(twoI), G_(0), Q_(Q) {}

    /// @}

    /// \name Getters
    /// @{

    /// \return spin * 2
    constexpr unsigned twoJ() const
    { return TwoJ_; }

    /// \return parity
    constexpr int P() const
    { return P_; }

    /// \return C-parity
    constexpr int C() const
    { return C_; }

    /// \return Isospin * 2
    constexpr unsigned twoI() const
    { return TwoI_; }

    /// \return G-parity
    constexpr int G() const
    { return G_; }

    /// \return Electric charge
    constexpr int Q() const
    { return Q_; }

    /// @}

private:

    /// Spin * 2
    unsigned TwoJ_;

    /// Parity
    int P_;

    /// C-parity
    int C_;

    /// Isospin * 2
    unsigned TwoI_;

    /// G-parity
    int G_;

    /// Electric charge
    int Q_;

};

/// equality operator
inline const bool operator==(const QuantumNumbers& lhs, const QuantumNumbers& rhs)
{ return lhs.twoJ() == rhs.twoJ() and lhs.P() == rhs.P() and lhs.C() == rhs.C()
  and lhs.twoI() == rhs.twoI() and lhs.G() == rhs.G() and lhs.Q() == rhs.Q(); }

/// inequality operator
inline const bool operator!=(const QuantumNumbers& lhs, const QuantumNumbers& rhs)
{ return !(lhs == rhs); }

/// convert to string
inline std::string to_string(const QuantumNumbers& Q)
{ return spin_to_string(Q.twoJ()) + (Q.P() > 0 ? "+" : "-") + (Q.C() == 0 ? "" : (Q.C() > 0 ? "+" : "-")); }

}

#endif
