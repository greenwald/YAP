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

#include "fwd/Spin.h"

#include <ostream>
#include <string>

namespace yap {

/// \class QuantumNumbers
/// \brief Quantum numbers of a Particle
/// \author Johannes Rauch

class QuantumNumbers
{
public:

    /// \name Constructors
    /// @{

    /// Constructor
    constexpr QuantumNumbers(unsigned twoJ, int P, int C, unsigned twoI, int G, int Q)
        : TwoJ_(twoJ), P_(P), C_(C), TwoI_(twoI), G_(G), Q_(Q) {}

    /// IJPQ constructor
    constexpr QuantumNumbers(unsigned twoI, unsigned twoJ, int P, int Q)
        : QuantumNumbers(twoJ, P, 0, twoI, 0, Q) {}

    /// JPQ constructor
    constexpr QuantumNumbers(unsigned twoJ, int P, int Q)
        : QuantumNumbers(twoJ, P, 0, 0, 0, Q) {}

    /// JQ constructor
    constexpr QuantumNumbers(unsigned twoJ, int Q)
        : QuantumNumbers(twoJ, 0, 0, 0, 0, Q) {}

    /// Default constructor
    /// is inconsistent
    QuantumNumbers()
        : QuantumNumbers(0, 0, 0, 0) {}

    /// @}

    /// check consistency
    virtual bool consistent() const;

    /// \name Getters
    /// @{

    /// \return spin * 2
    constexpr unsigned twoJ() const
    { return TwoJ_; }

    /// \return spin
    constexpr double J() const
    { return TwoJ_ * 0.5; }

    /// \return parity
    constexpr int P() const
    { return P_; }

    /// \return C-parity
    constexpr int C() const
    { return C_; }

    /// \return Isospin * 2
    constexpr unsigned twoI() const
    { return TwoI_; }

    /// \return Isospin
    constexpr double I() const
    { return TwoI_ * 0.5; }

    /// \return G-parity
    constexpr int G() const
    { return G_; }

    /// \return Electric intge
    constexpr int Q() const
    { return Q_; }

    /// @}

    /// \name Setters
    /// @{

    /// Set Spin
    void setJ(double J)
    { TwoJ_ = 2 * J; }

    /// Set 2 * Spin
    void setTwoJ(unsigned J)
    { TwoJ_ = J; }

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
bool operator==(const QuantumNumbers& lhs, const QuantumNumbers& rhs);

/// inequality operator
inline bool operator!=(const QuantumNumbers& lhs, const QuantumNumbers& rhs)
{ return !(lhs == rhs); }

/// convert to string
inline std::string to_string(const QuantumNumbers& Q)
{ return spin_to_string(Q.twoJ()) + (Q.P() > 0 ? "+" : "-") + (Q.C() == 0 ? "" : (Q.C() > 0 ? "+" : "-")); }

/// convert to string
inline std::string debug_string(const QuantumNumbers& Q)
{ return (std::to_string(Q.twoJ()) + " " + std::to_string(Q.P()) + " " + std::to_string(Q.C()) + " " + std::to_string(Q.twoI()) + " " + std::to_string(Q.G()) + " " + std::to_string(Q.Q())); }

/// Overload << operator
inline std::ostream& operator<< (std::ostream& os, const QuantumNumbers& Q)
{ os << "JP" << ((Q.C() == 0) ? "" : "C") << " = " << to_string(Q); return os; }

}

#endif
