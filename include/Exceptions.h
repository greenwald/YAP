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

#ifndef yap_Exceptions_h
#define yap_Exceptions_h

#include <stdexcept>
#include <string>

namespace yap {

namespace exceptions {

/// \class Exception
/// \brief Base class for handling YAP exceptions
/// \defgroup Exceptions YAP exceptions
///
/// Use this class for creating exceptions with clarifying text, which
/// will not be caught by specific type. Use the other exception classes
/// inheriting from this, when an exception must be caught by specific
/// type.
class Exception : public std::exception {
public:
    /// Constructor
    /// \param what_arg String descripting exception
    /// \param func_nam Name of function originating exception
    Exception(const std::string& what_arg, const std::string& func_name) : std::exception(), What_(what_arg)
    { addFunc(func_name); }

    /// add to func name
    void addFunc(const std::string& s) noexcept
    {
        if (s.empty()) return;
        if (!What_.empty()) What_ += " ";
        What_ += "from "  + s;
    }

    /// \return what string (data)
    const char* what() const noexcept override
    { return What_.data(); }

protected:
    Exception() : std::exception() {}
    std::string What_;
};

/// \class AngularMomentumNotConserved
/// \ingroup Exceptions
struct AngularMomentumNotConserved : public Exception {
AngularMomentumNotConserved(const std::string& func_name)
    : Exception("Angular momentum not conserved", func_name) {}
};

/// \class InconsistentSpinProjection
/// \ingroup Exceptions
struct InconsistentSpinProjection : public Exception {
InconsistentSpinProjection(const std::string& what_arg, const std::string& func_name)
    : Exception(what_arg, func_name) {}
};

/// \class OutsidePhaseSpace
/// \ingroup Exceptions
struct OutsidePhaseSpace : public Exception {
OutsidePhaseSpace(const std::string& what_arg = "", const std::string& func_name = "FourMomenta::setSquaredMasses")
    : Exception(what_arg, func_name) {}
};

/// \class NotTwoBodyParticleCombination
/// \ingroup Exceptions
struct NotTwoBodyParticleCombination : public Exception {
NotTwoBodyParticleCombination(const std::string& what_arg = "", const std::string& func_name = "")
    : Exception(what_arg, func_name) {}
};

/// \class ParameterIsFixed
/// \ingroup Exceptions
struct ParameterIsFixed : public Exception {
ParameterIsFixed(const std::string& what_arg = "", const std::string& func_name = "")
    : Exception(what_arg, func_name) {}
};

/// \class NonfiniteResult
/// \ingroup Exceptions
struct NonfiniteResult : public Exception {};

/// \class InconsistentDataPoint
/// \ingroup Exceptions
struct InconsistentDataPoint : public Exception {
InconsistentDataPoint(const std::string& func_name = "")
    : Exception("Inconsistent DataPoint", func_name) {}
};

/// \class EmptyFourMomentaVector
/// \ingroup Exceptions
struct EmptyFourMomentaVector : public Exception {
EmptyFourMomentaVector(const std::string& func_name = "")
    : Exception("Empty FourMomenta vector", func_name) {}
};

}

}

#endif
