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

#ifndef yap_CachedDataValue_h
#define yap_CachedDataValue_h

#include "fwd/DataPointBase.h"

#include "FourVector.h"
#include "Parameter.h"

#include <memory>
#include <set>
#include <utility>
#include <vector>

namespace yap {

class CachedDataValue;
enum class CalculationStatus : bool;
class DataAccessor;
class ParticleCombination;
class StatusManager;
enum class VariableStatus;

/// \typedef CachedDataValueSet
/// \ingroup Data
/// \ingroup Cache
using CachedDataValueSet = std::set<std::shared_ptr<CachedDataValue> >;

/// \typedef CachedDataValueVector
/// \ingroup Data
/// \ingroup Cache
using CachedDataValueVector = std::vector<std::shared_ptr<CachedDataValue> >;

/// \struct DaughterCachedDataValue
/// \brief Stores a shared_ptr to a CachedDataValue and the index of
/// the daughter within a Particle Combination to pass to it
/// \ingroup Data
/// \ingroup Cache
struct DaughterCachedDataValue {
    /// constructor
    DaughterCachedDataValue(std::shared_ptr<CachedDataValue> cdv, size_t index)
        : CDV(cdv), Daughter(index) {}

    /// CachedDataValue
    std::shared_ptr<CachedDataValue> CDV;

    /// index of daughter to pass to it
    size_t Daughter;

    friend bool operator<(const DaughterCachedDataValue& A, const DaughterCachedDataValue& B)
    { return (A.CDV < B.CDV) or (A.Daughter < B.Daughter); }
};

/// \typedef DaughterCachedDataValueSet
/// \ingroup Data
/// \ingroup Cache
using DaughterCachedDataValueSet = std::set<DaughterCachedDataValue>;

/// \class CachedDataValueBase
/// \brief Class for managing cached values inside a #DataPoint
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class CachedDataValue : public std::enable_shared_from_this<CachedDataValue>
{
public:

    /// \struct stores calculation and variable statuses for a CachedDataValue
    struct Status {
        /// constructor
        Status();

        /// Calculation status
        CalculationStatus Calculation;

        /// Variable status
        VariableStatus Variable;

        /// assignment of Calculation
        void operator=(const CalculationStatus& s)
        { Calculation = s; }

        /// assignment of Variable
        /// Does not change Variable if it is fixed!
        void operator=(const VariableStatus& s);

        /// equality operator for checking the CalculationStatus
        friend bool operator==(const Status& S, const CalculationStatus& s)
        { return S.Calculation == s; }

        /// inequality operator for checking the CalculationStatus
        friend bool operator!=(const Status& S, const CalculationStatus& s)
        { return S.Calculation != s; }

        /// equality operator for checking the VariableStatus
        friend bool operator==(const Status& S, const VariableStatus& s)
        { return S.Variable == s; }

        /// inequality operator for checking the VariableStatus
        friend bool operator!=(const Status& S, const VariableStatus& s)
        { return S.Variable != s; }

        friend std::ostream& operator<<(std::ostream& str, const CachedDataValue::Status& S);
    };

    /// \name Managing dependencies
    /// @{

    /// add Parameters this CachedDataValue depends on
    void addDependencies(ParameterSet deps)
    { for (auto& dep : deps) addDependency(dep); }

    /// add Parameter this CachedDataValue depends on
    void addDependency(std::shared_ptr<ParameterBase> dep)
    { ParametersItDependsOn_.insert(dep); }

    /// add CachedDataValue's this CachedDataValue depends on
    void addDependencies(CachedDataValueSet deps)
    { for (auto& dep : deps) addDependency(dep); }

    /// add CachedDataValue's of a daughter this CachedDataValue depends on
    void addDependencies(DaughterCachedDataValueSet deps)
    { for (auto& dep : deps) addDependency(dep); }

    /// add CachedDataValue this CachedDataValue depends on
    void addDependency(std::shared_ptr<CachedDataValue> dep)
    { CachedDataValuesItDependsOn_.insert(dep); }

    /// add CachedDataValue of a daughter this CachedDataValue depends on
    void addDependency(DaughterCachedDataValue dep)
    { DaughterCachedDataValuesItDependsOn_.insert(dep); }

    /// remove dependency
    void removeDependency(std::shared_ptr<ParameterBase> dep);

    /// remove dependencies
    void removeDependencies(ParameterSet deps)
    { for (auto& dep : deps) removeDependency(dep); }

    /// check for dependency
    bool dependsOn(std::shared_ptr<ParameterBase> dep) const
    { return ParametersItDependsOn_.find(dep) != ParametersItDependsOn_.end(); }

    /// check for dependency
    bool dependsOn(std::shared_ptr<CachedDataValue> dep) const
    { return CachedDataValuesItDependsOn_.find(dep) != CachedDataValuesItDependsOn_.end(); }

    /// @}

    /// \name Getters
    /// @{

    /// \return raw pointer owning DataAccessor
    DataAccessor* owner() const
    { return Owner_; }

    /// \return index within owner
    int index() const
    { return Index_; }

    /// Get value from #DataPointBase for particular symmetrization
    /// \param index index of value to get from within cached value (must be less than #Size_)
    /// \param d #DataPointBase to get value from
    /// \param sym_index index of symmetrization to grab from
    /// \return Value of CachedDataValue inside the data point
    double value(unsigned index, const DataPointBase& d, unsigned sym_index) const;

    /// \return Size of cached value (number of real elements)
    virtual unsigned size() const
    { return Size_; }

    /// \return set of Parameters on which this object depends
    const ParameterSet& parameterDependencies() const
    { return ParametersItDependsOn_; }

    /// \return set of CachedDataValues on which this object depends
    const CachedDataValueSet& cachedDataValueDependencies() const
    { return CachedDataValuesItDependsOn_; }

    /// \return set of DaughterCachedDataValues on which this object depends
    const DaughterCachedDataValueSet& daughterCachedDataValueDependencies() const
    { return DaughterCachedDataValuesItDependsOn_; }

    /// @}

    /// \name Setters
    /// @{

    /// Set value into #DataPointBase for particular symmetrization
    /// (No update to VariableStatus or CalculationStatus is made!)
    /// \param index index of value to get from within cached value (must be less than #Size_)
    /// \param val Value to set to
    /// \param d #DataPointBase to update
    /// \param sym_index index of symmetrization to apply to
    void setValue(unsigned index, double val, DataPointBase& d, unsigned sym_index) const;

    /// @}

    /// grant friend status to DataAccessor to set itself owner
    friend class DataAccessor;

protected:

    /// Constructor (protected)
    CachedDataValue(unsigned size, ParameterSet pars = {}, CachedDataValueSet vals = {});

    /// set the owning DataAccessor
    void setDataAccessor(DataAccessor* owner);

    /// set index
    void setIndex(int i)
    { Index_ = i; }

    /// set position
    void setPosition(int p)
    { Position_ = p; }

private:

    /// Owning DataAccessor
    DataAccessor* Owner_;

    /// index within owner
    int Index_;

    /// Position of first element of cached value within data vector
    int Position_;

    /// Size of cached value (number of real elements)
    unsigned Size_;

    ParameterSet ParametersItDependsOn_;
    CachedDataValueSet CachedDataValuesItDependsOn_;
    DaughterCachedDataValueSet DaughterCachedDataValuesItDependsOn_;

};

/// \class RealCachedDataValue
/// \brief Class for managing a single real cached value inside a #DataPointBase
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class RealCachedDataValue : public CachedDataValue
{
public:

    /// create shared_ptr to RealCachedDataValue
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    static std::shared_ptr<RealCachedDataValue> create(DataAccessor* da, ParameterSet pars = {}, CachedDataValueSet vals = {});

    /// Set value into #DataPointBase for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val Value to set to
    /// \param d #DataPointBase to update
    /// \param sym_index index of symmetrization to apply to
    /// \param sm StatusManager
    void setValue(double val, DataPointBase& d, unsigned sym_index, StatusManager& sm) const;

    /// Get value from #DataPointBase for particular symmetrization
    /// \param d #DataPointBase to get value from
    /// \param sym_index index of symmetrization to grab from
    /// \return Value of CachedDataValue inside the data point
    double value(const DataPointBase& d, unsigned sym_index) const
    { return CachedDataValue::value(0, d, sym_index); }

private:

    /// Constructor
    RealCachedDataValue(ParameterSet pars = {}, CachedDataValueSet vals = {})
        : CachedDataValue(1, pars, vals) {}

};

/// \class ComplexCachedDataValue
/// \brief Class for managing a complex cached value inside a #DataPointBase
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class ComplexCachedDataValue : public CachedDataValue
{
public:

    /// create shared pointer to ComplexCachedDataValue
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    static std::shared_ptr<ComplexCachedDataValue> create(DataAccessor* da, ParameterSet pars = {}, CachedDataValueSet vals = {});

    /// Set value into #DataPointBase for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val Value to set to
    /// \param d #DataPointBase to update
    /// \param sym_index index of symmetrization to apply to
    /// \param sm StatusManager
    void setValue(std::complex<double> val, DataPointBase& d, unsigned sym_index, StatusManager& sm) const
    { setValue(real(val), imag(val), d, sym_index, sm); }

    /// Set value into #DataPointBase for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val_re real part of value to set to
    /// \param val_im imaginary part of value to set to
    /// \param d #DataPointBase to update
    /// \param sym_index index of symmetrization to apply to
    /// \param sm StatusManager
    void setValue(double val_re, double val_im, DataPointBase& d, unsigned sym_index, StatusManager& sm) const;

    /// Get value from #DataPointBase for particular symmetrization
    /// \param d #DataPointBase to get value from
    /// \param sym_index index of symmetrization to grab from
    /// \return Value of CachedDataValue inside the data point
    std::complex<double> value(const DataPointBase& d, unsigned  sym_index) const
    { return std::complex<double>(CachedDataValue::value(0, d, sym_index), CachedDataValue::value(1, d, sym_index)); }

private:

    /// Constructor (protected)
    /// see #create for details
    ComplexCachedDataValue(ParameterSet pars = {}, CachedDataValueSet vals = {})
        : CachedDataValue(2, pars, vals) {}


};

/// \class FourVectorCachedDataValue
/// \brief Class for managing a four-vector cached value inside a #DataPointBase
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
/// \ingroup Cache
class FourVectorCachedDataValue : public CachedDataValue
{
public:

    /// create shared pointer to ComplexCachedDataValue
    /// \param owner #DataAccessor to which this cached value belongs
    /// \param pars set of shared pointers to Parameters cached value depends on
    /// \param vals set of shared pointers to CachedValues cached value depends on
    static std::shared_ptr<FourVectorCachedDataValue> create(DataAccessor* da, ParameterSet pars = {}, CachedDataValueSet vals = {});

    /// Set value into #DataPointBase for particular symmetrization, and
    /// update VariableStatus for symm. and partition index
    /// \param val Value to set to
    /// \param d #DataPointBase to update
    /// \param sym_index index of symmetrization to apply to
    /// \param sm StatusManager
    void setValue(FourVector<double> val, DataPointBase& d, unsigned sym_index, StatusManager& sm) const;

    /// Get value from #DataPointBase for particular symmetrization
    /// \param d #DataPointBase to get value from
    /// \param sym_index index of symmetrization to grab from
    /// \return Value of CachedDataValue inside the data point
    FourVector<double> value(const DataPointBase& d, unsigned  sym_index) const
    {
        return FourVector<double>( { CachedDataValue::value(0, d, sym_index),
                                     CachedDataValue::value(1, d, sym_index),
                                     CachedDataValue::value(2, d, sym_index),
                                     CachedDataValue::value(3, d, sym_index)
                                   });
    }

private:

    /// Constructor (protected)
    /// see #create for details
    FourVectorCachedDataValue(ParameterSet pars = {}, CachedDataValueSet vals = {})
        : CachedDataValue(4, pars, vals) {}


};

}

#endif
