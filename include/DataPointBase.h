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

#ifndef yap_DataSetBase_h
#define yap_DataSetBase_h

#include "fwd/DataSetBase.h"

#include "FourVector.h"

namespace yap {

class Model;
class StatusManager;

class DataPointBase
{
protected:

    virtual const double& operator()(size_t i, size_t j, size_t k) const = 0;

    /// non-const access
    double& operator()(size_t i, size_t j, size_t k)
    { return const_cast<double&>(static_cast<const DataPointBase*>(this)->operator()(i, j, k)); }
    
    /// protected constructor
    DataPointBase() = default;
    
    /// \param dataSet DataSet this DataPoint belongs to
    DataPointBase(DataSetBase& dataSetBase);
    
public:
    
    /// set four momenta of data point
    /// \param P vector of FourVectors of final-state momenta
    /// \param sm StatusManager to update
    void setFinalStateMomenta(const std::vector<FourVector<double> >& P, StatusManager& sm);
    
    /// set four momenta of data point
    /// \param P vector of FourVectors of final-state momenta
    void setFinalStateMomenta(const std::vector<FourVector<double> >& P);
    
    /// \return size of data point [bytes]
    virtual unsigned dataSize() const = 0;

    /// check that two DataPoint's have same internal structure
    virtual bool equalStructure(const DataPointBase* B) const = 0;

    /// \return model to which DataPoint belongs
    const Model* model() const;

    /// grant friend status to CachedDataValue to access operator()
    friend class CachedDataValue;

    /// grant friend status to DataSetBase to set itself owner
    friend class DataSetBase;

private:

    /// raw pointer to owning DataSetBase
    DataSetBase* DataSetBase_;

};

/// \typedef DataPointVector
/// \brief stl vector of DataPoint's
using DataPointVector = std::vector<DataPointBase*>;

}

#endif
