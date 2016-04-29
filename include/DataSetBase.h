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

#include "fwd/DataPointBase.h"

#include "DataPartition.h"
#include "FourVector.h"

namespace yap {

class Model;
class StatusManager;

class DataSetBase : public DataPartitionBlock
{
protected:

    /// Constructor
    DataSetBase(const Model& m);

public:

    /// Copy constructor
    DataSetBase(const DataSetBase& other);

    /// Move constructor
    DataSetBase(DataSetBase&& other);

    /// Copy assignment operator
    DataSetBase& operator=(const DataSetBase& other);

    /// Move assignment operator
    DataSetBase& operator=(DataSetBase&& other);

    /// Destructor
    virtual ~DataSetBase();

    /// Swap
    void swap(DataSet& other);

    /// \return raw pointer to associated model
    const Model* model() const
    { return Model_; }

    /// \return global status manager
    StatusManager& globalStatusManager()
    { return GlobalStatusManager_; }

    /// Check if data point is consisent with data set
    bool consistent(const DataPointBase* d) const;

    /// access by index
    DataPointBase* operator[](size_t i)
    { return DataPoints_[i]; }

    /// access by index (with check)
    DataPointBase* at(size_t i)
    { return DataPoints_.at(i); }

    /// access back
    DataPointBase* back()
    { return DataPoints_.back(); }

    size_t size() const
    { return DataPoints_.size(); }

    /// const access to DataPoints_
    const DataPointVector& points() const
    { return DataPoints_; }

    /// reserve storage space
    void reserve(size_t n)
    { DataPoints_.reserve(n); }

    /// call shrink to fit on DataPoints_
    void shrink_to_fit()
    { DataPoints_.shrink_to_fit(); }

    /// add single empty data point
    virtual void addEmptyPoint() = 0;

    /// add empty data points
    /// \param n number of points to add
    /// \tparam T type of storage inside data point
    void addEmptyPoints(size_t n)
    { for (size_t i = 0; i < n; ++i) addEmptyPoint(); }

    /// Add data point
    /// \param P vector of four momenta
    void add(const std::vector<FourVector<double> >& P);

    /// grant friend status to DataPartition to access non-const dataPoints()
    friend class DataPartition;

protected:

    /// non-const access to DataPoints_
    DataPointVector& dataPoints()
    { return DataPoints_; }

private:

    /// sets this as owner of all its data points
    void assertDataPointOwnership();

    /// Associated model
    const Model* Model_;

    /// Global StatusManager
    StatusManager GlobalStatusManager_;

    /// vector of data points contained in set
    DataPointVector DataPoints_;

};

#endif
