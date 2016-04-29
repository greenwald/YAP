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

#ifndef yap_DataPartition_h
#define yap_DataPartition_h

#include "fwd/DataPointBase.h"
#include "fwd/DataSetBase.h"

#include "StatusManager.h"

#include <algorithm>
#include <iterator>
#include <memory>
#include <vector>

namespace yap {

class DataPartition;

/// \class DataIterator
/// \brief Class for iterating over a #DataPartition
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class DataIterator : virtual public std::forward_iterator_tag
{
public:

    /// increment operator
    DataIterator& operator++();

    /// dereference operator
    DataPointBase* operator*()
    { return *Iterator_; }

    /// dereference operator (const)
    const DataPointBase* operator*() const
    { return *Iterator_; }

    /// inequality operator
    bool operator!=(const DataIterator& it) const
    { return Iterator_ != it.Iterator_; }

    /// check ownership
    bool ownedBy(const DataPartition& dp) const
    { return Partition_ == &dp; }

    /// grant friend status to DataPartition to access Iterator_
    friend DataPartition;

protected:

    /// constructor with defaulted iterator
    /// \param p owning DataPartition
    DataIterator(const DataPartition& p)
        : Partition_(&p) {}

    /// constructor
    /// \param p owning DataPartition
    /// \param it vector<DataPoint> iterator to contain
    DataIterator(const DataPartition& p, DataPointVector::iterator it)
        : Partition_(&p), Iterator_(it) {}

    /// owning DataPartition
    const DataPartition* Partition_;

    /// iterator within vector<DataPoint>
    DataPointVector::iterator Iterator_;

};


/// \class DataPartition
/// \brief Class defining a partition of the DataSetBase
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
class DataPartition : public StatusManager
{
protected:

    /// Constructor
    /// \param sm StatusManager to copy StatusManager structure from
    /// \param begin vector<DataPoint>::iterator of start
    /// \param end vector<DataPoint>::iterator of end
    DataPartition(const StatusManager& sm, DataPointVector::iterator begin, DataPointVector::iterator end)
        : StatusManager(sm), Begin_(*this, begin), End_(*this, end) {}

    /// constructor taking a DataAccessorSet
    /// \param sDA DataAccessorSet to initialize StatusManager from
    DataPartition(const DataAccessorSet& sDA)
        : StatusManager(sDA), Begin_(*this), End_(*this) {}

public:

    /// virtual destructor (defaulted)
    virtual ~DataPartition() = default;

    /// copy constructor (defaulted)
    DataPartition(const DataPartition&) = default;

    /// move constructor (defaulted)
    DataPartition(DataPartition&&) = default;

    /// copy assignment operator (defaulted)
    DataPartition& operator=(const DataPartition&) = default;

    /// move assignment operator (defaulted)
    DataPartition& operator=(DataPartition&&) = default;

    /// \return begin iterator
    virtual const DataIterator& begin() const
    { return Begin_; }

    /// \return end iterator
    virtual const DataIterator& end() const
    { return End_; }

    /// grant friend status to DataIterator to call increment
    friend DataIterator;

    /// grant friend status to Model to call ...
    friend class Model;

protected:

    /// increment iterator;
    /// Must be overloaded in derived classes
    virtual void increment(DataIterator& it) const
    { it = End_; }

    /// \return vector<DataPoint> iterator inside DataIterator
    DataPointVector::iterator& rawIterator(DataIterator& it) const
    { return it.Iterator_; }

    /// set begin
    const DataIterator& setBegin(DataPointVector::iterator it)
    { Begin_ = DataIterator(*this, it); return Begin_; }

    /// set end
    const DataIterator& setEnd(DataPointVector::iterator it)
    { End_ = DataIterator(*this, it); return End_; }

    /// get non-const begin from DataSetBase
    static DataPointVector::iterator begin(DataSetBase& ds);

    /// get non-const end from DataSetBase
    static DataPointVector::iterator end(DataSetBase& ds);

private:

    /// begin DataIterator
    DataIterator Begin_;

    /// end DataIterator
    DataIterator End_;
    
};

/// \typedef DataPartitionVector
/// \brief Vector of unique_ptr's to DataPartition
using DataPartitionVector = std::vector<std::unique_ptr<DataPartition> >;

/// \class DataPartitionBlock
/// \brief A contiguous block of data
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
class DataPartitionBlock : public DataPartition
{
public:

    /// Constructor
    /// \param sm StatusManager to copy StatusManager structure from
    /// \param begin vector<DataPoint>::iterator of start
    /// \param end vector<DataPoint>::iterator of end
    DataPartitionBlock(const StatusManager& sm, DataPointVector::iterator begin, DataPointVector::iterator end)
        : DataPartition(sm, begin, end) {}

    /// \return DataParitionVector covering DataSetBase as contiguous blocks
    /// \param dataSetBase the DataSetBase
    /// \param n number of partitions to divide dataSetBase into
    static DataPartitionVector create(DataSetBase& dataSetBase, unsigned n);

    /// \return DataParitionVector covering DataSetBase as contiguous blocks of specified size
    /// \param dataSetBase The dataSetBase
    /// \param s maximum size of partitions to divide the dataSetBase into
    static DataPartitionVector createBySize(DataSetBase& dataSetBase, size_t s);

protected:

    /// constructor taking a DataAccessorSet
    /// \param sDA DataAccessorSet to initialize StatusManager from
    DataPartitionBlock(const DataAccessorSet& sDA)
        : DataPartition(sDA) {}

    /// increment DataIterator
    /// \param it DataIterator to iterate
    virtual void increment(DataIterator& it) const override
    { ++rawIterator(it); }
};

/// \class DataPartitionWeave
/// \brief A set of data spaced over the range [B,E) with spacing S = [B+0S, B+1S, B+2S, B+3S, ..., E)
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data
class DataPartitionWeave : public DataPartition
{
public:

    /// Constructor
    /// \param sm StatusManager to copy StatusManager structure from
    /// \param begin vector<DataPoint>::iterator of start
    /// \param end vector<DataPoint>::iterator of end
    /// \param spacing Spacing between consecutively evaluated points
    DataPartitionWeave(const StatusManager& sm, DataPointVector::iterator begin, DataPointVector::iterator end, unsigned spacing)
        : DataPartition(sm, begin, end), Spacing_(spacing) {}

    /// \return DataParitionVector covering DataSetBase as a weave
    /// \param dataSetBase The dataSetBase
    /// \param n number of partitions to divide dataSetBase into
    static DataPartitionVector create(DataSetBase& dataSetBase, unsigned n);

protected:

    /// constructor taking a DataAccessorSet
    /// \param sDA DataAccessorSet to initialize StatusManager from
    /// \param spacing for weave
    DataPartitionWeave(const DataAccessorSet& sDA, unsigned spacing)
        : DataPartition(sDA), Spacing_(spacing) {}

    /// increment DataIterator
    /// \param it DataIterator to iterate
    virtual void increment(DataIterator& it) const override;

    /// spacing between data points for the weaving
    unsigned Spacing_;
};

}

#endif
