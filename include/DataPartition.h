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

#include "DataPoint.h"
#include "DataSet.h"
#include "logging.h"

#include <algorithm>
#include <iterator>
#include <vector>

namespace yap {

class DataPartitionBase;

/// \class DataIterator
/// \brief Class for iterating over a #DataPartitionBase
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class DataIterator : virtual public std::forward_iterator_tag
{
public:

    DataIterator& operator++();

    DataPoint& operator*()
    { return *Iterator_; }

    const DataPoint& operator*() const
    { return *Iterator_; }

    bool operator!=(const DataIterator& it) const
    { return Iterator_ != it.Iterator_; }

    bool ownedBy(DataPartitionBase* dpb) const
    { return Owner_ == dpb; }

protected:

    friend DataPartitionBase;

    // Protected constructor
    DataIterator(DataPartitionBase* p, std::vector<DataPoint>::iterator& it)
        : Owner_(p), Iterator_(it)
    {}

    DataPartitionBase* Owner_;
    std::vector<DataPoint>::iterator Iterator_;

};


/// \class DataPartitionBase
/// \brief Class defining a partition of the DataSet
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class DataPartitionBase
{
public:

    /// Constructor
    DataPartitionBase(std::vector<DataPoint>::iterator begin, std::vector<DataPoint>::iterator end);

    const DataIterator& begin() const
    { return Begin_; }

    const DataIterator& end() const
    { return End_; }

    /// \return DataPartition's index
    unsigned index() const
    { return DataPartitionIndex_; }

    /// clone the DataPartition object
    virtual DataPartitionBase* clone() = 0;

protected:

    /// make copy constructor protected. use clone() instead
    DataPartitionBase(const DataPartitionBase&) = default;

    /// \name DataPartitionBase friends
    /// @{

    friend DataIterator;
    friend class InitialStateParticle;

    /// @}

    /// increment and
    /// \return if still in range
    virtual void increment(DataIterator& it) = 0;

    std::vector<DataPoint>::iterator& rawIterator(DataIterator& it)
    { return it.Iterator_; }

    void setIndex(unsigned i)
    { DataPartitionIndex_ = i; }

    DataIterator Begin_;
    DataIterator End_;

    unsigned DataPartitionIndex_;

};

/// \class DataPartitionWeave
/// \brief A set of data spaced over the range [B,E) with spacing S = [B+0S, B+1S, B+2S, B+3S, ..., E)
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class DataPartitionWeave : public DataPartitionBase
{
public:

    /// Constructor
    DataPartitionWeave(std::vector<DataPoint>::iterator begin, std::vector<DataPoint>::iterator end,
                       unsigned spacing)
        : DataPartitionBase(begin, end),
          Spacing_(spacing)
    {}

    /// clone the DataPartition object
    virtual DataPartitionBase* clone() override
    { return new DataPartitionWeave(*this); }

protected:

    virtual void increment(DataIterator& it) override;

    unsigned Spacing_;

};

/// \class DataPartitionBlock
/// \brief A contiguous block of data
/// \author Johannes Rauch, Daniel Greenwald
/// \ingroup Data

class DataPartitionBlock : public DataPartitionWeave
{
public:

    /// Constructor
    DataPartitionBlock(std::vector<DataPoint>::iterator begin, std::vector<DataPoint>::iterator end)
        : DataPartitionWeave(begin, end, 1)
    {}

    /// clone the DataPartition object
    virtual DataPartitionBase* clone() override
    { return new DataPartitionBlock(*this); }

};

/// DataPartitionCreators

/// create DataPartitionWeave'es
/// \param dataSet The dataSet
/// \param nPartitions number of partitions to divide the dataSet into
std::vector<std::unique_ptr<DataPartitionBase> > createDataPartitionsWeave(DataSet& dataSet, unsigned nPartitions);

/// create DataPartitionBlock's
/// \param dataSet The dataSet
/// \param maxBlockSize maximum size of DataPoints a block can have
std::vector<std::unique_ptr<DataPartitionBase> > createDataPartitionsBlockBySize(DataSet& dataSet, unsigned maxBlockSize);

/// create DataPartitionBlock's
/// \param dataSet The dataSet
/// \param nPartitions number of blocks to divide the dataSet into
std::vector<std::unique_ptr<DataPartitionBase> > createDataPartitionsBlock(DataSet& dataSet, unsigned nPartitions);


}

#endif
