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

#ifndef yap_DataSet_h
#define yap_DataSet_h

#include "DataPointBase.h"
#include "DataSetBase.h"
#include "Exceptions.h"

#include <vector>

namespace yap {

class Model;
class StatusManager;

template <typename T>
class DataSet : public DataSetBase
{
public:

    /// Constructor
    DataSet(const Model& m) : DataSetBase(m)
    {}

    /// add single empty data point
    void addEmptyPoint()
    {
        if (!model())
            throw exceptions::Exception("Model unset or deleted", "DataSet::add");
        DataPoints_.emplace_back(*this);

        auto& d = DataPoints_.back();

        if (!consistent(d)) {
            DataPoints_.pop_back();
            throw exceptions::Exception("produced inconsistent data point", "Model::addDataPoint");
        }
    }

    /// \return iterator to front of set
    const DataIterator& begin() const override
    { return const_cast<DataSet*>(this)->setBegin(const_cast<DataPointVector*>(&DataPoints_)->begin()); }

    /// \return iterator to end of set
    const DataIterator& end() const override
    { return const_cast<DataSet*>(this)->setEnd(const_cast<DataPointVector*>(&DataPoints_)->end()); }

protected:

};

/// swap
template <typename T>
inline void swap(DataSet<T>& A, DataSet<T>& B)
{ A.swap(B); }

}

#endif
