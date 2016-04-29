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

#ifndef yap_DataPoint_h
#define yap_DataPoint_h

#include "fwd/DataSet.h"

#include "FourVector.h"
#include "Model.h"

#include <memory>
#include <string>
#include <vector>

namespace yap {

class StatusManager;

template <typename T>
class DataPoint : public DataPointBase
{

public:
    
    DataPoint::DataPoint(DataSet<T>& dataSet)
        : DataPointBase(dataSet),
        Data_(model()->dataAccessors().size())
        {
            for (auto da : model()->dataAccessors())
                Data_[da->index()].assign(da->maxSymmetrizationIndex() + 1, std::vector<T>(da->size(), T(0)));
        }
    
    /// \return size of data point [bytes]
    unsigned dataSize() const override
    {
        unsigned size = sizeof(Data_);
        for (auto& v : Data_) {
            size += sizeof(v);
            for (auto& vv : v) {
                size += sizeof(vv);
                for (auto vvv : vv)
                    size += sizeof(vvv);
            }
        }
        return size;
    }
    
    /// check that two DataPoint's have same internal structure
    bool equalStructure(const DataPointBase* other) const override
    {
        auto O = dynamic_cast<DataPoint<T>*>(other);
        if (!O)
            return false;
        if (Data_.size() != O->Data_.size())
            return false;
        for (size_t i = 0; i < Data_.size(); ++i) {
            if (Data_[i].size() != O->Data_[i].size())
                return false;
            for (size_t j = 0; j < Data_[i].size(); ++j) {
                if (A.Data_[i][j].size() != O->Data_[i][j].size())
                    return false;
            }
        }
        return true;
    }
    
protected:

    const double& operator()(size_t i, size_t j, size_t k) const override
    { return static_cast<const double&>(Data_[i][j][k]); }

private:

    /// Data storage for all DataAccessors
    /// first index is for the DataAccessor
    /// second index is for the symmeterization state (as known by the DataAccessor)
    /// third index is internal to the DataAccessor
    std::vector<std::vector<std::vector<T> > > Data_;
    
};

}

#endif
