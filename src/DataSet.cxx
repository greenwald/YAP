#include "DataSet.h"

#include "DataPoint.h"
#include "FourMomenta.h"
#include "logging.h"

namespace yap {

//-------------------------
bool DataSet::consistent(const DataPoint& d) const
{
    bool result = true;

    if (!empty() && d.FourMomenta_.size() != (*this)[0].FourMomenta_.size()) {
        LOG(ERROR) << "DataSet::consistent(DataPoint) - DataPoint has wrong number of four momenta (" << d.FourMomenta_.size() << " != " << (*this)[0].FourMomenta_.size() << ")";
        result = false;
    }

    return result;
}

}
