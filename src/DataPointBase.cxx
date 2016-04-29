#include "DataPointBase.h"

#include "DataSetBase.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "StaticDataAccessor.h"

namespace yap {

//-------------------------
DataPointBase::DataPointBase(DataSetBase& dataSetBase) :
    DataSetBase_(&dataSetBase)
{
    if (!model())
        throw exceptions::Exception("Model unset", "DataPointBase::DataPointBase");
}

//-------------------------
const Model* DataPointBase::model() const
{
    return DataSetBase_->model();
}

//-------------------------
void DataPointBase::setFinalStateMomenta(const std::vector<FourVector<double> >& P, StatusManager& sm)
{
    if (!model())
        throw exceptions::Exception("Model unset", "DataPoint::setFinalStateMomenta");

    model()->fourMomenta()->setFinalStateMomenta(*this, P, sm);

    // call calculate on all static data accessors in model
    for (auto& sda : model()->staticDataAccessors())
        sda->calculate(*this, sm);

    FDEBUG("set final state momenta and calculated StaticDataAccessors")
}

//-------------------------
void DataPointBase::setFinalStateMomenta(const std::vector<FourVector<double> >& P)
{
    setFinalStateMomenta(P, *DataSet_);
}

}
