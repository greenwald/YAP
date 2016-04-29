#include "DataSetBase.h"

#include "DataPointBase.h"
#include "Model.h"

namespace yap {

//-------------------------
DataSetBase::DataSetBase(const Model& m) :
    DataPartitionBlock(m.dataAccessors()),
    Model_(&m),
    GlobalStatusManager_(m.dataAccessors())
{
}

//-------------------------
DataSetBase::DataSetBase(const DataSetBase& other) :
    DataPartitionBlock(other),
    DataPoints_(other.DataPoints_),
    Model_(other.Model_),
    GlobalStatusManager_(other.GlobalStatusManager_)
{
    assertDataPointOwnership();
}

//-------------------------
DataSetBase::DataSetBase(DataSetBase&& other) :
    DataPartitionBlock(std::move(other)),
    DataPoints_(std::move(other.DataPoints_)),
    Model_(std::move(other.Model_)),
    GlobalStatusManager_(std::move(other.GlobalStatusManager_))
{
    assertDataPointOwnership();
}

//-------------------------
DataSetBase& DataSetBase::operator=(const DataSetBase& other)
{
    DataPartitionBlock::operator=(other);
    Model_ = other.Model_;
    DataPoints_ = other.DataPoints_;
    GlobalStatusManager_ = other.GlobalStatusManager_;
    assertDataPointOwnership();
    return *this;
}

//-------------------------
DataSetBase& DataSetBase::operator=(DataSetBase&& other)
{
    DataPartitionBlock::operator=(std::move(other));
    Model_ = std::move(other.Model_);
    DataPoints_ = std::move(other.DataPoints_);
    GlobalStatusManager_ = std::move(other.GlobalStatusManager_);
    assertDataPointOwnership();
    return *this;
}

//-------------------------
DataSetBase::~DataSetBase()
{
    // delete data points
    for (auto& d : DataPoints_)
        delete d;
}

//-------------------------
void DataSetBase::swap(DataSetBase& other)
{
    using std::swap;
    swap(static_cast<DataPartitionBlock&>(*this), static_cast<DataPartitionBlock&>(other));
    swap(Model_, other.Model_);
    swap(DataPoints_, other.DataPoints_);
    assertDataPointOwnership();
    other.assertDataPointOwnership();
}

//-------------------------
bool DataSetBase::consistent(const DataPointBase* d) const
{ return DataPoints_.empty() or DataPoints_.front()->equalStructure(d); }

//-------------------------
void DataSetBase::add(const std::vector<FourVector<double> >& P)
{
    addEmptyPoint();
    DataPoints_.back()->setFinalStateMomenta(P);
}

//-------------------------
void DataSetBase::assertDataPointOwnership()
{
    for (auto& d : DataPoints_)
        d->DataSet_ = this;
}

}
