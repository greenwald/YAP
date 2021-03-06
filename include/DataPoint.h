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

#include <FourVector.h>

#include <complex>
#include <memory>
#include <set>
#include <vector>

namespace yap {

class CachedDataValue;
class DataAccessor;
class FourMomenta;
class HelicityAngles;
class MeasuredBreakupMomenta;
class ParticleCombination;

/// \class DataPoint
/// \brief Class for holding data and cached values per data point for fast calculation
/// \author Johannes Rauch, Daniel Greenwald
/// \defgroup Data Data-related classes

class DataPoint
{
public:

    /// \name Constructors
    /// @{

    /// 4-momenta constructor
    DataPoint(const std::vector<FourVector<double> >& P);

    // /// Invariant mass constructor
    //DataPoint(const std::map<std::shared_ptr<ParticleCombination>, double>& m2);

    /// @}

    const std::vector<FourVector<double> >& finalStateFourMomenta()
    { return FSPFourMomenta_; }

    bool setFinalStateFourMomenta(const std::vector<FourVector<double> >& fourMomenta);

    /// print information about the size of the DataPoint and its members
    void printDataSize();

    /// reserve space in Data_
    void allocateStorage(const FourMomenta& fourMom, const std::set<DataAccessor*> dataAccessors);

protected:

    /// \name DataPoint friends
    /// @{

    friend class CachedDataValue;
    friend class DataAccessor;
    friend class DataSet;
    friend class FourMomenta;

    /// @}

    /// vector of 4-momenta of final-state particles in event
    std::vector<FourVector<double> > FSPFourMomenta_;

    /// Vector of 4-momenta of non-final-state particles in event
    std::vector<FourVector<double> > FourMomenta_;

    /// Data storage for all DataAccessors
    /// first index is for the DataAccessor
    /// second index is for the symmeterization state (as known by the DataAccessor)
    /// third index is internal to the DataAccessor
    std::vector<std::vector<std::vector<double> > > Data_;

};

}

#endif
