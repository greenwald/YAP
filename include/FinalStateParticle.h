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

#ifndef yap_FinalStateParticle_h
#define yap_FinalStateParticle_h

#include "Particle.h"

namespace yap {

/// \ingroup Particle

class FinalStateParticle : public Particle
{
public:
    FinalStateParticle();
    ~FinalStateParticle();

    virtual Amp amplitude(DataPoint& d) {return Amp(1.);}
    //virtual bool consistent() const {return true;}

private:
    int PDGNumber_; /// PDG code of the particle
};

}

#endif