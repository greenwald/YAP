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

#ifndef yap_PHSP_h
#define yap_PHSP_h

#include "fwd/FourVector.h"

#include "FourMomenta.h"
#include "MassRange.h"
#include "Model.h"

#include <algorithm>
#include <random>
#include <vector>

namespace yap {

/// \return vector of four momenta for daughters uniformly randomly generated in phase space of model
/// \param M model to calculate with
/// \param initial_mass Initial mass to decay from
/// \param A mass axes
/// \param R2 vector of squared mass ranges of axes
/// \param g random generator to pass to uniform_real_distribution
/// \param max_attempts maximum number of attempts to make to find a point in phase space
template <class Generator>
const std::vector<FourVector<double> > phsp(const Model& M, double initial_mass, const MassAxes& A, const std::vector<MassRange>& R2, Generator& g, unsigned max_attempts)
{
    static std::uniform_real_distribution<double> uniform(0, std::nextafter(1., 2.));

    // create vector to store invariant masses in
    std::vector<double> m2(R2.size(), 0);
    std::vector<FourVector<double> > P;

    for (unsigned n = 0; n < max_attempts && P.empty(); ++n) {
        // generate random point in hypercube of mass ranges
        std::transform(R2.begin(), R2.end(), m2.begin(), [&](const MassRange & r2) {return r2[0] + (r2[1] - r2[0]) * uniform(g);});
        P = calculate_four_momenta(initial_mass, M, A, m2);
    }
    return P;
}

}

#endif
