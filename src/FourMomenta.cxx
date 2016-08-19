#include "FourMomenta.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "container_utils.h"
#include "DecayingParticle.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "FourVector.h"
#include "logging.h"
#include "MassAxes.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "StatusManager.h"

#include <iomanip>

namespace yap {

//-------------------------
FourMomenta::FourMomenta(Model& m) :
    StaticDataAccessor(m, equal_by_orderless_content),
    TotalIndex_(-1),
    P_(FourVectorCachedValue::create(*this)),
    M_(RealCachedValue::create(*this))
{
}

//-------------------------
void FourMomenta::setFinalStateMomenta(DataPoint& d, const std::vector<FourVector<double> >& P, StatusManager& sm) const
{
    if (P.empty())
        throw exceptions::EmptyFourMomentaVector("FourMomenta::setFinalStateMomenta");

    if (P.size() != FSPIndices_.size())
        throw exceptions::Exception("Wrong number of momenta provided (" + std::to_string(P.size()) + " != " + std::to_string(FSPIndices_.size()) + ")",
                                    "FourMomenta::setFinalStateMomenta");

    for (size_t i = 0; i < P.size(); ++i)
        P_->setValue(P[i], d, FSPIndices_[i], sm);
}

//-------------------------
void FourMomenta::addParticleCombination(std::shared_ptr<ParticleCombination> pc)
{
    StaticDataAccessor::addParticleCombination(pc);
    auto index = symmetrizationIndex(pc);

    // check for ISP
    if (TotalIndex_ < 0 and pc->indices().size() == model()->finalStateParticles().size())
        TotalIndex_ = index;

    /// check for FSP
    if (is_final_state_particle_combination(*pc)) {
        if (pc->indices()[0] + 1 > FSPIndices_.size())
            FSPIndices_.resize(pc->indices()[0] + 1, -1);
        if (FSPIndices_[pc->indices()[0]] < 0)
            FSPIndices_[pc->indices()[0]] = index;
    }
}

//-------------------------
bool FourMomenta::consistent() const
{
    bool C = StaticDataAccessor::consistent();

    if (TotalIndex_ < 0) {
        FLOG(ERROR) << "ISP symmetrization index has not been recorded.";
        C &= false;
    }

    if (FSPIndices_.size() != model()->finalStateParticles().size() or std::any_of(FSPIndices_.begin(), FSPIndices_.end(), [](int i) {return i < 0;})) {
        FLOG(ERROR) << "FSP symmetrization indices not all recorded";
        C &= false;
    }

    return C;
}

//-------------------------
const FourVector<double> FourMomenta::totalMomentum(const DataPoint& d) const
{
    if (TotalIndex_ < 0)
        throw exceptions::Exception("Initial-state particle unknown", "FourMomenta::totalMomentum");
    return P_->value(d, TotalIndex_);
}

//-------------------------
const std::vector<FourVector<double> > FourMomenta::finalStateMomenta(const DataPoint& d) const
{
    std::vector<FourVector<double> > P;
    P.reserve(FSPIndices_.size());

    for (size_t i = 0; i  < FSPIndices_.size(); ++i) {
        if (FSPIndices_[i] < 0)
            throw exceptions::Exception("Final-state particle " + std::to_string(i) + "unknown", "FourMomenta::finalStateFourMomenta");
        P.push_back(P_->value(d, FSPIndices_[i]));
    }
    return P;
}

//-------------------------
FourVector<double> FourMomenta::p(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    return P_->value(d, symmetrizationIndex(pc));
}

//-------------------------
double FourMomenta::m(const DataPoint& d, const std::shared_ptr<ParticleCombination>& pc) const
{
    return M_->value(d, symmetrizationIndex(pc));
}

//-------------------------
void FourMomenta::calculate(DataPoint& d, StatusManager& sm) const
{
    // set all masses as uncalculated
    sm.set(*M_, CalculationStatus::uncalculated);

    // get fsp four momenta:
    auto fsp = finalStateMomenta(d);

    for (auto& kv : symmetrizationIndices()) {

        // check if calculation unnecessary
        if (sm.status(*M_, kv.second) == CalculationStatus::calculated)
            continue;

        const auto P = std::accumulate(kv.first->indices().begin(), kv.first->indices().end(), FourVector_0,
        [&](const FourVector<double>& V, unsigned i) {return V + fsp[i];});

        P_->setValue(P, d, kv.second, sm);
        M_->setValue(abs(P), d, kv.second, sm);
    }
}

//-------------------------
std::ostream& print_mp_string(std::ostream& os, unsigned n, unsigned m_p, std::shared_ptr<ParticleCombination> pc, double m, FourVector<double> p, double M = -1)
{
    os << std::setw(n) << indices_string(*pc) << " : "
       << "m = " << std::setprecision(m_p) << m << " GeV/c^2";
    if (M >= 0)
        os << " (nominally " << std::setprecision(m_p) << M << " GeV/c^2)";
    else
        os << "            " << std::setw(m_p) << " "                << "         ";
    os << "\tp = " << p << " GeV";
    return os;
}

//-------------------------
std::ostream& FourMomenta::printMasses(const DataPoint& d, std::ostream& os) const
{
    os << "Invariant masses:" << std::endl;
    unsigned n_fsp = model()->finalStateParticles().size();
    unsigned n = n_fsp + 2;
    unsigned m_p = 6;

    std::set<unsigned> used;

    // print the ISP
    for (auto& pc_i : symmetrizationIndices())
        // if ISP and not yet printed (should only print once)
        if (pc_i.first->indices().size() == n_fsp and used.find(pc_i.second) == used.end()) {
            os << "    ISP : ";
            print_mp_string(os, n, m_p, pc_i.first, m(d, pc_i.first), p(d, pc_i.first));
            os << std::endl;
            used.insert(pc_i.second);
        }

    // print the FSP's
    for (size_t i = 0; i < FSPIndices_.size(); ++i)
        for (auto& pc_i : symmetrizationIndices())
            if (is_final_state_particle_combination(*pc_i.first) and pc_i.first->indices()[0] == i and used.find(pc_i.second) == used.end()) {
                os << "    FSP : ";
                print_mp_string(os, n, m_p, pc_i.first, m(d, pc_i.first), p(d, pc_i.first), model()->finalStateParticles()[i]->mass());
                os << std::endl;
                used.insert(pc_i.second);
            }


    // print the rest in increasing number of particle content
    for (unsigned i = 2; i < n_fsp; ++i)
        for (auto& pc_i : symmetrizationIndices())
            // if i-particle mass and not yet printed
            if (pc_i.first->indices().size() == i and used.find(pc_i.second) == used.end()) {
                os << "    " << i << "-p : ";
                print_mp_string(os, n, m_p, pc_i.first, m(d, pc_i.first), p(d, pc_i.first));
                os << std::endl;
                used.insert(pc_i.second);
            }

    // print unused
    for (auto& pc_i : symmetrizationIndices())
        if (used.find(pc_i.second) == used.end()) {
            os << "        : ";
            print_mp_string(os, n, m_p, pc_i.first, m(d, pc_i.first), p(d, pc_i.first));
            os << std::endl;
            used.insert(pc_i.second);
        }

    return os;
}

//-------------------------
bool check_invariant_masses(const MassAxes& axes, const std::vector<double>& squared_masses, const std::vector<FourVector<double> >& fourMomenta)
{
    for (size_t i = 0; i < axes.size(); ++i) {
        auto p = std::accumulate(axes[i]->indices().begin(), axes[i]->indices().end(), FourVector_0, [&](const FourVector<double>& p, unsigned j) {return p + fourMomenta[j];});
        if (fabs(norm(p) - squared_masses[i]) > 5. * std::numeric_limits<double>::epsilon() )
            return false;
    }
    return true;
}

//-------------------------
std::vector<FourVector<double> > calculate_four_momenta(double initial_mass, const FinalStateParticleVector& FSPs,
                                                        const MassAxes& axes, const std::vector<double>& squared_masses)
{
    if (FSPs.empty())
        throw exceptions::Exception("FSPs empty", "calculate_four_momenta");

    if (axes.size() != squared_masses.size())
        throw exceptions::Exception("Incorrect number of masses provided ("
                                    + std::to_string(squared_masses.size()) + " != " + std::to_string(axes.size()) + ")",
                                    "calculate_four_momenta");

    // check none are negative
    if (initial_mass < 0)
        throw exceptions::Exception("initial_mass is negative", "calculate_four_momenta");
    if (std::any_of(squared_masses.begin(), squared_masses.end(), std::bind(std::less<double>(), std::placeholders::_1, 0)))
        throw exceptions::Exception("negative squared mass given", "calculate_four_momenta");

    unsigned n_fsp = FSPs.size();

    /// \todo: check sign determination on x component for particles 5 and higher
    if (n_fsp > 4)
        throw exceptions::Exception("not yet supporting 5 or more particles", "calculate_four_momenta");

    // matrix of four-vector products
    std::vector<std::vector<double> > pp(n_fsp, std::vector<double>(n_fsp, -1));

    // set diagonal elements to squared masses, and store sum
    double m2_sum_1 = 0;
    double m_sum_1 = initial_mass;
    for (size_t i = 0; i < n_fsp; ++i) {
        pp[i][i] = pow(FSPs[i]->mass(), 2);
        m2_sum_1 += pp[i][i];
        m_sum_1 -= FSPs[i]->mass();
    }

    // add two-particle invariant masses for those provided
    double m2_sum_2 = 0;
    for (size_t i = 0; i < axes.size(); ++i) {
        if (axes[i]->indices()[0] < axes[i]->indices()[1])
            pp[axes[i]->indices()[0]][axes[i]->indices()[1]] = squared_masses[i];
        else
            pp[axes[i]->indices()[1]][axes[i]->indices()[0]] = squared_masses[i];
        m2_sum_2 += squared_masses[i];
    }

    //////////////////////////////////////////////////
    // find unset two-particle invariant mass
    // and calculate fsp squared masses
    for (unsigned i = 0; i < n_fsp; ++i)
        for (unsigned j = i + 1; j < n_fsp; ++j)
            if (pp[i][j] < 0)
                // calculate unset mass: m^2_missing = M^2_isp + (n_fsp - 2) * sum_fsp (m^2) - sum_given(m^2)
                pp[i][j] = pow(initial_mass, 2) + (n_fsp - 2.) * m2_sum_1 - m2_sum_2;
    // upper triangular elements are now two-particle squared masses

    double m_01 = sqrt(pp[0][1]);

    // finish calculation of off diagonal elements
    for (unsigned i = 0; i < n_fsp; ++i)
        for (unsigned j = i + 1; j < n_fsp; ++j) {
            // first check if m_ij > (M - sum_m + m_i + m_j)
            if (pp[i][j] > pow(m_sum_1 + FSPs[i]->mass() + FSPs[j]->mass(), 2))
                return std::vector<FourVector<double> >();

            // P_i * P_j = (m^2_ij - m^2_i - m^2_j) / 2
            pp[i][j] = (pp[i][j] - pp[i][i] - pp[j][j]) / 2;
            // equal to: if m^2_ij < (m_i + m_j)^2
            if (pp[i][j] < FSPs[i]->mass() * FSPs[j]->mass())
                return std::vector<FourVector<double> >();

            pp[j][i] = pp[i][j];
        }

    //////////////////////////////////////////////////
    // calculate all four momenta in m_01 rest frame:

    std::vector<FourVector<double> > P(n_fsp, FourVector_0);

    for (unsigned i = 0; i < n_fsp; ++i) {

        FourVector<double> p;

        if (i < 2) {

            P[i][0] = (pp[0][i] + pp[1][i]) / m_01; // E
            P[i][3] = pow_negative_one(i) * P[i][0] * sqrt(1. - pp[i][i] / pow(P[i][0], 2)); // Z

            if (!std::isfinite(P[i][3]))
                return std::vector<FourVector<double> >();

        } else {

            // Energy E_i = (P_0 * P_i + P_1 * P_i) / (E_1 + E_0)
            P[i][0] = (pp[0][i] + pp[1][i]) / (P[0][0] + P[1][0]);

            // if E^2 < m^2
            if (P[i][0] < sqrt(pp[i][i]))
                return std::vector<FourVector<double> >();

            // if p0 and p1 are at rest in m01 rest frame, Z_i := 0
            // else Z_i = (P_1 * P_i - P_0 * P_i - (E_1 - E_0) * E_i) / 2 / Z_0
            P[i][3] = (P[0][3] == 0) ? 0 : (pp[1][i] - pp[0][i] - (P[1][0] - P[0][0]) * P[i][0]) / 2. / P[0][3];

            if (i < 3) {

                // p2 in y-z plane, enforce P^2 = m^2
                P[i][2] = P[i][0] * sqrt(1 - pp[i][i] / pow(P[i][0], 2) - pow(P[i][3] / P[i][0], 2));

                if (!std::isfinite(P[i][2]))
                    return std::vector<FourVector<double> >();

                // phasespace check for special case: p0 and p1 are at rest in m01 rest frame
                if (n_fsp == 3 and P[0][3] == 0 and !check_invariant_masses(axes, squared_masses, P))
                    return std::vector<FourVector<double> >();

            } else {

                // if Y_2 == 0, Y_i := 0
                // else Y_i = (P_2 * P_i - P_2 * P_i) / Y_2
                P[i][2] = (P[2][2] == 0) ? 0 : (P[2][0] * P[i][0] - pp[2][i] - P[2][3] * P[i][3]) / P[2][2];

                if (i < 4) {

                    // enforce P^2 = m^2
                    P[i][1] = P[i][0] * sqrt(1 - pp[i][i] / pow(P[i][0], 2) - pow(P[i][2] / P[i][0], 2) - pow(P[i][3] / P[i][0], 2));

                    if (!std::isfinite(P[i][1]))
                        return std::vector<FourVector<double> >();

                } else
                    throw exceptions::Exception("not yet supporting 5 or more particles", "calculate_four_momenta");
            }
        }
    }

    return P;
}

//-------------------------
std::vector<FourVector<double> > calculate_four_momenta(double initial_mass, const Model& M,
                                                        const MassAxes& axes, const std::vector<double>& squared_masses)
{
    auto P = calculate_four_momenta(initial_mass, M.finalStateParticles(), axes, squared_masses);

    // adjust for user-provided coordinate system
    const auto& C = M.coordinateSystem();
    for (auto& p : P)
        p = FourVector<double>(p[0], p[1] * C[0] + p[2] * C[1] + p[3] * C[2]);

    return P;
}


}
