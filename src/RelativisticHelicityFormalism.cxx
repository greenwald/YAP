#include "RelativisticHelicityFormalism.h"

#include "CachedValue.h"
#include "ClebschGordan.h"
#include "HelicityAngles.h"
#include "logging.h"
#include "Model.h"
#include "Spin.h"
#include "WignerD.h"

namespace yap {

//-------------------------
RelativisticHelicitySpinAmplitude::RelativisticHelicitySpinAmplitude(Model& m, unsigned two_J, const SpinVector& two_j, unsigned l, unsigned two_s) :
    SpinAmplitude(m, two_J, two_j, l, two_s, equal_by_shared_pointer),
    RequiresHelicityAngles(two_J != 0)
{
    if (finalTwoJ().size() != 2)
        throw exceptions::Exception("Wrong number of daughter spins specified (" + std::to_string(finalTwoJ().size()) + " != 2)",
                                    "RelativisticHelicitySpinAmplitude::RelativisticHelicitySpinAmplitude");

    // check j1j2S triangle
    if (!triangle(finalTwoJ()[0], finalTwoJ()[1], twoS()))
        throw exceptions::AngularMomentumNotConserved("RelativisticHelicitySpinAmplitude::RelativisticHelicitySpinAmplitude");

    // angular momentum normalization factor
    /// \todo check which is the right one
    // double c = sqrt(2. * L() + 1);
    // double c = sqrt((2. * L() + 1) / 4. / pi<double>() );
    double c  = sqrt((2. * L() + 1) / (initialTwoJ() + 1.));

    auto finalJSum = std::accumulate(finalTwoJ().begin(), finalTwoJ().end(), 0u) / 2;
    auto MaxPsiInternal = *std::min_element(finalTwoJ().begin(), finalTwoJ().end()) / 2;
    auto MaxPsiPhi = std::min(initialTwoJ(), finalTwoJSum) / 2;
    auto MaxPsiChi = std::min(L(), finalJSum);
    auto MaxChiPhi = std::min(L(), initialTwoJ() / 2);
    auto totalRank = finalJSum + L() + initialTwoJ() / 2;
    auto MaxDelta = std::min(initialTwoJ(), twoS()) / 2;
    auto evenContraction = is_even(L() + finalJSum - initialTwoJ() / 2);
    auto IndexContractions = evenContraction ? totalRank / 2 : (totalRank - 3) / 2;

    struct ContractionKey
    {
        unsigned PsiInternal{0};
        unsigned cChiPhi{0};
        unsigned cPsiChi{0};
        unsigned cPsiPhi{0};
        unsigned cChiOmega{0};
        unsigned cPhiOmega{0};
        unsigned cPhiEps{0};
        unsigned cChiEps{0};

        unsigned ContractionNumber{0};
        
        const std::array<unsigned, 8> as_array() const
            { return {PsiInternal, cChiPhi cPsiChi, cPsiPhi, cChiOmega, cPhiOmega, cPhiEps, cChiEps}; }

        
        const bool operator<(const ContractionKey& rhs)
            {
                auto arr = as_array();
                auto rhs_arr = rhs.as_array();
                return lexicographical_compare(arr.begin(), arr.end(), rhs_arr.begin(), rhs_arr.end());
            }
    };

    std::map<ContractionKey, SOMETHING> contractions;
    
    ContractionKey K;
    
    unsigned MaxContractionNumber = 0;
    for (K.PsiInternal = 0; K.PsiInternal <= MaxPsiInternal; ++K.PsiInternal) {
        for (K.cChiPhi = 0; K.cChiPhi <= MaxChiPhi; ++K.cChiPhi) {
            for (K.cPsiChi = 0; K.cPsiChi <= MaxPsiChi; ++K.cPsiChi) {
                for (K.cPsiPhi = 0; K.cPsiPhi <= MaxPsiPhi; ++K.cPsiPhi) {

                    /////////////////////////
                    // check selection rules
                    if (K.PsiInternal + K.cPsiChi + K.cChiPhi + K.cPsiPhi != IndexContractions)
                        continue;

                    if (even_contraction) {
                        if (2 * K.PsiInternal + K.cPsiChi + K.cPsiPhi != finalJSum
                            or K.cPsiPhi + K.cChiPhi != L() or K.cChiPhi + K.cPsiPhi != initialTwoJ() / 2)
                            continue;
                    } else if (L() - K.cPsiPhi - K.cChiPhi > 1
                               or initialTwoJ() / 2- K.cPsiPhi - K.cChiPhi > 1
                               or finalJSum - 2 * K.PsiInternal - K.cPsiChi - K.cPsiChi < 0
                               or finalJSum - 2 * K.PsiInternal - K.cPsiChi - K.cPsiChi > 2)
                        continue;
                    /////////////////////////

                    /////////////////////////
                    // get contraction rank
                    std::vector<unsigned> rank;
                    rank.reserve(finalTwoJ().size());
                    std::transform(finalTwoJ().begin(), finalTwoJ().end(), std::back_inserter(rank),
                                   [&](SpinVector::value_type& twoj){return twoj / 2 - K.PsiInternal;});
                    
                    if (!evenContraction and finalJSum - 2 * K.PsiInternal - K.cPsiChi - K.cPsiChi == 2) {
                        // if any already zero, continue
                        if (std::any_of(rank.begin(), rank.end(), [](unsigned r){return r == 0;}))
                            continue;
                        // decrement each rank
                        std::for_each(rank.begin(), rank.end() [](unsigned& r){--r;});
                    }
                    /////////////////////////

                    for (K.cChiOmega = 0; K.cChiOmega <= rank[0] and K.cChiOmega <= K.cPsiChi; ++K.cChiOmega) {
                        for (K.cPhiOmega = 0; K.cPhiOmega <= rank[0] - K.cChiOmega and K.cPhiOmega <= K.cPsiPhi; ++K.cPhiOmega) {
                            K.cPhiEps = K.cPsiPhi - K.cPhiOmega;
                            K.cChiEps = K.cPsiChi - K.cChiOmega;
                            if (K.cPhiEps + K.cChiEps > rank[1])
                                continue;

                            for (unsigned delta = 0; delta <= MaxDelta; ++delta) {
                                auto it = contractions.find(K);
                                unsigned ContractionNumber = (it != contractions.end()) it->ContractionNumber : MaxContractionNumber + 1;
                                
                            }
                        }
                    }
                    
                }
            }
        }
    }
    
    // cache coefficients for each spin projection state
    for (const auto& two_m : projections(finalTwoJ()))



        try {
            double CG = ClebschGordan::couple(finalTwoJ()[0], two_m[0], finalTwoJ()[1], two_m[1], L(), twoS(), initialTwoJ());

            if (CG == 0)
                continue;

            Coefficients_[two_m] = c * CG;

            // add amplitudes for all initial spin projections
            for (auto two_M : projections(initialTwoJ()))
                // for J==0, only the Clebsch-Gordan coefficients are returned
                // ==> we don't need storage space in the DataPoint
                addAmplitude(two_M, two_m, initialTwoJ() == 0);

        } catch (const exceptions::InconsistentSpinProjection&) { /* ignore */ }

    if (Coefficients_.empty())
        throw exceptions::Exception("no valid nonzero Clebsch-Gordan coefficients stored", "RelativisticHelicitySpinAmplitude::RelativisticHelicitySpinAmplitude");
}

//-------------------------
const std::complex<double> RelativisticHelicitySpinAmplitude::calc(int two_M, const SpinProjectionVector& two_m,
        const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
{
    // helicity angles
    double phi   = model()->helicityAngles()->phi(d, pc);
    double theta = model()->helicityAngles()->theta(d, pc);

    std::vector<double> gamma;
    gamma.reserve(pc->daughters().size());
    std::transform(pc->daughters().begin(), pc->daughters().end(), std::back_inserter(gamma),
                   [&](const ParticleCombinationVector::value_type& daughter)
                   {return gamma(model()->fourMomenta()->p(d, daughter));});
    
    return std::conj(DFunction(initialTwoJ(), two_M, two_m[0] - two_m[1], phi, theta, 0))
           * Coefficients_.at(two_m);

    /// \todo Take a look at momentum-dependent Clebsch-Gordan
    /// coefficients by J. Friedrich and S.U. Chung implemented in
    /// rootPWA by C. Bicker
}

}
