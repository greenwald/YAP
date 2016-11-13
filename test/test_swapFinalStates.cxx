#include <catch.hpp>
#include <catch_capprox.hpp>

#include "helperFunctions.h"

#include <Attributes.h>
#include <BreitWigner.h>
#include <container_utils.h>
#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayTree.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <MassRange.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <Resonance.h>
#include <ZemachFormalism.h>

#include <cmath>

/**
 * Test that the amplitude remains the same after swapping the order of the final state particles
 */
std::complex<double> calculate_model(double isp_mass, yap::Model& M, const yap::MassAxes& A, std::vector<double> m2, yap::DataSet& data)
{
    auto ISPs = full_final_state_isp(M);
    if (ISPs.empty())
        throw;

    auto isp = ISPs[0];

    // calculate four-momenta
    auto P = calculate_four_momenta(isp_mass, M, A, m2);

    // if failed, outside phase space
    if (P.empty())
        return std::numeric_limits<double>::quiet_NaN();

    // reset data set
    data = M.createDataSet();
    // add point
    data.push_back(P);

    // return amplitude
    M.calculate(data);
    return amplitude(isp->decayTrees().at(0), data[0]);
}

TEST_CASE( "swapFinalStates" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    auto D_mass = F["D+"].mass();

    // create models
    std::vector<std::shared_ptr<yap::Model> > Z;     // Zemach
    std::vector<yap::MassAxes> mZ; // always (pi+, K-), (K-, K+)
    std::vector<yap::DataSet> dZ;

    std::vector<std::shared_ptr<yap::Model> > H;     // Helicity
    std::vector<yap::MassAxes> mH; // always (pi+, K-), (K-, K+)
    std::vector<yap::DataSet> dH;

    std::vector<int> FSP = {F["K-"].pdg(), F["pi+"].pdg(), F["K+"].pdg()};
    std::sort(FSP.begin(), FSP.end());
    do {

        unsigned i_piPlus = std::distance(FSP.begin(), std::find(FSP.begin(), FSP.end(), 211));
        unsigned i_kPlus  = std::distance(FSP.begin(), std::find(FSP.begin(), FSP.end(), 321));
        unsigned i_kMinus = std::distance(FSP.begin(), std::find(FSP.begin(), FSP.end(), -321));
        
        // Zemach
        Z.push_back(dkkp<yap::ZemachFormalism>(411, FSP));
        mZ.push_back(Z.back()->massAxes({{i_piPlus, i_kMinus}, {i_kMinus, i_kPlus}}));
        dZ.push_back(Z.back()->createDataSet(1));

        // Helicity
        H.push_back(dkkp<yap::HelicityFormalism>(411, FSP));
        mH.push_back(H.back()->massAxes({{i_piPlus, i_kMinus}, {i_kMinus, i_kPlus}}));
        dH.push_back(H.back()->createDataSet(1));

    } while (std::next_permutation(FSP.begin(), FSP.end()));

    // get piK and KK mass ranges
    auto m2_piK_range = yap::squared(yap::mass_range(D_mass, mZ[0][0], Z[0]->finalStateParticles()));
    auto m2_KK_range  = yap::squared(yap::mass_range(D_mass, mZ[0][1], Z[0]->finalStateParticles()));

    const unsigned N = 20;
    // loop over phase space
    for (double m2_piK = m2_piK_range[0]; m2_piK <= m2_piK_range[1]; m2_piK += (m2_piK_range[1] - m2_piK_range[0]) / N) {
        for (double m2_KK = m2_KK_range[0]; m2_KK <= m2_KK_range[1]; m2_KK += (m2_KK_range[1] - m2_KK_range[0]) / N) {

            std::vector<std::complex<double> > amps_Z(Z.size(), 0.);
            std::vector<std::complex<double> > amps_H(H.size(), 0.);

            for (size_t i = 0; i < Z.size(); ++i) {
                FDEBUG("Calculate Zemach for FinalState combination " << i);
                amps_Z[i] = calculate_model(D_mass, *Z[i], mZ[i], {m2_piK, m2_KK}, dZ[i]);

                FDEBUG("Calculate Helicity for FinalState combination " << i);
                amps_H[i] = calculate_model(D_mass, *H[i], mH[i], {m2_piK, m2_KK}, dH[i]);
            }



            // print
            DEBUG(m2_piK << ", " << m2_KK << " is " << ((std::isnan(real(amps_Z[0]))) ? "out" : "in") << " phase space");

            DEBUG("Zemach:                        Helicity:");
            //for (size_t i = 0; i < amps_Z.size(); ++i)
            for (size_t i = 0; i < 1; ++i) {
                double phaseDiff = arg(amps_Z[i]) - arg(amps_H[i]);
                DEBUG(amps_Z[i] << " " << norm(amps_Z[i]) << "     " << amps_H[i] << " " << norm(amps_H[i])
                      << "      ratio Z/H = " <<  norm(amps_Z[i]) / norm(amps_H[i])
                      << "      rel. phase = " << yap::deg(phaseDiff) << "°");
            }

            // if nan, check that all are nan
            if (amps_Z[0] != amps_Z[0]) {
                for (size_t i = 1; i < amps_Z.size(); ++i)
                    REQUIRE ( amps_Z[i] != amps_Z[i] );
                for (size_t i = 0; i < amps_H.size(); ++i)
                    REQUIRE ( amps_H[i] != amps_H[i] );
            }
            // otherwise check that amplitudes are approximately the same
            else {
                // // check equality for Zemach
                for (size_t i = 1; i < amps_Z.size(); ++i)
                    REQUIRE ( amps_Z[i - 1] == Catch::Detail::CApprox( amps_Z[i] ) );

                // check equality for Helicity
                for (size_t i = 1; i < amps_H.size(); ++i)
                    REQUIRE ( amps_H[i - 1] == Catch::Detail::CApprox( amps_H[i] ) );

                // todo check if Zemach and Helicity have the same phase
                double phaseDiff = yap::deg(arg(amps_Z[0]) - arg(amps_H[0]));
                /*// 180 degree are ok????????????
                if (phaseDiff >= 180)
                    phaseDiff -= 180;
                if (phaseDiff <= -180)
                    phaseDiff += 180;*/
                //REQUIRE( phaseDiff == Approx(0) );

                /// \todo check if Zemach and Helicity have same magnitude
                //REQUIRE( norm(amps_Z[i]) == Approx(norm(amps_H[i])) );
            }

        }
    }

}
