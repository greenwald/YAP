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
#include <logging.h>
#include <HelicityFormalism.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <Model.h>
#include <Parameter.h>
#include <ParticleFactory.h>
#include <PDL.h>
#include <Resonance.h>
#include <ZemachFormalism.h>

#include <cmath>

/**
 * Test that the amplitude remains the same after swapping the axes of the DalitzPlot
 */

TEST_CASE( "swapDalitzAxes" )
{

    // disable debug logs in test
    yap::disableLogs(el::Level::Debug);
    //yap::plainLogs(el::Level::Debug);

    std::array<double, 4> PDGs = {411, 321, -321, 211}; // D+ -> K+ K- pi+
    std::array<double, 2> m2_ab_range = {0.4, 1.9};
    std::array<double, 2> m2_bc_range = {0.9, 3.1};

    auto F = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");


    const unsigned N = 20;
    for (double m2_ab = m2_ab_range[0]; m2_ab <= m2_ab_range[1]; m2_ab += (m2_ab_range[1] - m2_ab_range[0]) / N) {
        for (double m2_bc = m2_bc_range[0]; m2_bc <= m2_bc_range[1]; m2_bc += (m2_bc_range[1] - m2_bc_range[0]) / N) {

            // calc 3rd inv mass square
            double m2_ac = pow(F[PDGs[0]].mass(), 2)
                + pow(F[PDGs[1]].mass(), 2)
                + pow(F[PDGs[2]].mass(), 2)
                + pow(F[PDGs[3]].mass(), 2)
                - (m2_ab + m2_bc);

            if (m2_ac < 0.) {
                //std::cout << "m2_ac < 0.\n";
                continue;
            }


            // loop over SpinFormalisms
            for (unsigned i_formalism = 0; i_formalism < 2; ++i_formalism) {

                std::vector<std::complex<double> > resultingAmplitudes(6, 0.);

                // loop over axis swaps
                for (unsigned i = 0; i < 6; ++i) {

                    auto M = (i_formalism == 0)
                        ? dkkp<yap::ZemachFormalism>(PDGs[0], std::vector<int>(PDGs.begin() + 1, PDGs.end()))
                        : dkkp<yap::HelicityFormalism>(PDGs[0], std::vector<int>(PDGs.begin() + 1, PDGs.end()));

                    // Dalitz coordinates
                    yap::MassAxes A;
                    std::vector<double> squared_masses;
                    switch (i) {
                        case 0:
                            // original
                            A = M->massAxes({{0, 1}, {1, 2}});
                            squared_masses = {m2_ab, m2_bc};
                            break;
                        case 1:
                            // 0 <-> 1
                            A = M->massAxes({{1, 0}, {0, 2}});
                            squared_masses = {m2_ab, m2_ac};
                            break;
                        case 2:
                            // 0 <-> 2
                            A = M->massAxes({{2, 1}, {1, 0}});
                            squared_masses = {m2_bc, m2_ab};
                            break;
                        case 3:
                            // 1 <-> 2
                            A = M->massAxes({{0, 2}, {2, 1}});
                            squared_masses = {m2_ac, m2_bc};
                            break;
                        case 4:
                            A = M->massAxes({{1, 2}, {2, 0}});
                            squared_masses = {m2_bc, m2_ac};
                            break;
                        case 5:
                            A = M->massAxes({{2, 0}, {0, 1}});
                            squared_masses = {m2_ac, m2_ab};
                            break;
                    }

                    // calculate four-momenta
                    auto P = calculate_four_momenta(F[PDGs[0]].mass(), *M, A, squared_masses);

                    // if failed, outside phase space
                    if (P.empty()) {
                        //std::cout<<"PhSp ";
                        continue;
                    }

                    // reset data set
                    auto data = M->createDataSet();
                    // add point
                    data.push_back(P);

                    M->calculate(data);
                    resultingAmplitudes[i] = amplitude(M->initialStates().begin()->first->decayTrees().at(0), data[0]);
                    //std::cout<<resultingAmplitudes[i]<<"   ";
                }

                REQUIRE( std::isfinite(resultingAmplitudes[0].real()) );
                REQUIRE( std::isfinite(resultingAmplitudes[0].imag()) );
                for (unsigned i = 1; i < resultingAmplitudes.size(); ++i)
                    REQUIRE( resultingAmplitudes[0] == Catch::Detail::CApprox(resultingAmplitudes[i]) );

                //std::cout<<"\n";

            }// end loop over SpinFormalisms

        }
    }




}
