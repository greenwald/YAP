#include "BreitWigner.h"
#include "DataSet.h"
#include "DecayChannel.h"
#include "DecayTree.h"
#include "FinalStateParticle.h"
#include "Flatte.h"
#include "FourMomenta.h"
#include "FreeAmplitude.h"
#include "HelicityAngles.h"
#include "HelicityFormalism.h"
#include "ImportanceSampler.h"
#include "logging.h"
#include "make_unique.h"
#include "MassAxes.h"
#include "MassRange.h"
#include "Model.h"
#include "ModelIntegral.h"
#include "Parameter.h"
#include "ParticleCombination.h"
#include "ParticleFactory.h"
#include "PDL.h"
#include "PoleMass.h"
#include "PHSP.h"
#include "RelativisticBreitWigner.h"
#include "Resonance.h"
#include "ZemachFormalism.h"

#include <memory>
#include <random>
#include <vector>

int main( int argc, char** argv)
{

    yap::plainLogs(el::Level::Debug);

    // use common radial size for all resonances
    double radialSize = 3.; // [GeV^-1]

    LOG(INFO) << "Start";

    yap::Model M(std::make_unique<yap::ZemachFormalism>());

    LOG(INFO) << "Model created";

    yap::ParticleFactory factory = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    LOG(INFO) << "factory created";

    // initial state particle
    auto D = factory.decayingParticle(factory.pdgCode("D+"), radialSize);

    LOG(INFO) << "D created";

    // final state particles
    auto piPlus = factory.fsp(211);
    auto piMinus = factory.fsp(-211);

    LOG(INFO) << "fsp's created";

    // set final state
    M.setFinalState(piPlus, piMinus, piPlus);

    LOG(INFO) << "final state set";

    // rho
    auto rho = factory.resonance(113, radialSize, std::make_shared<yap::RelativisticBreitWigner>());
    rho->addChannel(piPlus, piMinus);

    // f_2(1270)
    auto f_2 = factory.resonance(225, radialSize, std::make_shared<yap::BreitWigner>());
    f_2->addChannel(piPlus, piMinus);

    // f_0(980)
    auto f_0_980_flatte = std::make_shared<yap::Flatte>();
    f_0_980_flatte->addChannel(0.406, piPlus->mass()->value());
    f_0_980_flatte->addChannel(0.406 * 2, factory.particleTableEntry("K+").Mass);
    auto f_0_980 = yap::Resonance::create(yap::QuantumNumbers(0, 0), 0.965, "f_0_980", radialSize, f_0_980_flatte);
    f_0_980->addChannel(piPlus, piMinus);

    // f_0(1370)
    auto f_0_1370 = yap::Resonance::create(factory.quantumNumbers("f_0"), 1.350, "f_0_1370", radialSize, std::make_unique<yap::BreitWigner>(0.265));
    f_0_1370->addChannel(piPlus, piMinus);

    // f_0(1500)
    auto f_0_1500 = factory.resonance(factory.pdgCode("f_0(1500)"), radialSize, std::make_unique<yap::BreitWigner>());
    f_0_1500->addChannel(piPlus, piMinus);

    // sigma a.k.a. f_0(500)
    auto sigma = factory.resonance(factory.pdgCode("f_0(500)"), radialSize, std::make_unique<yap::PoleMass>(std::complex<double>(0.470, -0.220)));
    sigma->addChannel(piPlus, piMinus);

    // Add channels to D
    D->addChannel(rho,      piPlus)->freeAmplitudes().begin()->get()->setValue(std::polar(1., 0.));
    D->addChannel(f_0_980,  piPlus)->freeAmplitudes().begin()->get()->setValue(std::polar(1.4, yap::rad(12.)));
    D->addChannel(f_2,      piPlus)->freeAmplitudes().begin()->get()->setValue(std::polar(2.1, yap::rad(-123.)));
    D->addChannel(f_0_1370, piPlus)->freeAmplitudes().begin()->get()->setValue(std::polar(1.3, yap::rad(-21.)));
    D->addChannel(f_0_1500, piPlus)->freeAmplitudes().begin()->get()->setValue(std::polar(1.1, yap::rad(-44.)));
    D->addChannel(sigma,    piPlus)->freeAmplitudes().begin()->get()->setValue(std::polar(3.7, yap::rad(-3.)));
    // D->addChannel(piPlus, piMinus, piPlus);

    M.addInitialStateParticle(D);

    // check consistency
    if (M.consistent())
        LOG(INFO) << "consistent!";
    else
        LOG(INFO) << "inconsistent!";

    // print stuff
    //yap::ParticleCombination::printParticleCombinationSet();

    std::cout << "\n" << D->particleCombinations().size() << " D symmetrizations \n";

    std::cout << "\nFour momenta symmetrizations with " << M.fourMomenta()->nSymmetrizationIndices() << " indices \n";

    // std::cout << "\nHelicity angle symmetrizations with " << M.helicityAngles()->maxSymmetrizationIndex() + 1 << " indices \n";

    D->printDecayChain();
    std::cout << "\n";

    std::cout << *M.spinAmplitudeCache() << std::endl;
    M.printDataAccessors(false);

    // get default Dalitz axes
    auto A = M.massAxes();
    auto m2r = yap::squared(yap::mass_range(A, D, M.finalStateParticles()));

    // generate points randomly in phase space of model
    std::mt19937 g(0);

    // create data set
    auto data = M.createDataSet();

    // generate 10,000 phase-space-distributed data points
    std::generate_n(std::back_inserter(data), 10000,
                    std::bind(yap::phsp<std::mt19937>, std::cref(M), D->mass()->value(), A, m2r, g, std::numeric_limits<unsigned>::max()));

    LOG(INFO) << data.size() << " data points of " << data[0].bytes() << " bytes each = " << (data.size() * data[0].bytes()) * 1.e-6 << " MB";

    M.calculate(data);

    for (const auto& isp_b2 : M.initialStateParticles()) {
        for (const auto& m_dtv : isp_b2.first->decayTrees()) {

            LOG(INFO) << "\n" << to_string(m_dtv.second);

            auto A_DT = amplitude(m_dtv.second, data[0]);
            LOG(INFO) << "A_DT = " << A_DT;
            LOG(INFO) << "|A_DT|^2 = " << norm(A_DT);

            yap::ModelIntegral MI(m_dtv.second);
            yap::ImportanceSampler::calculate(MI, data);

            LOG(INFO) << "integral = " << to_string(MI.integral());
            auto ff = fit_fractions(MI);
            for (size_t i = 0; i < ff.size(); ++i)
                LOG(INFO) << "fit fraction " << 100. * ff[i] << "% for " << to_string(*MI.decayTrees()[i]->freeAmplitude());
            LOG(INFO) << "sum of fit fractions = " << std::accumulate(ff.begin(), ff.end(), 0.);

            LOG(INFO) << "cached integral components:";
            auto I_cached = cached_integrals(MI);
            for (const auto& row : I_cached)
                LOG(INFO) << std::accumulate(row.begin(), row.end(), std::string(""),
                                             [](const std::string & s, const std::complex<double>& c)
            { return s + "\t" + std::to_string(real(c)) + " + " + std::to_string(imag(c));}).erase(0, 1);

            LOG(INFO) << "integral components:";
            auto I = integrals(MI);
            for (const auto& row : I)
                LOG(INFO) << std::accumulate(row.begin(), row.end(), std::string(""),
                                             [](const std::string & s, const std::complex<double>& c)
            { return s + "\t" + std::to_string(real(c)) + " + " + std::to_string(imag(c));}).erase(0, 1);
        }
    }

    LOG(INFO) << "alright!";
}
