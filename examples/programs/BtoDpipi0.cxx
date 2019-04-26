#include <Attributes.h>
#include <BreitWigner.h>
#include <ConstantWidthBreitWigner.h>
#include <container_utils.h>
#include <DataSet.h>
#include <DecayChannel.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <DecayTreeVectorIntegral.h>
#include <FinalStateParticle.h>
#include <Flatte.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <FreeAmplitude.h>
#include <HelicityAngles.h>
#include <HelicityFormalism.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <make_unique.h>
#include <MassAxes.h>
#include <MassRange.h>
#include <Model.h>
#include <ModelIntegral.h>
#include <Parameter.h>
#include <ParticleCombination.h>
#include <ParticleTable.h>
#include <PDL.h>
#include <PHSP.h>
#include <PoleMass.h>
#include <Sort.h>
#include <ZemachFormalism.h>

#include <memory>
#include <random>
#include <vector>

int main( int argc, char** argv)
{

    yap::plainLogs(el::Level::Debug);

    //-------------------------
    // create a particle table
    yap::ParticleTable T;
    // insert particles: pdg code, name, quantum numbers (charge, 2*spin, parity), mass, Breit-Widgner parameters (width) where needed
    // initial states
    T.insert(-521, "B-",    QuantumNumbers(-1, 0, -1), 5279.32e-3);
    T.insert(-511, "B0bar", QuantumNumbers( 0, 0, -1), 5279.63e-3);
    // final states
    T.insert( 411, "D+",    QuantumNumbers(+1, 0, -1), 1869.65e-3);
    T.insert( 421, "D0",    QuantumNumbers( 0, 0, -1), 1864.83e-3);
    T.insert(-211, "pi-",   QuantumNumbers(-1, 0, -1), 139.57061e-3);
    T.insert( 111, "pi0",   QuantumNumbers( 0, 0, -1), 134.9770e-2);
    // resonances
    // rho's
    QuantumNumbers rho_quantum_numbers(-1, 2, -1);
    T.insert(-213,     "rho770-",  rho_quantum_numbers, 775.26e-3, 149.4e-3);
    T.insert(-100213,  "rho1450-", rho_quantum_numbers, 1465e-3, 400e-3);
    T.insert(-30213,   "rho1700-", rho_quantum_numbers, 1720e-3, 250e-3);
    // D's
    T.insert(, "D*2400z", QuantumNumbers( 0, 0, +1), 2318e-3, 267e-3);
    T.insert(, "D*2400-", QuantumNumbers(-1, 0, +1), 2351e-3, 230e-3);
    T.insert(, "D2420z",  QuantumNumbers( 0, 2, +1), 2420.8e-3, 31.7e-3);
    T.insert(, "D2420-",  QuantumNumbers(-1, 2, +1), 2432.2e-3, 25e-3);
    T.insert(, "D*2460z", QuantumNumbers( 0, 4, +1), 2460.7e-3, 47.5e-3);
    T.insert(, "D*2460-", QuantumNumbers(-1, 4, +1), 2465.4e-3, 46.7e-3);
    //-------------------------

    //-------------------------
    // make model of decay
    //-------------------------

    // choose which B decays
    int B_charge = -1;
    
    // use common radial size for all resonances
    double r = 3.; // [GeV^-1]

    // initial state
    auto B = yap::DecayingParticle::create(T[(B_charge == 0 ? "B0bar" : "B-")], r);
    
    // final state
    auto D = yap::FinalStateParticle::create(T[(B_charge == 0 ? "D+" : "D0")]);
    auto piminus = yap::FinalStateParticle::create(T["pi-"]);
    auto pi0     = yap::FinalStateParticle::create(T["pi0"]);

    //-------------------------
    // intermediate resonances
    //-------------------------

    // rho's
    auto rho770  = yap::DecayingParticle::create(T["rho770-"],  r, std::make_shared<yap::BreitWigner>(T["rho770-"]));
    auto rho1450 = yap::DecayingParticle::create(T["rho1450-"], r, std::make_shared<yap::BreitWigner>(T["rho1450-"]));
    auto rho1700 = yap::DecayingParticle::create(T["rho1700-"], r, std::make_shared<yap::BreitWigner>(T["rho1700-"]));

    auto D2400 = yap::DecayingParticle::create(T[(B_charge == 0 ? "D*2400z" : "D*2400-")], r, std::make_shared<yap::BreitWigner>(T[(B_charge ? : "D*2400z", "D*2400-")]));
    auto D2420 = yap::DecayingParticle::create(T[(B_charge == 0 ? "D2420z" : "D2420-")], r, std::make_shared<yap::BreitWigner>(T[(B_charge ? : "D2420z", "D2420-")]));
    auto D2460 = yap::DecayingParticle::create(T[(B_charge == 0 ? "D*2460z" : "D*2460-")], r, std::make_shared<yap::BreitWigner>(T[(B_charge ? : "D*2460z", "D*2460-")]));
    
    // create model and set final state
    yap::Model M(std::make_unique<yap::ZemachFormalism>());
    // yap::Model M(std::make_unique<yap::HelicityFormalism>());
    M.setFinalState(D, piminus, pizero);

    // for all rho's add rho -> pi- pi0; and B -> D rho
    for (auto& rho : {rho770, rho1700, rho1700}) {
        rho->addStrongDecay(piminus, pizero);
        B->addWeakDecay(D, rho);
    }

    // for all D resonances add D -> pi pi, and B -> D* pi with proper charge configuration
    for (auto& Dst : {D2400, D2420, D2460}) {
        Dst->addStrongDecay(D, (B_charge == 0 ? piminus : pizero));
        B->addWeakDecay(Dst, (B_charge == 0 ? pizero : piminus));
    }
    
    // nonresonant decay
    D->addWeakDecay(D, piminus, pizero);

    *free_amplitude(*D, yap::to(rho))                     = std::polar(1., 0.);
    *free_amplitude(*D, yap::to(f_0_980))                 = std::polar(1.4, yap::rad(12.));
    *free_amplitude(*D, yap::to(f_2))                     = std::polar(2.1, yap::rad(-123.));
    *free_amplitude(*D, yap::to(f_0_1370))                = std::polar(1.3, yap::rad(-21.));
    *free_amplitude(*D, yap::to(f_0_1500))                = std::polar(1.1, yap::rad(-44.));
    *free_amplitude(*D, yap::to(sigma))                   = std::polar(3.7, yap::rad(-3.));
    *free_amplitude(*D, yap::to(piPlus, piMinus, piPlus)) = std::polar(0.1, yap::rad(45.));

    M.lock();

    // check consistency
    if (!M.consistent()) {
        LOG(INFO) << "inconsistent!";
        return 1;
    }

    // print stuff

    LOG(INFO) << "Decays:";
    MULTILINE(LOG(INFO),to_decay_string(*D));

    LOG(INFO);
    LOG(INFO) << "SpinAmplitudeCache:";
    MULTILINE(LOG(INFO),to_string(*M.spinAmplitudeCache()));

    LOG(INFO);
    LOG(INFO) << "Free amplitudes (sorted by fixed/notfixed, parent_name, l):";
    for (const auto& fa : sort(free_amplitudes(M), yap::compare_by<yap::is_fixed>(),
                               yap::by_parent_name<>(), yap::by_l<>()))
        LOG(INFO) << yap::to_string(*fa);
    
    // get default Dalitz axes
    auto A = M.massAxes();
    auto m2r = yap::squared(yap::mass_range(T["D+"].mass(), A, M.finalStateParticles()));

    // generate points randomly in phase space of model
    std::mt19937 g(0);

    // create data set
    auto data = M.createDataSet();

    // generate 10,000 phase-space-distributed data points
    std::generate_n(std::back_inserter(data), 10000,
                    std::bind(yap::phsp<std::mt19937>, std::cref(M), T["D+"].mass(), A, m2r, g, std::numeric_limits<unsigned>::max()));

    LOG(INFO);
    LOG(INFO) << data.size() << " data points of " << data[0].bytes() << " bytes each = " << data.bytes() * 1.e-6 << " MB";

    M.calculate(data);

    yap::ModelIntegral MI(M);
    yap::ImportanceSampler::calculate(MI, data);

    for (const auto& mci : MI.integrals()) {

        LOG(INFO);
        // MULTILINE(LOG(INFO), to_string(mci.Integral.decayTrees()));

        auto A_DT = amplitude(mci.Integral.decayTrees(), data[0]);
        LOG(INFO) << "A_DT = " << A_DT;
        LOG(INFO) << "|A_DT|^2 = " << norm(A_DT);

        LOG(INFO) << "integral = " << to_string(integral(mci.Integral));
        auto ff = fit_fractions(mci.Integral);
        for (size_t i = 0; i < ff.size(); ++i)
            LOG(INFO) << "fit fraction " << to_string(100. * ff[i]) << "% for " << to_string(*mci.Integral.decayTrees()[i]->freeAmplitude());
        LOG(INFO) << "sum of fit fractions = " << to_string(std::accumulate(ff.begin(), ff.end(), yap::RealIntegralElement()));

        // LOG(INFO) << "cached integral components:";
        // auto I_cached = cached_integrals(mci.Integral);
        // for (const auto& row : I_cached)
        //     LOG(INFO) << std::accumulate(row.begin(), row.end(), std::string(""),
        //                                  [](std::string & s, const yap::ComplexIntegralElement & c)
        //                                  { return s += "\t" + to_string(c);}).erase(0, 1);

        LOG(INFO) << "integral components:";
        auto I = integrals(mci.Integral);
        for (const auto& row : I)
            LOG(INFO) << std::accumulate(row.begin(), row.end(), std::string(""),
                                         [](std::string & s, const yap::ComplexIntegralElement & c)
                                         { return s += "\t" + to_string(c);}).erase(0, 1);
    }

    // LOG(INFO) << std::endl << "Fixed amplitudes: ";
    // for (const auto& fa : free_amplitudes(M, yap::is_fixed()))
    //     LOG(INFO) << yap::to_string(*fa);

    // LOG(INFO) << std::endl << "Free amplitudes: ";
    // for (const auto& fa : free_amplitudes(M, yap::is_not_fixed()))
    //     LOG(INFO) << yap::to_string(*fa);

    // for (unsigned l = 0; l <= 2; ++l) {
    //     LOG(INFO) << std::endl << "Amplitudes with l = " << l << ": ";
    //     for (const auto& fa : free_amplitudes(M, yap::l_equals(l)))
    //         LOG(INFO) << yap::to_string(*fa);
    // }

    LOG(INFO) << "alright!";

    return 0;
}
