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
    // insert particles: pdg code, name, quantum numbers (charge, 2*spin, parity), mass [GeV], Breit-Widgner parameters (width) [GeV] where needed
    // initial states
    T.insert(521, "B+", QuantumNumbers(1, 0, -1), 5279.32e-3);
    T.insert(511, "B0", QuantumNumbers(0, 0, -1), 5279.63e-3);
    // final states
    T.insert(411, "D+",  QuantumNumbers(1, 0, -1), 1869.65e-3);
    T.insert(421, "D0",  QuantumNumbers(0, 0, -1), 1864.83e-3);
    T.insert(211, "pi+", QuantumNumbers(1, 0, -1), 139.57061e-3);
    T.insert(111, "pi0", QuantumNumbers(0, 0, -1), 134.9770e-2);
    // resonances
    // rho's
    QuantumNumbers rho_quantum_numbers(-1, 2, -1);
    T.insert(213,     "rho(770)+",  rho_quantum_numbers, 775.26e-3, 149.4e-3);
    T.insert(100213,  "rho(1450)+", rho_quantum_numbers, 1465e-3, 400e-3);
    T.insert(30213,   "rho(1700)+", rho_quantum_numbers, 1720e-3, 250e-3);
    // D*'s
    // TODO : Put in PDG codes
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

    // for all D* resonances add D* -> D pi, and B -> D* pi with proper charge configuration
    for (auto& Dst : {D2400, D2420, D2460}) {
        Dst->addStrongDecay(D, (B_charge == 0 ? piminus : pizero));
        B->addWeakDecay(Dst, (B_charge == 0 ? pizero : piminus));
    }

    // set amplitudes
    // You need to change these numbers
    *free_amplitude(*B, yap::to(rho770))  = std::polar(1., 0.);
    *free_amplitude(*B, yap::to(rho1450)) = std::polar(1.4, yap::rad(12.));
    *free_amplitude(*B, yap::to(rho1700)) = std::polar(1.4, yap::rad(12.));
    *free_amplitude(*B, yap::to(D2400)) = std::polar(1.4, yap::rad(12.));
    *free_amplitude(*B, yap::to(D2420)) = std::polar(1.4, yap::rad(12.));
    *free_amplitude(*B, yap::to(D2460)) = std::polar(1.4, yap::rad(12.));

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
    LOG(INFO) << "Free amplitudes (sorted by fixed/notfixed, parent_name, l):";
    for (const auto& fa : sort(free_amplitudes(M), yap::compare_by<yap::is_fixed>(),
                               yap::by_parent_name<>(), yap::by_l<>()))
        LOG(INFO) << yap::to_string(*fa);
    

    // // get default Dalitz axes
    // double B_mass = T[(B_charge == 0 ? "B0bar" : "B-")].mass();
    // auto A = M.massAxes();
    // auto m2r = yap::squared(yap::mass_range(B_mass, A, M.finalStateParticles()));

    // // phase-space generator
    // std::mt19937 g(0);
    // auto phsp_gen = std::bind(yap::phsp<std::mt19937>, std::cref(M), B_mass, A, m2r, g, std::numeric_limits<unsigned>::max());
    
    // // find maximum of model:
    // unsigned maxcalc_n_data_points = 10000;
    // auto maxcalc_data = M.createDataSet();
    // std::generate_n(std::back_inserter(maxcalc_data), maxcalc_n_data_points, phsp_gen);
    // M.calculate(maxcalc_data);
    // double max_intensity = -1;
    // for (const auto& d : maxcalc_data)
    //     max_intensity = std::max(max_intensity, yap::intensity(M, d));
    // delete maxcal_cdata;

    // // generate data
    // // auto data = M.createDataSet();
    // // unsigned n_data_points = 10000;
    // // unsigned max_number_of_tries = 1000000;
    // // for (unsigned n = 0; n < max_number_of_tries && data.size() < n_data_points; ++n) {
        
    // // }
       
    LOG(INFO) << "alright!";

    return 0;
}
