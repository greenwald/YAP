#include <ParticleTable.h>

using namespace yap;

ParticleTable B_particle_list()
{

    ParticleTable T;

    // B's
    T.insert(ParticleTableEntry(521, "B+", QuantumNumbers(1, 0, -1), 5279.32e-3));
    T.insert(ParticleTableEntry(511, "B0", QuantumNumbers(0, 0, -1), 5279.63e-3));

             // D's
    T.insert(ParticleTableEntry(411, "D+",  QuantumNumbers(1, 0, -1), 1869.65e-3));
    T.insert(ParticleTableEntry(421, "D0",  QuantumNumbers(0, 0, -1), 1864.83e-3));

    // pi's
    T.insert(ParticleTableEntry(211, "pi+", QuantumNumbers(1, 0, -1), 139.57061e-3));
    T.insert(ParticleTableEntry(111, "pi0", QuantumNumbers(0, 0, -1), 134.9770e-3));

    // rho's
    T.insert(ParticleTableEntry(213,     "rho(770)+",  QuantumNumbers(1, 2, -1), 775.26e-3, {149.4e-3}));
    T.insert(ParticleTableEntry(100213,  "rho(1450)+", QuantumNumbers(1, 2, -1), 1465e-3, {400e-3}));
    T.insert(ParticleTableEntry(30213,   "rho(1700)+", QuantumNumbers(1, 2, -1), 1720e-3, {250e-3}));

    // D*'s
    T.insert(ParticleTableEntry(10421, "D_0*(2400)0", QuantumNumbers(0, 0, +1), 2318e-3, {267.e-3}));
    T.insert(ParticleTableEntry(10411, "D_0*(2400)+", QuantumNumbers(1, 0, +1), 2351e-3, {230.e-3}));
    T.insert(ParticleTableEntry(10423, "D_1(2420)0",  QuantumNumbers(0, 2, +1), 2420.8e-3, {31.7e-3}));
    T.insert(ParticleTableEntry(10413, "D_1(2420)+",  QuantumNumbers(1, 2, +1), 2432.2e-3, {25e-3}));
    T.insert(ParticleTableEntry(425,   "D_2*(2460)0", QuantumNumbers(0, 4, +1), 2460.7e-3, {47.5e-3}));
    T.insert(ParticleTableEntry(415,   "D_2*(2460)+", QuantumNumbers(1, 4, +1), 2465.4e-3, {46.7e-3}));

    for (const auto& t : T)
        T.insert(cp_conjugate(t));
    
    return T;
}
