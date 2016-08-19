// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_gen.h"

#include <CalculationStatus.h>
#include <DataSet.h>
#include <DecayingParticle.h>
#include <DecayTree.h>
#include <FinalStateParticle.h>
#include <FourMomenta.h>
#include <FreeAmplitude.h>
#include <MassRange.h>
#include <Model.h>
#include <ParticleCombination.h>
#include <SpinAmplitudeCache.h>
#include <VariableStatus.h>

#include <complex>
#include <limits>

// -----------------------
bat_gen::bat_gen(std::string name, std::unique_ptr<yap::Model> M, double initial_mass, std::vector<std::vector<unsigned> > pcs)
: bat_yap_base(name, std::move(M)), InitialMass_(initial_mass)
{
    for (auto& kv : model()->initialStateParticles()) {

        std::cout << "Initial state particle " << to_string(*kv.first) << " with " << to_string(kv.second) << std::endl;

        for (const auto& m_dtv : kv.first->decayTrees())
            for (const auto& dt : m_dtv.second)
                if (dt->freeAmplitude()->variableStatus() != yap::VariableStatus::fixed)
                    std::cout << to_string(*dt->freeAmplitude())
                              << "  =  (" << abs(dt->freeAmplitude()->value()) << ", "
                              << yap::deg(arg(dt->freeAmplitude()->value())) << " deg)"
                              << std::endl;
        std::cout << std::endl;
    }

    axes() = model()->massAxes(pcs);
    auto m2r = yap::squared(yap::mass_range(InitialMass_, axes(), model()->finalStateParticles()));

    for (size_t i = 0; i < axes().size(); ++i) {
        std::string axis_label = "m2_" + indices_string(*axes()[i]).substr(1, 2);
        AddParameter(axis_label, m2r[i][0], m2r[i][1], axis_label, "[GeV]");
        std::cout << "Added parameter " << axis_label
                  << " with range = [" << m2r[i][0] << ", " << m2r[i][1] << "]"
                  << std::endl;
    }
}

// ---------------------------------------------------------
void bat_gen::MCMCUserInitialize()
{
    bat_yap_base::MCMCUserInitialize();
    Data_.assign(GetNChains(), model()->createDataSet(1));
}

// ---------------------------------------------------------
double bat_gen::LogLikelihood(const std::vector<double>&)
{
    unsigned c = GetCurrentChain();
    double L = sum_of_log_intensity(*model(), Data_[c]);
    // model()->setParameterFlagsToUnchanged();
    increaseLikelihoodCalls(c);
    return L;
}

// ---------------------------------------------------------
double bat_gen::LogAPrioriProbability(const std::vector<double>& parameters)
{
    // calculate four-momenta
    auto P = calculate_four_momenta(InitialMass_, *model(), axes(), parameters);

    // if failed, outside phase space
    if (P.empty())
        return -std::numeric_limits<double>::infinity();

    unsigned c = GetCurrentChain();
    Data_[c].setAll(yap::VariableStatus::changed);
    Data_[c].setAll(yap::CalculationStatus::uncalculated);
    model()->setFinalStateMomenta(Data_[c][0], P, Data_[c]);
    return 0;
}
