#ifndef __BAT__BAT_GEN__H
#define __BAT__BAT_GEN__H

#include "bat_yap_base.h"

#include <fwd/DataSet.h>
#include <fwd/Model.h>

#include <memory>
#include <string>
#include <vector>

class bat_gen : public bat_yap_base
{

public:

    bat_gen(std::string name, std::unique_ptr<yap::Model> M,
            double initial_mass,
            std::vector<std::vector<unsigned> > pcs = {});

    void MCMCUserInitialize() override;

    double LogLikelihood(const std::vector<double>&) override;

    double LogAPrioriProbability(const std::vector<double>& parameters) override;

private:

    /// mass of isp
    double InitialMass_;

    std::vector<yap::DataSet> Data_;

};

#endif
