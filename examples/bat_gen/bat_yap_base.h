#ifndef __BAT__BAT_YAP_BASE__H
#define __BAT__BAT_YAP_BASE__H

#include <BAT/BCModel.h>

#include <fwd/DecayingState.h>
#include <fwd/Model.h>
#include <MassAxes.h>

#include <memory>
#include <numeric>
#include <string>
#include <vector>

class bat_yap_base : public BCModel
{

public:

    bat_yap_base(std::string name, std::unique_ptr<yap::Model> M);

    virtual void MCMCUserInitialize() override
    { LikelihoodCalls_.assign(GetNChains(), 0); }

    unsigned likelihoodCalls() const
    { return std::accumulate(LikelihoodCalls_.begin(), LikelihoodCalls_.end(), 0); }

    const std::shared_ptr<yap::DecayingState> isp() const
    { return ISP_; }

    const std::unique_ptr<yap::Model>& model() const
    { return Model_; }

    const yap::MassAxes& axes() const
    { return Axes_; }

    std::shared_ptr<yap::DecayingState> isp()
    { return ISP_; }

    std::unique_ptr<yap::Model>& model()
    { return Model_; }

    yap::MassAxes& axes()
    { return Axes_; }

protected:

    void increaseLikelihoodCalls(unsigned c)
    { ++LikelihoodCalls_[c]; }

    void increaseLikelihoodCalls()
    { increaseLikelihoodCalls(GetCurrentChain()); }

private:
    yap::MassAxes Axes_;
    std::unique_ptr<yap::Model> Model_;
    std::shared_ptr<yap::DecayingState> ISP_;

    std::vector<unsigned> LikelihoodCalls_;

};

#endif
