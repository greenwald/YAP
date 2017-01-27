#ifndef __BAT__MI_FIT__H
#define __BAT__MI_FIT__H

#include "fit_base.h"

#include <fwd/Model.h>

#include <memory>
#include <string>
#include <vector>

class TTree;

class mi_fit : public fit_base
{

public:

    /// constructor
    /// \param name name of bat model
    /// \param M yap::model
    /// \param pcs vector<vector<unsigned>> defining mass axes
    mi_fit(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs = {});

    /// set parameters into model
    virtual void setParameters(const std::vector<double>& p) override;

    void printSwave(const std::string& filename, const std::vector<double>& pars = {}, const std::vector<double>& uncs = {});

    std::vector<double> truth() const;
    
    virtual double LogAPrioriProbability(const std::vector<double>& P) override;
    
protected:

    std::shared_ptr<yap::FreeAmplitude> FAtoPiPiSwave_;
    std::shared_ptr<yap::FreeAmplitude> FAtoF2_;

    yap::FreeAmplitudeVector FAfromPiPiSwave_;
    
};

#endif
