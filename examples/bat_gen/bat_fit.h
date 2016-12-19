#ifndef __BAT__BAT_FIT__H
#define __BAT__BAT_FIT__H

#include "fit_base.h"

#include <fwd/FreeAmplitude.h>
#include <fwd/Parameter.h>

#include <memory>
#include <string>
#include <vector>

class bat_fit : public fit_base
{

public:

    /// constructor
    /// \param name name of bat model
    /// \param M yap::model
    /// \param pcs vector<vector<unsigned>> defining mass axes
    bat_fit(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs = {});

    /// add a complex parameter from the model to bat
    void addParameter(std::string name, std::shared_ptr<yap::ComplexParameter> P, std::complex<double> low, std::complex<double> high, std::string latex = "", std::string units = "");

    /// ada a real parameter from the model to bat
    void addParameter(std::string name, std::shared_ptr<yap::RealParameter> P, double low, double high, std::string latex = "", std::string units = "");

    /// set the priors for a FreeAmplitude's amplitude and phase
    void setPriors(std::shared_ptr<yap::FreeAmplitude> fa, BCPrior* amp_prior, BCPrior* arg_prior);
    
    /// set the range for a FreeAmplitude's real and imaginary parts
    void setRealImagRanges(std::shared_ptr<yap::FreeAmplitude> fa, double real_low, double real_high, double imag_low, double imag_high);
    
    /// set the range for a FreeAmplitude's abs and arg obvservables
    void setAbsArgRanges(std::shared_ptr<yap::FreeAmplitude> fa, double abs_low, double abs_high, double arg_low, double arg_high);

    /// fix a FreeAmplitude
    void fix(std::shared_ptr<yap::FreeAmplitude> A, double amp, double phase);

    /// log prior
    virtual double LogAPrioriProbability(const std::vector<double>& p) override;

    /// calculate  observables
    virtual void CalculateObservables(const std::vector<double>& p) override;

    /// init size of CalculatedFitFractions_
    virtual void MCMCUserInitialize() override;

    /// set parameters into model
    virtual void setParameters(const std::vector<double>& p) override;

    /// find the position in the parameter list of the first element of a free amplitude
    size_t findFreeAmplitude(std::shared_ptr<yap::FreeAmplitude> A) const;

    /// \return free amplitudes
    const yap::FreeAmplitudeVector& freeAmplitudes() const
    { return FreeAmplitudes_; }

protected:

    /// vector of parameters to set in model
    yap::ParameterVector Parameters_;

    /// offset of where first user-set parameter is
    int FirstParameter_{-1};

    /// offset of where first user-set observable is
    int FirstObservable_{-1};

    /// list of decay trees integrated over
    yap::DecayTreeVector DecayTrees_;

    /// Free amplitudes of model to set
    yap::FreeAmplitudeVector FreeAmplitudes_;

    /// BCPrior on abs(FreeAmplitude)
    std::vector<BCPrior*> AbsPriors_;

    /// BCPrior on arg(FreeAmplitude)
    std::vector<BCPrior*> ArgPriors_;

    /// Calculated fit fractions (for observables)
    std::vector<yap::RealIntegralElementVector> CalculatedFitFractions_;

};

#endif
