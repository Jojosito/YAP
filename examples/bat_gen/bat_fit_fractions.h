#ifndef __BAT__BAT_FIT_FRACTIONS__H
#define __BAT__BAT_FIT_FRACTIONS__H

#include "bat_fit.h"

class bat_fit_fractions : public bat_fit
{

public:

    /// constructor
    /// \param name name of bat model
    /// \param M yap::model
    /// \param pcs vector<vector<unsigned>> defining mass axes
    bat_fit_fractions(std::string name, std::unique_ptr<yap::Model> M, const std::vector<std::vector<unsigned> >& pcs = {}) :
        bat_fit(name, std::move(M), pcs)
    { }

    /// log likelihood
    double LogLikelihood(const std::vector<double>& p) override
    {
        bat_fit::setParameters(p);
        return 0.;
    }

};

#endif
