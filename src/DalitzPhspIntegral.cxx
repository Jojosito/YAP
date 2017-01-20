#include "DalitzPhspIntegral.h"

#include "FastNumericalIntegration/DEIntegrator.h"
#include "FinalStateParticle.h"

namespace yap {

//-------------------------
DalitzIntegrand::DalitzIntegrand(double isp_mass, const FinalStateParticleVector& fsps) :
    M_(isp_mass)
{
    if (fsps.size() != 3)
        throw exceptions::Exception("Can only calculate Dalitz phasespace volume for 3 final state particles", "DalitzIntegrand::DalitzIntegrand");

    ma_ = fsps[0]->mass();
    mb_ = fsps[1]->mass();
    mc_ = fsps[2]->mass();
}

//-------------------------
double DalitzIntegrand::operator()(double mab2) const
{
    double mab = sqrt(mab2);
    double Eb = (mab2 - ma_*ma_ + mb_*mb_)/(2.*mab);
    double Ec = (M_*M_ - mab2 - mc_*mc_)/(2.*mab);
    double pb = sqrt(Eb*Eb - mb_*mb_);
    double pc = sqrt(Ec*Ec - mc_*mc_);
    return 4.*pb*pc;
}

//-------------------------
double dalitz_phasespace_volume(double isp_mass, const FinalStateParticleVector& fsps, const double relErr) {
    if (fsps.size() != 3)
        throw exceptions::Exception("Can only calculate Dalitz phasespace volume for 3 final state particles", "dalitz_phasespace_volume");

    double absErr = relErr;

    DalitzIntegrand f(isp_mass, fsps);
    double lowerBound = pow(fsps[0]->mass() + fsps[1]->mass(), 2);
    double upperBound = pow(isp_mass - fsps[2]->mass(), 2);

    int evaluations(0);
    double errorEstimate(1.e99); // absolute
    double result(0.);

    for (unsigned i = 0; i < 10; ++i) {
        result = DEIntegrator<DalitzIntegrand>::Integrate(f, lowerBound, upperBound, absErr, evaluations, errorEstimate);

        if (errorEstimate/result <= relErr)
            break;

        absErr = result * relErr;
    }

    return result;
}

}
