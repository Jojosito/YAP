#include "GounarisSakurai.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
#include "MathUtilities.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleCombination.h"

namespace yap {

//-------------------------
void GounarisSakurai::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    // if no calculation necessary, exit
    if (D.status(*T(), si) != CalculationStatus::uncalculated)
        return;

    /////////////////////////
    // common factors:

    const double m0 = mass()->value();
    const double m02 = m0 * m0;

    // check channel is only to two daughters
    if (pc->daughters().size() != 2)
        throw exceptions::Exception("Wrong number of daughters (" + std::to_string(pc->daughters().size()) + " != 2",
                                    "GounarisSakurai::calculate");

    /// \todo also implement for non fsp
    double s0 = 0;
    for (auto d : pc->daughters()) {
        if (not is_final_state_particle_combination(*d))
            throw exceptions::Exception("Daughter is not a final state particle",
                                        "GounarisSakurai::calculate");

        s0 += model()->finalStateParticles()[d->indices()[0]]->mass();
    }
    s0 = s0 * s0;

    const double m02ms0 = m02 - s0;

    const double Gamma0 = width()->value();

    for (auto& d : D) {
        const double s = model()->fourMomenta()->m2(d, pc);

        const double Gamma = Gamma0 * m0/sqrt(s) * pow((s - s0) / m02ms0, 3./2.);
        const double d0 = 1./pi() * m0 / sqrt(m02ms0) * (1. - s0 / m02ms0 * (2. - 3. * h(m02, s0)));
        const double M2 = m02 * (1. + Gamma0 * (s - s0) / (pi() * pow(m02ms0, 3./2.)) * (2. * h(s, s0) - (2. + s0 / (m02) * (s - m02) / (s - s0)) * h(m02, s0) - m02ms0 / (m02) * (s - m02) / (s - s0)));

        T()->setValue((m02 + d0 * m0 * Gamma0) / (M2 - s - 1_i * m0 * Gamma), d, si, D);
    }

    D.status(*T(), si) = CalculationStatus::calculated;
}

//-------------------------
double GounarisSakurai::h(double x, double s0) const
{
    return sqrt((x - s0) / x) * log((sqrt(x) + sqrt(x - s0)) / sqrt(s0));
}

//-------------------------
double GounarisSakurai::h_prime(double x, double s0) const
{
    return 1. / (2. * x) * (1. + s0 / (x - s0) * h(x, s0));
}

}




