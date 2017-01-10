#include "a1MassShape.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "DataPartition.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"

namespace yap {

//-------------------------
void a1MassShape::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    // if no calculation necessary, exit
    if (D.status(*T(), si) != CalculationStatus::uncalculated)
        return;
        
    // T := 1 / (M^2 - m^2 - i * M * Gamma)
    for (auto& d : D) {
        //// 3 pion invariant mass
        //double s = model()->fourMomenta()->m2(d, pc);
        //
        //// partial width of the a_1 to 2pi+ pi-
        //auto integral_3pi = integral(Integral_)
        //// normalize: Gamma_a1_2piplus_piminus(s = Mass_) = width_
        //double norm = Width_->value() * pow(Mass_->value(), 3./2.) / integral_3pi;
        //
        //double Gamma_a1_2piplus_piminus = norm / pow(s, 3./2.) * integral_3pi;
        //
        //
        //
        //double width = 2.* Gamma_a1_2piplus_piminus + /* g^2 * K*K(s) */;
        //
        //auto M2_iMG = pow(mass()->value(), 2) - 1_i * mass()->value() * width;
        //T()->setValue(1. / (M2_iMG - s), d, si, D);
    }

    D.status(*T(), si) = CalculationStatus::calculated;
}

}




