#include "PiPiSWave.h"

#include <CalculationStatus.h>
#include <DataPartition.h>
#include <FourMomenta.h>
#include <logging.h>
#include <Matrix.h>
#include <MeasuredBreakupMomenta.h>
#include <Model.h>
#include <PDL.h>

#include <assert.h>

namespace yap {

//-------------------------
PiPiSWaveAuMorganPennington::PiPiSWaveAuMorganPennington() :
    a_(2, SquareMatrix<std::complex<double>, 2>(0)),
    c_(5, SquareMatrix<std::complex<double>, 2>(0)),
    vesSheet_(false),
    CachedMassShape_(ComplexCachedValue::create(*this))
{
    const double f[2] = {0.1968, -0.0154};  // AMP Table 1, M solution: f_1^1 and f_2^1

    a_[0][0][0] =  0.1131;  // AMP Table 1, M solution: f_2^2
    a_[0][0][1] =  0.0150;  // AMP Table 1, M solution: f_1^3
    a_[0][1][0] =  0.0150;  // AMP Table 1, M solution: f_1^3
    a_[0][1][1] = -0.3216;  // AMP Table 1, M solution: f_2^3
    a_[1][0][0] = f[0] * f[0];
    a_[1][0][1] = f[0] * f[1];
    a_[1][1][0] = f[1] * f[0];
    a_[1][1][1] = f[1] * f[1];

    c_[0][0][0] =  0.0337;                // AMP Table 1, M solution: c_11^0
    c_[1][0][0] = -0.3185;                // AMP Table 1, M solution: c_11^1
    c_[2][0][0] = -0.0942;                // AMP Table 1, M solution: c_11^2
    c_[3][0][0] = -0.5927;                // AMP Table 1, M solution: c_11^3
    c_[4][0][0] =  0.1957;                // AMP Table 1, M solution: c_11^4
    c_[0][0][1] = c_[0][1][0] = -0.2826;  // AMP Table 1, M solution: c_12^0
    c_[1][0][1] = c_[1][1][0] =  0.0918;  // AMP Table 1, M solution: c_12^1
    c_[2][0][1] = c_[2][1][0] =  0.1669;  // AMP Table 1, M solution: c_12^2
    c_[3][0][1] = c_[3][1][0] = -0.2082;  // AMP Table 1, M solution: c_12^3
    c_[4][0][1] = c_[4][1][0] = -0.1386;  // AMP Table 1, M solution: c_12^4
    c_[0][1][1] =  0.3010;                // AMP Table 1, M solution: c_22^0
    c_[1][1][1] = -0.5140;                // AMP Table 1, M solution: c_22^1
    c_[2][1][1] =  0.1176;                // AMP Table 1, M solution: c_22^2
    c_[3][1][1] =  0.5204;                // AMP Table 1, M solution: c_22^3
    c_[4][1][1] = -0.3977;                // AMP Table 1, M solution: c_22^4

    sP_[0][0] = -0.0074;  // AMP Table 1, M solution: s_0
    sP_[0][1] =  0.9828;  // AMP Table 1, M solution: s_1

    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    piChargedMass_   = T["pi+"].mass();
    piNeutralMass_   = T["pi0"].mass();
    kaonChargedMass_ = T["K+" ].mass();
    kaonNeutralMass_ = T["K0" ].mass();
    kaonMeanMass_    = (kaonChargedMass_ + kaonNeutralMass_) / 2.;

    assert(piChargedMass_ > 0);
    assert(piNeutralMass_ > 0);
    assert(kaonChargedMass_ > 0);
    assert(kaonNeutralMass_ > 0);
    assert(kaonMeanMass_ > 0);
}

//-------------------------
void PiPiSWaveAuMorganPennington::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    const std::complex<double> imag(0., 1.);

    for (auto& d : D) {

        double mass = model()->fourMomenta()->m(d, pc);
        double s    = mass * mass;
        if (fabs(s - sP_[0][1]) < 1.e-6) {
            mass += 1.e-6;
            s     = mass * mass;
        }

        const std::complex<double> qPiPi   = measured_breakup_momenta::q_complex(s, piChargedMass_,   piChargedMass_  );
        const std::complex<double> qPi0Pi0 = measured_breakup_momenta::q_complex(s, piNeutralMass_,   piNeutralMass_  );
        const std::complex<double> qKK     = measured_breakup_momenta::q_complex(s, kaonChargedMass_, kaonChargedMass_);
        const std::complex<double> qK0K0   = measured_breakup_momenta::q_complex(s, kaonNeutralMass_, kaonNeutralMass_);
        std::complex<double>       qKmKm   = measured_breakup_momenta::q_complex(s, kaonMeanMass_,    kaonMeanMass_   );

        SquareMatrix<std::complex<double>, 2> rho(0);
        if (vesSheet_) {
            if (qKmKm.imag() > 0)
                qKmKm *= -1;
            rho[0][0] = (2. * qPiPi) / mass;
            rho[1][1] = (2. * qKmKm) / mass;
        } else {
            rho[0][0] = ((2. * qPiPi) / mass + (2. * qPi0Pi0) / mass) / 2.;
            rho[1][1] = ((2. * qKK)   / mass + (2. * qK0K0)   / mass) / 2.;
        }

        const double scale = (s / (4. * kaonMeanMass_ * kaonMeanMass_)) - 1.;

        SquareMatrix<std::complex<double>, 2> M(0);
        for (unsigned int i = 0; i < 2 /*sP_.size2()*/; ++i) {
            const std::complex<double> fa = 1. / (s - sP_[0][i]);
            M += fa * a_[i];
        }
        for (unsigned int i = 0; i < c_.size(); ++i) {
            const std::complex<double> sc = pow(scale, (int)i);
            M += sc *c_[i];
        }

        // modification: off-diagonal terms set to 0
        M[0][1] = 0;
        M[1][0] = 0;

        const SquareMatrix<std::complex<double>, 2> matrix = M - imag * rho;
        const std::complex<double> det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        CachedMassShape_->setValue(matrix[1][1] / det, d, si, D);
    }

    D.status(*CachedMassShape_, si) = CalculationStatus::calculated;
}

//-------------------------
const std::complex<double> PiPiSWaveAuMorganPennington::value(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
{
    return CachedMassShape_->value(d, symmetrizationIndex(pc));
}

//-------------------------
PiPiSWaveAuMorganPenningtonKachaev::PiPiSWaveAuMorganPenningtonKachaev() :
    PiPiSWaveAuMorganPennington()
{
    // change parameters according to Kachaev's prescription
    c_[4][0][0] = 0; // was 0.1957;
    c_[4][1][1] = 0; // was -0.3977;

    a_[0][0][1] = 0; // was 0.0150
    a_[0][1][0] = 0; // was 0.0150

    // a_[1] are the f's from the AMP paper
    a_[1][0][0] = 0;
    a_[1][0][1] = 0;
    a_[1][1][0] = 0;
    a_[1][1][1] = 0;
}

}
