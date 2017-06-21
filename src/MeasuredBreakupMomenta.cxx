#include "MeasuredBreakupMomenta.h"

#include "DataPoint.h"
#include "Exceptions.h"
#include "FourMomenta.h"
#include "Model.h"
#include "ParticleCombination.h"

namespace yap {

namespace measured_breakup_momenta {

//-------------------------
double q2(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc, const Model& m)
{
    if (pc->daughters().size() != 2)
        throw exceptions::Exception("invalid number of daughters (" + std::to_string(pc->daughters().size()) + ")",
                                    "measured_breakup_momenta::q2");

    return q2(m.fourMomenta()->m2(d, pc),
        m.fourMomenta()->m(d, pc->daughters()[0]),
        m.fourMomenta()->m(d, pc->daughters()[1]));
}

//-------------------------
double q2(double m2_R, double m_a, double m_b)
{
    return (m_a == m_b)
            ?
                    m2_R / 4. - m_a * m_a
                    :
                    (m2_R - pow(m_a + m_b, 2)) * (m2_R - pow(m_a - m_b, 2)) / m2_R / 4.;
}

//-------------------------
std::complex<double> q_complex(double m2_R, double m_a, double m_b)
{
    const double q_2 = q2(m2_R, m_a, m_b) ;
    const double q  = sqrt(fabs(q_2));
    if (q_2 < 0.)
        return std::complex<double>(0., q);
    return std::complex<double>(q, 0.);
}

}

}
