#include "SmearedFlatte.h"

#include "CalculationStatus.h"
#include "DataPoint.h"
#include "DataPartition.h"
#include "DecayChannel.h"
#include "Exceptions.h"
#include "FinalStateParticle.h"
#include "FourMomenta.h"
#include "logging.h"
#include "MeasuredBreakupMomenta.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleTable.h"

namespace yap {

//-------------------------
void SmearedFlatte::calculate(DataPartition& D, const std::shared_ptr<const ParticleCombination>& pc, unsigned si) const
{
    // if no calculation necessary, exit
    if (D.status(*T(), si) != CalculationStatus::uncalculated)
        return;

    /////////////////////////
    // common factors:

    // mass^2
    const double m  = mass()->value();
    const double m2 = m*m;

    // width at nominal mass
    std::complex<double> w = 0;
    for (const auto& fc : channels())
        w += fc.Coupling->value() * measured_breakup_momenta::q_complex(m2, fc.Particles[0]->mass(), fc.Particles[1]->mass());

    // normalization factor
    auto w_o_m = 2. * w / m;


    // calculate amplitude for masses
    double m_min = std::max(0., m_min_ - 4. * sigma_);
    double m_max = m_max_ + 4. * sigma_;
    double m_step = 0.01 * sigma_;
    std::vector<double> shape_real;
    std::vector<double> shape_imag;
    shape_real.reserve( (m_max - m_min) / m_step );
    shape_imag.reserve(shape_real.size());

    for (double sqrt_s = m_min; sqrt_s < m_max; sqrt_s += m_step) {
        double s = sqrt_s * sqrt_s;

        // calculate width term := sum of coupling * complex-breakup-momentum
        std::complex<double> ws = 0;
        for (const auto& fc : channels())
            ws += fc.Coupling->value() * measured_breakup_momenta::q_complex(s, fc.Particles[0]->mass(), fc.Particles[1]->mass());

        auto amp = w_o_m / (m2 - s - 1_i * 2. * ws / m);
        shape_real.push_back(amp.real());
        shape_imag.push_back(amp.imag());
    }

    // smear
    gaussianiir1d(shape_real.data(), shape_real.size(), sigma_/m_step, 10);
    gaussianiir1d(shape_imag.data(), shape_imag.size(), sigma_/m_step, 10);

    /////////////////////////
    for (auto& d : D) {
        double sqrt_s = model()->fourMomenta()->m(d, pc);
        std::complex<double> amp(interpolate(sqrt_s, shape_real, m_min, m_step), interpolate(sqrt_s, shape_imag, m_min, m_step));
        T()->setValue(amp, d, si, D);
    }

    D.status(*T(), si) = CalculationStatus::calculated;
}

double SmearedFlatte::interpolate(double m, std::vector<double>& shape, double m_min, double m_step) const
{
    const double m_max = m_min + shape.size() * m_step;
    if (m < m_min || m > m_max)
        throw exceptions::Exception("Out of boundaries", "SmearedFlatte::interpolate");

    unsigned i = (m - m_min) / m_step;

    if (i >= shape.size())
        return shape.back();

    double mL = m_min + i * m_step;
    double mR = mL + m_step;
    double yL = shape.at(i);
    double yR = shape.at(i+1);

    double dydx = ( yR - yL ) / ( mR - mL ); // gradient

    return yL + dydx * ( m - mL ); // linear interpolation
}


/**
 * \file gaussianiir1d.c
 * \brief Fast 1D Gaussian convolution IIR approximation
 * \author Pascal Getreuer <getreuer@gmail.com>
 *
 * Copyright (c) 2011, Pascal Getreuer
 * All rights reserved.
 *
 * This program is free software: you can redistribute it and/or modify it
 * under, at your option, the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version, or the terms of the
 * simplified BSD license.
 *
 * You should have received a copy of these licenses along with this program.
 * If not, see <http://www.gnu.org/licenses/> and
 * <http://www.opensource.org/licenses/bsd-license.html>.
 */

/**
 * \brief Fast 1D Gaussian convolution IIR approximation
 * \param data the data to be convolved, modified in-place
 * \param length number of elements
 * \param sigma the standard deviation of the Gaussian in pixels
 * \param numsteps number of timesteps, more steps implies better accuracy
 *
 * Implements the fast Gaussian convolution algorithm of Alvarez and Mazorra,
 * where the Gaussian is approximated by a cascade of first-order infinite
 * impulsive response (IIR) filters.  Boundaries are handled with half-sample
 * symmetric extension.
 *
 * Gaussian convolution is approached as approximating the heat equation and
 * each timestep is performed with an efficient recursive computation.  Using
 * more steps yields a more accurate approximation of the Gaussian.  A
 * reasonable default value for \c numsteps is 4.
 *
 * Reference:
 * Alvarez, Mazorra, "Signal and Image Restoration using Shock Filters and
 * Anisotropic Diffusion," SIAM J. on Numerical Analysis, vol. 31, no. 2,
 * pp. 590-605, 1994.
 */
void SmearedFlatte::gaussianiir1d(double* data, long length, double sigma, unsigned numsteps) const
{

    if(!data || length < 1 || sigma <= 0 || numsteps < 0)
        return;

    const double lambda = (sigma*sigma)/(2.0*numsteps);
    const double nu = (1.0 + 2.0*lambda - sqrt(1.0 + 4.0*lambda))/(2.0*lambda);
    const double boundaryscale = 1./(1. - nu);
    const double postscale = pow(nu/lambda, numsteps);

    for(unsigned step = 0; step < numsteps; step++)
    {
        long i;

        data[0] *= boundaryscale;

        /* Filter rightwards (causal) */
        for(i = 1; i < length; i++)
            data[i] += nu * data[i - 1];

        data[i = length - 1] *= boundaryscale;

        /* Filter leftwards (anti-causal) */
        for(; i > 0; i--)
            data[i - 1] += nu*data[i];
    }

    for(long i = 0; i < length; i++)
        data[i] *= postscale;

    return;
}

}
