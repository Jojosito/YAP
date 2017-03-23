#ifndef __BAT__D4PI_SCALES__H
#define __BAT__D4PI_SCALES__H

#include <AmplitudeBasis.h>
#include <MathUtilities.h>

#include <complex>

// which waves to include in the model
const bool a_rho_pi_S  = true;
const bool a_rho_pi_D  = true;
const bool a_sigma_pi  = true;
const bool rho_rho     = true;
const bool omega_omega = false;
const bool f_0_pipi    = true;
const bool f_2_pipi    = true;
const bool sigma_pipi  = true;

const bool flat_4pi    = false;

const bool a1_bowler   = true;

const bool bg_flat_4pi = false;
const bool bg_rho_2pi  = false;

// scaling to reproduce (approximately) the fit fractions of the FOCUS model
const double scale_a_rho_pi_D  = a1_bowler ? 10.5088  : 9.83252 ;
const double scale_a_sigma_pi  = a1_bowler ? 1.08726  : 1.0887  ;
const double scale_rho_rho     = a1_bowler ? 0.72479  : 0.732132;
const double scale_f_0_pipi    = a1_bowler ? 0.374681 : 0.377872;
const double scale_f_2_pipi    = a1_bowler ? 9.25244  : 9.36535 ;
const double scale_sigma_pipi  = a1_bowler ? 0.205471 : 0.2074  ;

// scaled FOCUS amplitudes
const double omega_rel_amp = -1.e-5;
const std::complex<double> amp_a_rho_pi_D   = std::polar(scale_a_rho_pi_D * 0.241, yap::rad(82.));
const std::complex<double> amp_a_omega_pi_S = omega_rel_amp;
const std::complex<double> amp_a_omega_pi_D = omega_rel_amp * amp_a_rho_pi_D;
const std::complex<double> amp_a_sigma_pi   = std::polar(scale_a_sigma_pi * 0.439, yap::rad(193.));
yap::amplitude_basis::canonical<double> amp_rho_rho(yap::amplitude_basis::transversity<double>(
        std::polar(scale_rho_rho * 0.624, yap::rad(357.)),    // longitudinal
        std::polar(scale_rho_rho * 0.157, yap::rad(120.)),    // parallel
        std::polar(scale_rho_rho * 0.384, yap::rad(163.)) )); // perpendicular
yap::amplitude_basis::canonical<double> amp_omega_omega(yap::amplitude_basis::transversity<double>(
        std::polar(omega_rel_amp * scale_rho_rho * 0.624, yap::rad(357.)),    // longitudinal
        std::polar(omega_rel_amp * scale_rho_rho * 0.157, yap::rad(120.)),    // parallel
        std::polar(omega_rel_amp * scale_rho_rho * 0.384, yap::rad(163.)) )); // perpendicular
const std::complex<double> amp_f_0_pipi     = std::polar(scale_f_0_pipi * 0.233, yap::rad(261.));
const std::complex<double> amp_f_2_pipi     = std::polar(scale_f_2_pipi * 0.338, yap::rad(317.));
const std::complex<double> amp_sigma_pipi   = std::polar(scale_sigma_pipi * 0.432, yap::rad(254.));
const std::complex<double> amp_pipiFlat     = std::polar(0.566, yap::rad(-44.9));

/* Scaled FOCUS amplitudes
D0 --> f_0, pi+, pi-, L = 0, S = 0       (mag, phase) = (0.0742105,   -99°)         (real, imag) = (-0.0116091, -0.0732968)
D0 --> f_2, pipiFlat, L = 2, S = 2       (mag, phase) = (2.66812,     -43°)         (real, imag) = (1.95134, -1.81965)
D0 --> f_0(500), pi+, pi-, L = 0, S = 0  (mag, phase) = (0.0755192,  -106°)         (real, imag) = (-0.0208159, -0.0725937)
D0 --> rho0, rho0, L = 0, S = 0          (mag, phase) = (0.27357,     162.965°)     (real, imag) = (-0.261568, 0.0801429)
D0 --> rho0, rho0, L = 2, S = 2          (mag, phase) = (0.287792,      6.38147°)   (real, imag) = (0.286009, 0.0319874)
D0 --> rho0, rho0, L = 1, S = 1          (mag, phase) = (0.236966,    163°)         (real, imag) = (-0.226611, 0.069282)
a_1+ --> rho0, pi+, L = 2, S = 1         (mag, phase) = (3.18388,      82°)         (real, imag) = (0.44311, 3.15289)
a_1+ --> f_0(500), pi+, L = 1, S = 0     (mag, phase) = (0.474726,   -167°)         (real, imag) = (-0.462559, -0.10679)
*/

// YAP fitted amplitudes
/*
const std::complex<double> amp_a_rho_pi_D   = std::polar(, yap::rad());
const std::complex<double> amp_a_omega_pi_S = std::polar(, yap::rad());
const std::complex<double> amp_a_omega_pi_D = std::polar(, yap::rad());
const std::complex<double> amp_a_sigma_pi   = std::polar(, yap::rad());
yap::amplitude_basis::canonical<double> amp_rho_rho(
        std::polar(, yap::rad()),   // S
        std::polar(, yap::rad()),   // P
        std::polar(, yap::rad()) ); // D
yap::amplitude_basis::canonical<double> amp_omega_omega(
        std::polar(, yap::rad()),   // S
        std::polar(, yap::rad()),   // P
        std::polar(, yap::rad()) ); // D
const std::complex<double> amp_f_0_pipi     = std::polar(, yap::rad());
const std::complex<double> amp_f_2_pipi     = std::polar(, yap::rad());
const std::complex<double> amp_sigma_pipi   = std::polar(, yap::rad());
const std::complex<double> amp_pipiFlat     = std::polar(, yap::rad());
*/

// YAP fitted amplitudes
/*
const std::complex<double> amp_a_rho_pi_D   = std::polar(3.18388, yap::rad(82));
const std::complex<double> amp_a_omega_pi_S = std::polar(0.01, yap::rad(180));
const std::complex<double> amp_a_omega_pi_D = std::polar(3.18388, yap::rad(82));
const std::complex<double> amp_a_sigma_pi   = std::polar(0.474726, yap::rad(-167));
yap::amplitude_basis::canonical<double> amp_rho_rho(
        std::polar(0.27357 , yap::rad(162.965)),   // S
        std::polar(0.287792, yap::rad(6.38147)),   // P
        std::polar(0.236966, yap::rad(163    )) ); // D
yap::amplitude_basis::canonical<double> amp_omega_omega(
        std::polar(0.0027357 , yap::rad(-17.0348)),   // S
        std::polar(0.00236966, yap::rad(-17     )),   // P
        std::polar(0.00287792, yap::rad(-173.619)) ); // D
const std::complex<double> amp_f_0_pipi     = std::polar(0.0742105  , yap::rad(-99 ));
const std::complex<double> amp_f_2_pipi     = std::polar(2.66812    , yap::rad(-43 ));
const std::complex<double> amp_sigma_pipi   = std::polar(0.0755192  , yap::rad(-106));
const std::complex<double> amp_pipiFlat     = std::polar(1.11022e-16, yap::rad(-90 ));
*/

#endif
