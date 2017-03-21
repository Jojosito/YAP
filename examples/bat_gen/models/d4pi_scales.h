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
const bool omega_omega = true;
const bool f_0_pipi    = true;
const bool f_2_pipi    = true;
const bool sigma_pipi  = true;

const bool flat_4pi    = true;

const bool a1_bowler   = true;

const bool bg_flat_4pi = false;
const bool bg_rho_2pi  = false;

// scaling to reproduce (approximately) the fit fractions of the FOCUS model
const double scale_a_rho_pi_D  = a1_bowler ? 13.2111  : 9.83252 ;
const double scale_a_sigma_pi  = a1_bowler ? 1.08138  : 1.0887  ;
const double scale_rho_rho     = a1_bowler ? 0.617098 : 0.732132;
const double scale_f_0_pipi    = a1_bowler ? 0.3185   : 0.377872;
const double scale_f_2_pipi    = a1_bowler ? 7.89385  : 9.36535 ;
const double scale_sigma_pipi  = a1_bowler ? 0.174813 : 0.2074  ;

// scaled FOCUS amplitudes
const double omega_rel_amp = -0.01;
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
const std::complex<double> amp_pipiFlat     = 0.;


// YAP fitted amplitudes
/*const std::complex<double> amp_a_rho_pi_D   = std::polar(4.24525, yap::rad(-4.78052));
const std::complex<double> amp_a_omega_pi_S = std::polar(0.0553893, yap::rad(166.538));
const std::complex<double> amp_a_omega_pi_D = std::polar(0.166574, yap::rad(121.108));
const std::complex<double> amp_a_sigma_pi   = std::polar(0.776635, yap::rad(-27.7332));
yap::amplitude_basis::canonical<double> amp_rho_rho(
        std::polar(1.06742, yap::rad(-85.108)),   // S
        std::polar(16.1126, yap::rad(174.119)),   // P
        std::polar(3.46506, yap::rad(-90.1618)) ); // D
yap::amplitude_basis::canonical<double> amp_omega_omega(
        std::polar(0.084335, yap::rad(47.2535 )),   // S
        std::polar(0.343503, yap::rad(-37.8345)),   // P
        std::polar(0.238281, yap::rad(47.2729 )) ); // D
const std::complex<double> amp_f_0_pipi     = std::polar(0.245876, yap::rad(41.8381));
const std::complex<double> amp_f_2_pipi     = std::polar(3.32129, yap::rad(136.129));
const std::complex<double> amp_sigma_pipi   = std::polar(0.545923, yap::rad(111.971));
*/


/*
a_1+ --> rho0, pi+, L = 0, S = 1  =  (1, 0 deg)
a_1+ --> rho0, pi+, L = 2, S = 1  =  (4.24525, -4.78052 deg)
a_1+ --> omega, pi+, L = 0, S = 1  =  (0.0553893, 166.538 deg)
a_1+ --> omega, pi+, L = 2, S = 1  =  (0.166574, 121.108 deg)
a_1+ --> f_0(500), pi+, L = 1, S = 0  =  (0.776635, -27.7332 deg)

D0 --> rho0, rho0, L = 0, S = 0  =  (1.06742, -85.108 deg)
D0 --> rho0, rho0, L = 1, S = 1  =  (16.1126, 174.119 deg)
D0 --> rho0, rho0, L = 2, S = 2  =  (3.46506, -90.1618 deg)

D0 --> rho0, rho0, L = 0, S = 0  =  (1.48856, -94.4085 deg)
D0 --> rho0, rho0, L = 1, S = 1  =  (1.83628e-09, 112.138 deg)
D0 --> rho0, rho0, L = 2, S = 2  =  (4.00409, -120.074 deg)

D0 --> omega, omega, L = 0, S = 0  =  (0.084335, 47.2535 deg)
D0 --> omega, omega, L = 1, S = 1  =  (0.343503, -37.8345 deg)
D0 --> omega, omega, L = 2, S = 2  =  (0.238281, 47.2729 deg)

D0 --> f_0, pi+, pi-, L = 0, S = 0  =  (0.245876, 41.8381 deg)
D0 --> f_2, pipiFlat, L = 2, S = 2  =  (3.32129, 136.129 deg)
D0 --> f_0(500), pi+, pi-, L = 0, S = 0  =  (0.545923, 111.971 deg)
D0 --> pipiFlat, pipiFlat, L = 0, S = 0  =  (3.57214, 6.42922 deg)
*/










#endif
