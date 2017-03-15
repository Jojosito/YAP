// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D4PI_SCALES__H
#define __BAT__D4PI_SCALES__H


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
const double scale_a_rho_pi_D  = a1_bowler ? 14.0539  : 4.2395231085;
const double scale_a_sigma_pi  = a1_bowler ? 1.09234  : 1.0739853519;
const double scale_rho_rho     = a1_bowler ? 0.668474 : 0.7158375483;
const double scale_f_0_pipi    = a1_bowler ? 0.345016 : 0.3762733548;
const double scale_f_2_pipi    = a1_bowler ? 8.55105  : 12.687283302;
const double scale_sigma_pipi  = a1_bowler ? 0.189367 : 0.204794164 ;

#endif
