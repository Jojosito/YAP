#ifndef __BAT__D4PI_SCALES__H
#define __BAT__D4PI_SCALES__H

#include <AmplitudeBasis.h>
#include <MathUtilities.h>

#include <complex>

// which waves to include in the model
// m_D = 1.86484
// m_D - m_pi   = 1.72527
// m_D - 2 m_pi = 1.58570
const bool d4pi_omega = false;
//
const bool d4pi_a_1            = true;
//
const bool d4pi_a_1_1420       = false;
const bool d4pi_a_1_1640       = true;
//
const bool d4pi_pi1300         = false; // todo: too broad for BreitWigner???
const bool d4pi_pi1800         = true; // heavy and overlapping with other resonances
//
const bool d4pi_a_2_1320       = true; // -> rho pi // seems to have low ff in data
const bool d4pi_pi_2_1670      = true; // -> rho pi, f_2(1270) pi, pipiS pi // seems to have low ff in data
//
const bool d4pi_pipiS_pipiS    = true;
const bool d4pi_pipiS_f_0      = true; // large overlap with pipiS_pipiS and f_0_f_0
const bool d4pi_pipiS_rho      = false; // low support in data
const bool d4pi_pipiS_omega    = true and d4pi_omega;
const bool d4pi_pipiS_f_2      = true;
//
const bool d4pi_f_0_f_0        = true;
const bool d4pi_f_0_rho        = false; // low support in data
const bool d4pi_f_0_omega      = true and d4pi_omega;
const bool d4pi_f_0_f_2        = false; // very large overlap with pipiS, f2
//
const bool d4pi_rho_rho        = true;
const bool d4pi_rho_omega      = true and d4pi_omega;
const bool d4pi_rho_f_2        = false;
//
const bool d4pi_omega_omega    = true and d4pi_omega;
const bool d4pi_omega_f_2      = false and d4pi_omega; // eher kleine ff
//
const bool d4pi_f_2_f_2        = false; // up to d4pi_max_L // only S wave
//
const bool d4pi_f_0_1370_pipiS = false; // large overlap with pipiS, pipiS

//const bool d4pi_flat


// to add

// f_0(1500) into pipiS wave

// omega(1420) ??? -> rho pi
// h_1(1170) ??? -> rho pi

// all f0_1370 combinations ???

// background
const bool d4pi_bg_flat_4pi = true; // ~45% of BG
const bool d4pi_bg_rho      = false;  // ~5% of BG
const bool d4pi_bg_sigma    = false; // ~50 % of BG

const bool d4pi_bg_f_0      = false;
const bool d4pi_bg_a_1      = false;
const bool d4pi_bg_K        = false; //

/* not in bg
const bool bg_pipiS    = false;
const bool bg_omega    = false;  // 0% of BG
*/

// configuration
const bool d4pi_fix_a_rho_pi_S = true; // default: true

const bool d4pi_a1_bowler   = true and (d4pi_a_1);
const bool d4pi_pm_shared   = true;  // share free amplitudes of 3pi + and - states

const bool d4pi_a1_plus     = true;
const bool d4pi_a1_minus    = true;

const bool d4pi_fix_amplitudes  = false;
const bool d4pi_free_parameters = false; // widths, coupling constants, ...
/*
 * a1 mass, width, K*K coupling
 * f_0 mass, Flatte parameters
 * pi(1300) mass, width
 *
 */

const int  d4pi_max_L       = 2; // D
//const int  d4pi_max_L       = 3; // F

const bool d4pi_bg_only     = false; // fix D admixture to 0
const bool d4pi_fix_bg      = false; // fix bg admixtures, free D0 admixture


// YAP fitted amplitudes

/*
const std::complex<double> amp_a_rho_pi_D       = ;
const std::complex<double> amp_a_pipiS_pi       = ;
const std::complex<double> amp_a_minus_rho_pi_S = ;
const std::complex<double> amp_a_minus_rho_pi_D = ;
const std::complex<double> amp_a_minus_pipiS_pi = ;

const std::complex<double> amp_pipiS_pipiS      = ;
const std::complex<double> amp_pipiS_f_0        = ;
const std::complex<double> amp_pipiS_rho        = ;
const std::complex<double> amp_pipiS_omega      = ;
const std::complex<double> amp_pipiS_f_2        = ;

const std::complex<double> amp_f_0_f_0          = ;
const std::complex<double> amp_f_0_rho          = ;
const std::complex<double> amp_f_0_omega        = ;
const std::complex<double> amp_f_0_f_2          = ;

yap::amplitude_basis::canonical<double> amp_rho_rho(
        ,   // S
        ,   // P
         ); // D
yap::amplitude_basis::canonical<double> amp_rho_omega(
        ,   // S
        ,   // P
         ); // D
const std::vector<std::complex<double>> amp_rho_f_2 = {
        std::polar(0., yap::rad(0.)),   // S (not used)
        ,   // P
        ,   // D
         }; // F

yap::amplitude_basis::canonical<double> amp_omega_omega(
        ,   // S
        ,   // P
         ); // D
const std::vector<std::complex<double>> amp_omega_f_2 = {
        std::polar(0., yap::rad(0.)),   // S (not used)
        ,   // P
        ,   // D
         }; // F

const std::vector<std::complex<double>> amp_f_2_f_2 = {
        ,   // S
        ,   // P
        ,   // D
        ,   // F
         }; // G

const std::complex<double> amp_f_0_1370_pipiS   = ;
const std::complex<double> amp_pi1300_pi_pi_pi  = ;

const double admixture_bg_flat_4pi = ;
const double admixture_bg_rho      = ;
const double admixture_bg_omega    = ;
 */


const std::complex<double> d4pi_amp_a_rho_pi_D       = std::polar(-1.406, 1.471);
const std::complex<double> d4pi_amp_a_pipiS_pi       = std::polar(-0.115, 0.0306);
const std::complex<double> d4pi_amp_a_minus          = std::polar(1., yap::rad(0.));

const std::complex<double> d4pi_amp_pi1300_plus      = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_pi1300_minus     = std::polar(0.61 * abs(d4pi_amp_pi1300_plus), yap::rad(0.)); // from CLEO
const std::complex<double> d4pi_amp_pi1300_rho       = std::polar(1., yap::rad(0.));

const std::complex<double> d4pi_amp_pi_1800_plus     = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_pi_1800_minus    = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_pi_1800_f_0      = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_pi_1800_f_0_1500  = std::polar(1., yap::rad(0.));

const std::complex<double> d4pi_amp_a_1_1420_plus    = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_a_1_1420_minus   = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_a_1_1420_rho_S   = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_a_1_1420_rho_D   = std::polar(1., yap::rad(0.));

const std::complex<double> d4pi_amp_a_1_1640_plus    = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_a_1_1640_minus   = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_a_1_1640_rho_S   = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_a_1_1640_rho_D   = std::polar(1., yap::rad(0.));

const std::complex<double> d4pi_amp_a_2_1320_plus    = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_a_2_1320_minus   = std::polar(1., yap::rad(0.));

const std::complex<double> d4pi_amp_pi_2_1670_plus   = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_pi_2_1670_minus  = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_f_2_D            = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_rho              = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_pipiS            = std::polar(1., yap::rad(0.));

const std::complex<double> d4pi_amp_pipiS_pipiS      = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_pipiS_f_0        = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_pipiS_rho        = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_pipiS_omega      = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_pipiS_f_2        = std::polar(1., yap::rad(0.));

const std::complex<double> d4pi_amp_f_0_f_0          = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_f_0_rho          = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_f_0_omega        = std::polar(1., yap::rad(0.));
const std::complex<double> d4pi_amp_f_0_f_2          = std::polar(1., yap::rad(0.));

double d4pi_rho_rho_magnitude = 1.;
yap::amplitude_basis::canonical<double> d4pi_amp_rho_rho(
        std::polar(d4pi_rho_rho_magnitude, yap::rad(0.)),   // S              // from CLEO
        std::polar(d4pi_rho_rho_magnitude, yap::rad(0.)),   // P        // from CLEO
        std::polar(d4pi_rho_rho_magnitude, yap::rad(0.)) ); // D      // from CLEO
yap::amplitude_basis::canonical<double> d4pi_amp_rho_omega(
        std::polar(d4pi_rho_rho_magnitude, yap::rad(0.)),   // S
        std::polar(d4pi_rho_rho_magnitude, yap::rad(0.)),   // P
        std::polar(d4pi_rho_rho_magnitude, yap::rad(0.)) ); // D
const std::vector<std::complex<double>> d4pi_amp_rho_f_2 = {
        std::polar(0., yap::rad(0.)),   // S (not used)
        std::polar(1., yap::rad(0.)),   // P
        std::polar(1., yap::rad(0.)),   // D
        std::polar(0., yap::rad(0.)) }; // F

yap::amplitude_basis::canonical<double> d4pi_amp_omega_omega(
        std::polar(1., yap::rad(0.)),   // S
        std::polar(1., yap::rad(0.)),   // P
        std::polar(1., yap::rad(0.)) ); // D
const std::vector<std::complex<double>> d4pi_amp_omega_f_2 = {
        std::polar(0., yap::rad(0.)),   // S (not used)
        std::polar(1., yap::rad(0.)),   // P
        std::polar(1., yap::rad(0.)),   // D
        std::polar(1., yap::rad(0.)) }; // F

const std::vector<std::complex<double>> d4pi_amp_f_2_f_2 = {
        std::polar(1., yap::rad(0.)),   // S
        std::polar(1., yap::rad(0.)),   // P
        std::polar(1., yap::rad(0.)),   // D
        std::polar(1., yap::rad(0.)),   // F
        std::polar(1., yap::rad(0.)) }; // G

const std::complex<double> d4pi_amp_f_0_1370_pipiS     = std::polar(1., yap::rad(0.));

const double d4pi_admixture_bg_flat_4pi = 6.18;
const double d4pi_admixture_bg_rho      = 1.043/4000.;
const double d4pi_admixture_bg_sigma    = 4.431/4000.;


/*
All phases 0;
a_1 rho pi S FF ~ 40%
other a_1 FF according to FOCUS/CLEO
rho rho      FF ~ 25 %
pi 1300 ~ according to CLEO

const double scale = 0.0086;
const double scale_omega = 1.;

const std::complex<double> amp_a_rho_pi_D       = std::polar(1.33, yap::rad(0.)); // from CLEO
const std::complex<double> amp_a_pipiS_pi       = std::polar(0.22, yap::rad(0.)); // from CLEO
const std::complex<double> amp_a_minus          = std::polar(0.33, yap::rad(0.)); // from CLEO

const std::complex<double> amp_pipiS_pipiS      = std::polar(scale*1./sqrt(11.0258), yap::rad(0.));
const std::complex<double> amp_pipiS_f_0        = std::polar(scale*1./sqrt(9.12748), yap::rad(0.));
const std::complex<double> amp_pipiS_rho        = std::polar(scale*1./sqrt(0.248371), yap::rad(0.));
const std::complex<double> amp_pipiS_omega      = std::polar(scale*1./sqrt(0.0123395), yap::rad(0.));
const std::complex<double> amp_pipiS_f_2        = std::polar(scale*1./sqrt(0.000587534), yap::rad(0.));

const std::complex<double> amp_f_0_f_0          = std::polar(scale*1./sqrt(0.533048), yap::rad(0.));
const std::complex<double> amp_f_0_rho          = std::polar(scale*1./sqrt(0.0460055), yap::rad(0.));
const std::complex<double> amp_f_0_omega        = std::polar(scale*1./sqrt(0.00275747), yap::rad(0.));
const std::complex<double> amp_f_0_f_2          = std::polar(scale*1., yap::rad(0.));

double rho_rho_magnitude = 0.77;
yap::amplitude_basis::canonical<double> amp_rho_rho(
        std::polar(rho_rho_magnitude, yap::rad(0.)),   // S              // from CLEO
        std::polar(rho_rho_magnitude * 2.4, yap::rad(0.)),   // P        // from CLEO
        std::polar(rho_rho_magnitude * 11.84, yap::rad(0.)) ); // D      // from CLEO
yap::amplitude_basis::canonical<double> amp_rho_omega(
        std::polar(rho_rho_magnitude, yap::rad(0.)),   // S
        std::polar(rho_rho_magnitude * 2.4, yap::rad(0.)),   // P
        std::polar(rho_rho_magnitude * 11.84, yap::rad(0.)) ); // D
const std::vector<std::complex<double>> amp_rho_f_2 = {
        std::polar(0., yap::rad(0.)),   // S (not used)
        std::polar(scale*1./sqrt(9.71746e-08/0.0310407), yap::rad(0.)),   // P
        std::polar(scale*1./sqrt(3.02382e-08/0.0310407), yap::rad(0.)),   // D
        std::polar(0., yap::rad(0.)) }; // F

yap::amplitude_basis::canonical<double> amp_omega_omega(
        std::polar(scale_omega*rho_rho_magnitude, yap::rad(0.)),   // S
        std::polar(scale_omega*rho_rho_magnitude * 2.4, yap::rad(0.)),   // P
        std::polar(scale_omega*rho_rho_magnitude * 11.84, yap::rad(0.)) ); // D
const std::vector<std::complex<double>> amp_omega_f_2 = {
        std::polar(0., yap::rad(0.)),   // S (not used)
        std::polar(scale*1./sqrt(1.24276e-07), yap::rad(0.)),   // P
        std::polar(scale*1./sqrt(2.92872e-08), yap::rad(0.)),   // D
        std::polar(scale*1./sqrt(0.167107), yap::rad(0.)) }; // F

const std::vector<std::complex<double>> amp_f_2_f_2 = {
        std::polar(scale*1., yap::rad(0.)),   // S
        std::polar(scale*1., yap::rad(0.)),   // P
        std::polar(scale*1., yap::rad(0.)),   // D
        std::polar(scale*1., yap::rad(0.)),   // F
        std::polar(scale*1., yap::rad(0.)) }; // G

const std::complex<double> amp_f_0_1370_pipiS     = std::polar(scale*1., yap::rad(0.));

const std::complex<double> amp_pi1300_plus        = std::polar(0.035, yap::rad(0.)); // from CLEO
const std::complex<double> amp_pi1300_minus       = std::polar(0.61 * abs(amp_pi1300_plus), yap::rad(0.)); // from CLEO
const std::complex<double> amp_pi1300_plus_rho    = std::polar(14., yap::rad(0.)); // from CLEO

const double admixture_bg_flat_4pi = 0.;
const double admixture_bg_rho      = 0.;
const double admixture_bg_omega    = 0.;
 */


/*
const std::complex<double> amp_a_rho_pi_D       = std::polar(1.98035, yap::rad(-90.4228));
const std::complex<double> amp_a_pipiS_pi       = std::polar(0.09, yap::rad(0.));
const std::complex<double> amp_a_minus          = std::polar(0.115935, yap::rad(169.948));

const std::complex<double> amp_pipiS_pipiS      = std::polar(0.025386, yap::rad(139.838));
const std::complex<double> amp_pipiS_f_0        = std::polar(scale/sqrt(5.13339), yap::rad(0.));
const std::complex<double> amp_pipiS_rho        = std::polar(scale/sqrt(0.139686), yap::rad(0.));
const std::complex<double> amp_pipiS_omega      = std::polar(scale_omega*scale/sqrt(0.00693984), yap::rad(0.));
const std::complex<double> amp_pipiS_f_2        = std::polar(2.41647, yap::rad(127.021));

const std::complex<double> amp_f_0_f_0          = std::polar(0.154851, yap::rad(6.76005));
const std::complex<double> amp_f_0_rho          = std::polar(scale/sqrt(0.025874), yap::rad(0.));
const std::complex<double> amp_f_0_omega        = std::polar(scale_omega*scale/sqrt(0.00155083), yap::rad(0.));
const std::complex<double> amp_f_0_f_2          = std::polar(scale/sqrt(3.61567e-05), yap::rad(0.));

yap::amplitude_basis::canonical<double> amp_rho_rho(
        std::polar(2.68908, yap::rad(46.304)),   // S
        std::polar(0.201868, yap::rad(46.6202)),   // P
        std::polar(5.56805, yap::rad(-175.414)) ); // D
yap::amplitude_basis::canonical<double> amp_rho_omega(
        std::polar(scale_omega*scale/sqrt(0.00070011), yap::rad(0.)),   // S
        std::polar(scale_omega*scale/sqrt(0.000750721), yap::rad(0.)),   // P
        std::polar(scale_omega*scale/sqrt(8.67845e-05), yap::rad(0.)) ); // D
const std::vector<std::complex<double>> amp_rho_f_2 = {
        std::polar(0., yap::rad(0.)),   // S (not used)
        std::polar(scale/sqrt(1.76065e-06), yap::rad(0.)),   // P
        std::polar(scale/sqrt(5.47869e-07), yap::rad(0.)),   // D
        std::polar(scale/sqrt(3.29061e-08), yap::rad(0.)) }; // F

yap::amplitude_basis::canonical<double> amp_omega_omega(
        std::polar(scale_omega*scale/sqrt(2.16263e-05), yap::rad(0.)),   // S
        std::polar(scale_omega*scale/sqrt(2.41374e-05), yap::rad(0.)),   // P
        std::polar(scale_omega*scale/sqrt(3.46258e-06), yap::rad(0.)) ); // D
const std::vector<std::complex<double>> amp_omega_f_2 = {
        std::polar(0., yap::rad(0.)),   // S (not used)
        std::polar(scale_omega*scale/sqrt(6.98941e-08), yap::rad(0.)),   // P
        std::polar(scale_omega*scale/sqrt(1.64714e-08), yap::rad(0.)),   // D
        std::polar(scale_omega*scale/sqrt(1.22732e-09), yap::rad(0.)) }; // F

const std::vector<std::complex<double>> amp_f_2_f_2 = {
        std::polar(scale/sqrt(8.8749e-11), yap::rad(0.)),   // S
        std::polar(scale/sqrt(1.8876e-11), yap::rad(0.)),   // P
        std::polar(scale/sqrt(2.25569e-12), yap::rad(0.)),   // D
        std::polar(scale/sqrt(3.02117e-13), yap::rad(0.)),   // F
        std::polar(scale/sqrt(2.81809e-14), yap::rad(0.)) }; // G

const std::complex<double> amp_f_0_1370_pipiS     = std::polar(scale/sqrt(3.47548), yap::rad(0.));

const std::complex<double> amp_pi1300_plus        = std::polar(scale*0.35298920392, yap::rad(0.));
const std::complex<double> amp_pi1300_plus_rho    = std::polar(scale/sqrt(0.045107), yap::rad(0.));
const std::complex<double> amp_pi1300_minus       = std::polar(scale/sqrt(8.04123), yap::rad(0.));
 */

/*
const std::complex<double> amp_a_rho_pi_D       = std::polar(0.565685, yap::rad(-135));
const std::complex<double> amp_a_pipiS_pi       = std::polar(0.449188, yap::rad(-62.9355));
const std::complex<double> amp_a_minus_rho_pi_S = std::polar(0.226954, yap::rad(-92.3392));
const std::complex<double> amp_a_minus_rho_pi_D = std::polar(0.400695, yap::rad(86.6251));
const std::complex<double> amp_a_minus_pipiS_pi = std::polar(0.565685, yap::rad(-135));

const std::complex<double> amp_pipiS_pipiS      = std::polar(0.342089, yap::rad(-140.879));
const std::complex<double> amp_pipiS_f_0        = std::polar(0.175663, yap::rad(-16.6121));
const std::complex<double> amp_pipiS_rho        = std::polar(0.516082, yap::rad(-140.811));
const std::complex<double> amp_pipiS_omega      = std::polar(0.0139159, yap::rad(-174.283));
const std::complex<double> amp_pipiS_f_2        = std::polar(0.565685, yap::rad(-135));

const std::complex<double> amp_f_0_f_0          = std::polar(0.228112, yap::rad(-178.884));
const std::complex<double> amp_f_0_rho          = std::polar(0.145945, yap::rad(28.187));
const std::complex<double> amp_f_0_omega        = std::polar(0.0109275, yap::rad(31.0855));
const std::complex<double> amp_f_0_f_2          =  std::polar(0.565685, yap::rad(-45));

yap::amplitude_basis::canonical<double> amp_rho_rho(
        std::polar(0.542731, yap::rad(47.4776)),   // S
        std::polar(0.162253, yap::rad(71.6219)),   // P
        std::polar(0.565685, yap::rad(135)) ); // D
yap::amplitude_basis::canonical<double> amp_rho_omega(
        std::polar(0.00456069, yap::rad(148.651)),   // S
        std::polar(0.00280598, yap::rad(-170.202)),   // P
        std::polar(0.0274816, yap::rad(-162.588)) ); // D
const std::vector<std::complex<double>> amp_rho_f_2 = {
        std::polar(0., yap::rad(0.)),   // S (not used)
        std::polar(0.439789, yap::rad(114.56)),   // P
        std::polar(0.565685, yap::rad(135)),   // D
        std::polar(0.518774, yap::rad(50.4482)) }; // F

yap::amplitude_basis::canonical<double> amp_omega_omega(
        std::polar(0.000992408, yap::rad(-138.596)),   // S
        std::polar(0.000477543, yap::rad(-149.603)),   // P
        std::polar(0.00481807, yap::rad(63.758)) ); // D
const std::vector<std::complex<double>> amp_omega_f_2 = {
        std::polar(0., yap::rad(0.)),   // S (not used)
        std::polar(0.1348, yap::rad(103.006)),   // P
        std::polar(0.0885956, yap::rad(-87.6605)),   // D
        std::polar(0.361273, yap::rad(33.7057)) }; // F

const std::vector<std::complex<double>> amp_f_2_f_2 = {
        std::polar(0.565685, yap::rad(45))  ,   // S
        std::polar(0.565685, yap::rad(-45)) ,   // P
        std::polar(0.565685, yap::rad(-135)),   // D
        std::polar(0.565685, yap::rad(-45)) ,   // F
        std::polar(0.565685, yap::rad(135))  }; // G

const std::complex<double> amp_f_0_1370_pipiS   = std::polar(0., yap::rad(0.));
const std::complex<double> amp_pi1300_pi_pi_pi  = std::polar(0., yap::rad(0.));

const double admixture_bg_flat_4pi = 4.591;
const double admixture_bg_rho      = 0.1947;
const double admixture_bg_omega    = 0.00166;
*/

/*
const std::complex<double> amp_a_rho_pi_D   = std::polar(1.89514, yap::rad(-18.8152));
const std::complex<double> amp_a_pipiS_pi   = std::polar(5.69742, yap::rad(86.5461));
const std::complex<double> amp_a_minus_rho_pi_S   = std::polar(0.197753, yap::rad(-120.447));
const std::complex<double> amp_a_minus_rho_pi_D   = std::polar(1.31253, yap::rad(131.705));
const std::complex<double> amp_a_minus_pipiS_pi   = std::polar(1.44997, yap::rad(-93.5425));

yap::amplitude_basis::canonical<double> amp_rho_rho(
        std::polar(0.719214, yap::rad(121.281)),   // S
        std::polar(0.340345, yap::rad(142.623)),   // P
        std::polar(1.43405 , yap::rad(-150.721)) ); // D
yap::amplitude_basis::canonical<double> amp_omega_omega(
        std::polar(0.00325468, yap::rad(-113.889)),   // S
        std::polar(0.000948174, yap::rad(-139.279)),   // P
        std::polar(0.00709255, yap::rad(125.987)) ); // D

const std::complex<double> amp_pipiS_f_0     = std::polar(0.169688, yap::rad(-136.757));
const std::complex<double> amp_pipiS_f_2     = std::polar(7.31826, yap::rad(176.726));
const std::vector<std::complex<double>> amp_f_2_f_2 = {
        std::polar(0.1, yap::rad(0.)), // S
        std::polar(0.1, yap::rad(0.)), // P
        std::polar(0.1, yap::rad(0.)), // D
        std::polar(0.1, yap::rad(0.)) }; // F

const std::complex<double> amp_pipiS_pipiS    = std::polar(0.914355, yap::rad(67.9631));
const std::complex<double> amp_f_0_1370_pipiS = std::polar(0.1, yap::rad(0.));
const std::complex<double> amp_pipiS_rho      = std::polar(0.1, yap::rad(0.));
const std::complex<double> amp_pi1300_pi_pi_pi = std::polar(0.1, yap::rad(0.));

const double admixture_bg_flat_4pi = 13.44;
const double admixture_bg_rho      = 0.7464;
const double admixture_bg_omega    = 0.006562;
*/


/*const std::complex<double> amp_a_rho_pi_D   = std::polar(7.44092, yap::rad(-176.918));
const std::complex<double> amp_a_pipiS_pi   = std::polar(9.39864, yap::rad(177.563));
const std::complex<double> amp_a_minus_rho_pi_S   = std::polar(0.587248, yap::rad(24.4885));
const std::complex<double> amp_a_minus_rho_pi_D   = std::polar(1.61761, yap::rad(-165.684));
const std::complex<double> amp_a_minus_pipiS_pi   = std::polar(4.40066, yap::rad(174.786));

yap::amplitude_basis::canonical<double> amp_rho_rho(
        std::polar(3.38646, yap::rad(6.78347)),   // S
        std::polar(1.52758, yap::rad(9.46767)),   // P
        std::polar(0.432964, yap::rad(112.502)) ); // D
yap::amplitude_basis::canonical<double> amp_omega_omega(
        std::polar(0.00743428, yap::rad(145.577)),   // S
        std::polar(0.0105618, yap::rad(64.6185)),   // P
        std::polar(0.0369229, yap::rad(27.8847)) ); // D

const std::complex<double> amp_pipiS_f_0     = std::polar(0.870275, yap::rad(-155.194));
const std::complex<double> amp_pipiS_f_2     = std::polar(6.7082, yap::rad(176.582));
const std::vector<std::complex<double>> amp_f_2_f_2 = {
        std::polar(0.1, yap::rad(0.)), // S
        std::polar(0.1, yap::rad(0.)), // P
        std::polar(0.1, yap::rad(0.)), // D
        std::polar(0.1, yap::rad(0.)) }; // F

const std::complex<double> amp_pipiS_pipiS    = std::polar(2.8256, yap::rad(1.94176));
const std::complex<double> amp_f_0_1370_pipiS = std::polar(0.1, yap::rad(0.));
const std::complex<double> amp_pipiS_rho      = std::polar(0.1, yap::rad(0.));
const std::complex<double> amp_pi1300_pi_pi_pi = std::polar(0.1, yap::rad(0.));

const double admixture_bg_flat_4pi = 76.54;
const double admixture_bg_rho      = 7.741;
const double admixture_bg_omega    = 0.06622;*/


#endif
