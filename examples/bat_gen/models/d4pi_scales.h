#ifndef __BAT__D4PI_SCALES__H
#define __BAT__D4PI_SCALES__H

#include <AmplitudeBasis.h>
#include <MathUtilities.h>

#include <complex>

// which waves to include in the model
const bool a_rho_pi_S  = true;
const bool a_rho_pi_D  = true;
const bool a_pipiS_pi  = true;

const bool pi1300      = true;

const bool pipiS_pipiS = true;
const bool pipiS_f_0   = false;
const bool pipiS_rho   = false;
const bool pipiS_omega = false;
const bool pipiS_f_2   = true;

const bool f_0_f_0     = true;
const bool f_0_rho     = false;
const bool f_0_omega   = false;
const bool f_0_f_2     = false;

const bool rho_rho     = true;
const bool rho_omega   = false;
const bool rho_f_2     = false;

const bool omega_omega = true;
const bool omega_f_2   = false;

const bool f_2_f_2     = false;


const bool f_0_1370_pipiS = false; // lieber erstmal weglassen, highly ambiguous!


// to add
// pi1300_rho_pi
// all f0_1370 combinations
// a_1_1640
// pi_2_1670

// background
const bool bg_flat_4pi = true;
const bool bg_rho      = true;  // ~10% of BG

// configuration
const bool a1_bowler   = true;//true and (a_rho_pi_S or a_rho_pi_D);
const bool a1_shared   = false;  // share a1+ and a1- free amplitudes

const bool a1_plus     = true;
const bool a1_minus    = true;

const bool free_parameters = false; // widths, coupling constants, ...

const int max_L        = 2; // D

const bool bg_only     = false; // fix D admixture to 0


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

const double scale = sqrt(0.194167)/sqrt(60.);
const double scale_omega = 0.1;

const std::complex<double> amp_a_rho_pi_D       = std::polar(1.98035, yap::rad(-90.4228));
const std::complex<double> amp_a_pipiS_pi       = std::polar(0.0334687, yap::rad(0.));
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

const double admixture_bg_flat_4pi = 0.;
const double admixture_bg_rho      = 0.;
const double admixture_bg_omega    = 0.;

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
