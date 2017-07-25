// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************
#include "models/a1.h"
#include "tools.h"

#include <DalitzPhspIntegral.h>
#include <Exceptions.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <make_unique.h>
#include <MassRange.h>
#include <PHSP.h>

#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <random>

int main()
{

    /*yap::plainLogs(el::Level::Info);

    const unsigned n_integrationPoints = 1e5;

    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    const double m_pi = T["pi+"].mass();
    const double m_D0 = T["D0"].mass();

    std::mt19937 g(164852419);

    const unsigned nBins = 150;
    const double low_m = 3.*m_pi;
    const double hi_m =  m_D0 - m_pi;
    //const double hi_m = sqrt(3.1);
    const double a1_mass = T["a_1+"].mass();
    const double a1_width = 0.430;

    TGraph g_phsp;
    TGraph g_density;
    TGraph g_int;
    TGraph g_w;
    TGraph g_a1_re;
    TGraph g_a1_im;
    TGraph g_a1_norm;
    TGraph g_a1_arg;
    TGraph g_a1_argand;
    TGraph g_KK;

    const unsigned n_threads = 4;

    // create model
    auto m = a1_fit();

    ImportanceSamplerGenerator impSampGen(*m.model(), n_threads);



    // Get normalizing width
    const double phsp = dalitz_phasespace_volume(a1_mass, m.model()->finalStateParticles());
    double norm_width = phsp
                        * impSampGen(a1_mass, phsp * n_integrationPoints)
                        / pow(a1_mass, 3);

    if (isnan(norm_width) or norm_width == 0)
        LOG(ERROR) << "norm_width invalid";
    LOG(INFO) << "norm_width = " << norm_width;

    for (unsigned i = 1; i <= nBins; ++i) {

        static const double scaling = 1.5; // gives higher density of samples at lower m2, where the width is changing more rapidly
        const double m2 = low_m*low_m + (hi_m*hi_m - low_m*low_m) * pow(i, scaling)/pow(nBins-1, scaling);

        double mass = sqrt(m2);

        const double phsp = dalitz_phasespace_volume(mass, m.model()->finalStateParticles());
        const int nPoints = std::max(unsigned(100), unsigned(phsp * n_integrationPoints));
        const double density = impSampGen(mass, nPoints, nPoints / 2 + 1);

        const double value = phsp * density;
        if (not isnan(value)) {
            g_phsp.SetPoint(g_int.GetN(), m2, phsp);
            g_density.SetPoint(g_int.GetN(), m2, density);
            g_int.SetPoint(g_int.GetN(), m2, value);

            double w = value / pow(mass, 3) * a1_width / norm_width;


            // K*K threshold
            double mKK2 = m2 - pow(8.9166000e-01 + 4.9367700e-01, 2);
            if (mKK2 > 0) {
                //w += 0.06 * sqrt(mKK2);
                g_KK.SetPoint(g_KK.GetN(), m2, sqrt(mKK2));
            }
            else
                g_KK.SetPoint(g_KK.GetN(), m2, 0.);

            //w = a1_width;
            g_w.SetPoint(g_w.GetN(), m2, w);

            std::complex<double> a = 1. / std::complex<double>(a1_mass * a1_mass - mass * mass, -a1_mass * w);

            g_a1_re.SetPoint(g_a1_re.GetN(), m2, real(a));
            g_a1_im.SetPoint(g_a1_re.GetN(), m2, imag(a));
            g_a1_norm.SetPoint(g_a1_norm.GetN(), m2, norm(a));
            g_a1_norm.SetPoint(g_a1_arg.GetN(), m2, arg(a));
        }

        LOG(INFO) << "Integral(" << mass << " GeV) = " << value;
    }


    TFile file("a_1_3pi_integral.root", "RECREATE");
    file.cd();
    file.WriteTObject(&g_phsp, "g_phsp_volume_vs_m2");
    file.WriteTObject(&g_density, "g_density_vs_m2");
    file.WriteTObject(&g_int, "g_integral_vs_m2");
    file.WriteTObject(&g_w, "g_width_vs_m2");
    file.WriteTObject(&g_a1_re, "g_a1_re_vs_m2");
    file.WriteTObject(&g_a1_im, "g_a1_im_vs_m2");
    file.WriteTObject(&g_a1_norm, "g_a1_norm_vs_m2");
    file.WriteTObject(&g_KK, "g_KK_vs_m2");


    file.Close();


    return 0;*/
}
