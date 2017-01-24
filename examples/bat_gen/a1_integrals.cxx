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
#include <MassRange.h>
#include <PHSP.h>
#include <logging.h>
#include <make_unique.h>

#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TTree.h>

#include <algorithm>
#include <chrono>
#include <random>

int main()
{
    yap::plainLogs(el::Level::Info);

    const unsigned n_integrationPoints = 1e4;

    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    const double m_pi = T["pi+"].mass();
    const double m_D0 = T["D0"].mass();

    std::mt19937 g(164852419);

    const unsigned nBins = 150;
    const double low_m = 3.*m_pi;
    const double hi_m =  m_D0 - m_pi;
    //const double hi_m = 3;
    const double a1_mass = T["a_1+"].mass();
    const double a1_width = 0.560;

    TH1D h("3piIntegral", "3piIntegral", nBins, low_m*low_m, hi_m*hi_m);

    TGraph g_phsp;
    TGraph g_density;
    TGraph g_int;
    TGraph g_w;
    TGraph g_a1_re;
    TGraph g_a1_im;
    TGraph g_a1_norm;

    const unsigned n_threads = 4;

    // create model
    auto m = a1_fit();

    // Get normalizing width
    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(a1_mass, m.axes(), m.model()->finalStateParticles()));
    m.integrationPointGenerator() = std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), a1_mass, m.axes(), m2r, g, std::numeric_limits<unsigned>::max());
    m.setNIntegrationPoints(n_integrationPoints, 1e5, n_threads);
    m.integrate();

    double norm_width = dalitz_phasespace_volume(a1_mass, m.model()->finalStateParticles()) * integral(m.modelIntegral()).value() / pow(a1_mass, 3. / 2);
    if (isnan(norm_width) or norm_width == 0)
        LOG(ERROR) << "norm_width invalid";
    LOG(INFO) << "norm_width = " << norm_width;

    for (unsigned i = 1; i <= nBins; ++i) {

        double m2 = h.GetXaxis()->GetBinCenter(i);
        double mass = sqrt(m2);

        // get FSP mass ranges
        m2r = yap::squared(mass_range(mass, m.axes(), m.model()->finalStateParticles()));

        m.integrationPointGenerator() = std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), mass, m.axes(), m2r, g, std::numeric_limits<unsigned>::max());
        m.setNIntegrationPoints(n_integrationPoints, 1e5, n_threads);
        m.integrate();

        const double phsp = dalitz_phasespace_volume(mass, m.model()->finalStateParticles());
        const double density = integral(m.modelIntegral()).value();
        const double value = phsp * density;
        if (not isnan(value)) {
            h.SetBinContent(i, value);
            g_phsp.SetPoint(g_int.GetN(), m2, phsp);
            g_density.SetPoint(g_int.GetN(), m2, density);
            g_int.SetPoint(g_int.GetN(), m2, value);

            double w = value / pow(mass, 3./2) * a1_width / norm_width;
            //w = a1_width;
            g_w.SetPoint(g_w.GetN(), m2, w);

            std::complex<double> a = 1. / std::complex<double>(a1_mass * a1_mass - mass * mass, -a1_mass * w);

            g_a1_re.SetPoint(g_a1_re.GetN(), m2, real(a));
            g_a1_im.SetPoint(g_a1_re.GetN(), m2, imag(a));
            g_a1_norm.SetPoint(g_a1_norm.GetN(), m2, norm(a));
        }

        LOG(INFO) << "Integral(" << mass << " GeV) = " << value;
    }

    h.Draw();
    g_int.Draw("same");

    TFile file("a_1_3pi_integral.root", "RECREATE");
    file.cd();
    h.Write();
    file.WriteTObject(&g_phsp, "g_phsp_volume_vs_m2");
    file.WriteTObject(&g_density, "g_density_vs_m2");
    file.WriteTObject(&g_int, "g_integral_vs_m2");
    file.WriteTObject(&g_w, "g_width_vs_m2");
    file.WriteTObject(&g_a1_re, "g_a1_re_vs_m2");
    file.WriteTObject(&g_a1_im, "g_a1_im_vs_m2");
    file.WriteTObject(&g_a1_norm, "g_a1_norm_vs_m2");

    file.Close();


    return 0;
}
