// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "models/a1.h"
#include "tools.h"

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

    const unsigned n_integrationPoints = 1e3;

    auto T = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/evt.pdl");

    const double m_pi = T[211].mass();
    const double m_D0 = T[421].mass();

    std::mt19937 g(0);

    const unsigned nBins = 10;
    const double low = 3.*m_pi;
    const double hi =  m_D0 - m_pi;
    const double a1_mass = 1.240;
    const double a1_width = 0.560;

    TH1D h("3piIntegral", "3piIntegral", nBins, low, hi);

    TGraph g_int;
    TGraph g_w;
    TGraph g_a1_re;
    TGraph g_a1_im;
    TGraph g_a1_norm;

    // create model
    auto m = a1_fit();

    // Get normalizing width
    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(a1_mass, m.axes(), m.model()->finalStateParticles()));
    m.integrationPointGenerator() = std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), a1_mass, m.axes(), m2r, g, std::numeric_limits<unsigned>::max());
    m.setNIntegrationPoints(n_integrationPoints, 1e5, 4);
    try {
        m.integrate();
    }
    catch (std::exception& e) {
        LOG(ERROR) << "Failed to integrate: " << e.what();
        throw;
    }
    double norm_width = integral(m.modelIntegral()).value() / pow(a1_mass, 3. / 2);
    if (isnan(norm_width) or norm_width == 0)
        LOG(ERROR) << "norm_width invalid";
    LOG(INFO) << "norm_width = " << norm_width;

    for (unsigned i = 1; i <= nBins; ++i) {

        double mass = h.GetXaxis()->GetBinCenter(i);

        // get FSP mass ranges
        m2r = yap::squared(mass_range(mass, m.axes(), m.model()->finalStateParticles()));

        m.integrationPointGenerator() = std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), mass, m.axes(), m2r, g, std::numeric_limits<unsigned>::max());
        m.setNIntegrationPoints(n_integrationPoints, 1e5, 4);

        try {
            m.integrate();
        }
        catch (std::exception& e) {
            LOG(ERROR) << "Failed to integrate: " << e.what();
            continue;
        }

        double value = integral(m.modelIntegral()).value();
        if (not isnan(value)) {
            h.SetBinContent(i, value);
            g_int.SetPoint(g_int.GetN(), mass, value);

            double w = value / pow(mass, 3./2) * a1_width / norm_width;
            g_w.SetPoint(g_w.GetN(), mass, w);

            std::complex<double> a = 1. / std::complex<double>(a1_mass * a1_mass - mass * mass, -a1_mass * w);

            g_a1_re.SetPoint(g_a1_re.GetN(), mass, real(a));
            g_a1_im.SetPoint(g_a1_re.GetN(), mass, imag(a));
            g_a1_norm.SetPoint(g_a1_norm.GetN(), mass, norm(a));
        }

        LOG(INFO) << "Integral(" << mass << ") = " << integral(m.modelIntegral()).value();
    }

    h.Draw();
    g_int.Draw("same");

    TFile file("a_1_3pi_integral.root", "RECREATE");
    file.cd();
    h.Write();
    file.WriteTObject(&g_int, "g_int");
    file.WriteTObject(&g_w, "g_w");
    file.WriteTObject(&g_a1_re, "g_a1_re");
    file.WriteTObject(&g_a1_im, "g_a1_im");
    file.WriteTObject(&g_a1_norm, "g_a1_norm");

    file.Close();


    return 0;
}
