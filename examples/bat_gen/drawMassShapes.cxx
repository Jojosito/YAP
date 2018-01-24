/*
 * drawMassShapes.cxx
 *
 *  Created on: Jul 17, 2017
 *      Author: ne53mad
 */

#include "models/d4pi.h"

#include <DataSet.h>
#include <MassAxes.h>
#include <FourMomenta.h>
#include <FourVector.h>
#include <logging.h>

#include <TFile.h>
#include <TGraph.h>
#include <TH1D.h>
#include <TTree.h>

int main()
{
    // generate model
    auto M = d4pi();

    auto T = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    //assert(m.model()->initialStates()[0]->name() == "D0");
    auto isp_mass = T["D0"].mass();

    // generate data
    auto A = M->massAxes();
    auto m2r = yap::squared(yap::mass_range(isp_mass, A, M->finalStateParticles()));

    std::mt19937 g(0);
    // fill data set with nPoints points
    yap::DataSet data(*M);
    std::generate_n(std::back_inserter(data), 10000,
            std::bind(yap::phsp<std::mt19937>, std::cref(*M), isp_mass, A, m2r, g, std::numeric_limits<unsigned>::max()));

    auto partitions = yap::DataPartitionBlock::create(data, 1);

    yap::ModelIntegral integral(*M);

    yap::ImportanceSampler::calculate(integral, data, true);

    for (auto part : particles(*M)) {

        auto decayingPart = std::dynamic_pointer_cast<yap::DecayingParticle>(part);

        if (not decayingPart)
            continue;

        auto shape = decayingPart->massShape();

        if (not shape)
            continue;

        LOG(INFO) << decayingPart->name();

        std::map<double, std::complex<double> > massShape; // complex value vs mass
        const auto pc = shape->symmetrizationIndices().begin()->first;
        for (auto d : data) {
            double mass = M->fourMomenta()->m(d, pc);
            auto value = shape->value(d, pc);
            massShape[mass] = value;
        }

        // graphs
        TGraph g_re;
        TGraph g_im;
        TGraph g_norm;
        TGraph g_arg;
        TGraph g_argand;


        double delta_m = 0.005;
        double m_prev = 0;

        for (auto mv : massShape) {
            double m = mv.first;
            if (fabs(m - m_prev) < delta_m)
                continue;

            auto value = mv.second;

            double argument = yap::deg(arg(value));
            if (argument < 0.)
                argument += 180.;

            g_re.SetPoint(g_re.GetN(), m, real(value));
            g_im.SetPoint(g_im.GetN(), m, imag(value));
            g_norm.SetPoint(g_norm.GetN(), m, norm(value));
            g_arg.SetPoint(g_arg.GetN(), m, argument);
            g_argand.SetPoint(g_argand.GetN(), real(value), imag(value));
        }

        TString filename =  "mass_shape_" + decayingPart->name() + ".root";

        filename.ReplaceAll("(", "_");
        filename.ReplaceAll(")", "");
        filename.ReplaceAll("+", "p");
        filename.ReplaceAll("-", "m");

        TFile file(filename, "RECREATE");
        file.cd();

        for (auto g : {g_re, g_im, g_norm, g_arg, g_argand})
            g.Draw("AC");

        g_re.GetXaxis()->SetTitle("m [GeV/c^{2}]");
        g_im.GetXaxis()->SetTitle("m [GeV/c^{2}]");
        g_norm.GetXaxis()->SetTitle("m [GeV/c^{2}]");
        g_arg.GetXaxis()->SetTitle("m [GeV/c^{2}]");
        g_argand.GetXaxis()->SetTitle("Re(#Delta)");

        g_re.GetYaxis()->SetTitle("Re(#Delta)");
        g_im.GetYaxis()->SetTitle("Im(#Delta)");
        g_norm.GetYaxis()->SetTitle("#left|#Delta#right|^{2}");
        g_arg.GetYaxis()->SetTitle("Arg(#Delta) [deg]");
        g_argand.GetYaxis()->SetTitle("Im(#Delta)");

        file.WriteTObject(&g_re, "real_vs_m");
        file.WriteTObject(&g_im, "imag_vs_m");
        file.WriteTObject(&g_norm, "norm_vs_m");
        file.WriteTObject(&g_arg, "arg_vs_m");
        file.WriteTObject(&g_argand, "argand");

        file.Close();
    }

    return 0;
}



