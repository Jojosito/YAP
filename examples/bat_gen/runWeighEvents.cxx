// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include <BAT/BCLog.h>
#include <BAT/BCAux.h>

#include <FourMomenta.h>
#include <HelicityFormalism.h>
#include <logging.h>
#include <make_unique.h>
#include <MassAxes.h>

#include "bat_gen.h"
#include "models/d3pi.h"
#include "models/d3pi_phsp.h"
#include "models/d4pi.h"
#include "models/dkkpi.h"
#include "models/D_K0pi0pi0.h"
#include "tools.h"

#include <TFile.h>
#include <TTree.h>

#include <chrono>

using namespace std;
using namespace yap;

int main()
{
    plainLogs(el::Level::Info);

    const double D0_mass = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/d4pi.pdl")["D0"].mass();

    auto m = d4pi();

    TFile* f = TFile::Open("/home/ne53mad/YAPfork/buildRelease/output/D4pi_phsp__mcmc.root", "READ");
    TTree* tree = static_cast<TTree*>(f->Get("D4pi_phsp__mcmc"));

    double m2_01, m2_12, m2_23, m2_02, m2_13;
    tree->SetBranchAddress("m2_01", &m2_01);
    tree->SetBranchAddress("m2_12", &m2_12);
    tree->SetBranchAddress("m2_23", &m2_23);
    tree->SetBranchAddress("m2_02", &m2_02);
    tree->SetBranchAddress("m2_13", &m2_13);

    DataSet dataSet(*m);
    dataSet.addEmptyDataPoints(1);
    DataPoint& dp = dataSet.at(0);

    const auto massAxes = m->massAxes();
    //const MassAxes massAxes{{0,1}, {1,2}, {2,3}, {0,2}, {1,3}};
    const auto FSPs = m->finalStateParticles();
    const auto components = m->components();
    assert(components.size() == 1);

    LOG(INFO) << "mass axes: " << to_string(massAxes);

    std::vector<FourVector<double> > P;



    //unsigned N = tree->GetEntries();
    unsigned N = 20000000;

    TFile* f_weight = TFile::Open("/home/ne53mad/YAPfork/buildRelease/output/D4pi_phsp_mcmc_weight_a1_rho_pi_S.root", "RECREATE");

    double weight;
    TTree* weightTree = new TTree("weight", "weight");
    weightTree->Branch("weight", &weight, "weight/D");

    for (unsigned i = 0; i < N; ++i) {

        if (i%10000==0)
            LOG(INFO) << i;

        // lag
        if (i%5 != 0) {
            weight = 0;
            weightTree->Fill();
            continue;
        }

        tree->GetEntry(i);

        P = calculate_four_momenta(D0_mass, FSPs, massAxes,
                                        {m2_01, m2_12, m2_23, m2_02, m2_13});

        // check
        assert ( norm(FourVector<double>(P[0] + P[1])) - m2_01 < 1e-5 );
        assert ( norm(FourVector<double>(P[1] + P[2])) - m2_12 < 1e-5 );
        assert ( norm(FourVector<double>(P[2] + P[3])) - m2_23 < 1e-5 );
        assert ( norm(FourVector<double>(P[0] + P[2])) - m2_02 < 1e-5 );
        assert ( norm(FourVector<double>(P[1] + P[3])) - m2_13 < 1e-5 );
        //

        m->setFinalStateMomenta(dp, P, dataSet);
        m->calculate(dataSet);

        weight = intensity(components[0], dp);
        weightTree->Fill();

        // check
        const auto amp = amplitude(components[0].decayTrees(), dp);

        m->setFinalStateMomenta(dp, {P[2], P[3], P[0], P[1]}, dataSet);
        m->calculate(dataSet);
        assert (weight - intensity(components[0], dp) < 1e-7);
        assert (abs(amp - amplitude(components[0].decayTrees(), dp)) < 1e-7);

        m->setFinalStateMomenta(dp, {P[1], P[0], P[3], P[2]}, dataSet);
        m->calculate(dataSet);
        assert (weight - intensity(components[0], dp) < 1e-7);
        assert (abs(amp - amplitude(components[0].decayTrees(), dp)) < 1e-7);

        m->setFinalStateMomenta(dp, {P[0], P[3], P[2], P[1]}, dataSet);
        m->calculate(dataSet);
        assert (weight - intensity(components[0], dp) < 1e-7);
        assert (abs(amp - amplitude(components[0].decayTrees(), dp)) < 1e-7);
        //


        //LOG(INFO) << "intensity = " << intens;
    }

    weightTree->Write();
    f_weight->Close();


    return 0;
}
