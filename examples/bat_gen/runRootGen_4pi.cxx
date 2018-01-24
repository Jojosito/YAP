#include <logging.h>

#include <DTo4piStructs.h>
#include <PDL.h>

#include <TFile.h>
#include <TGenPhaseSpace.h>
#include <TRandom3.h>
#include <TTree.h>

#include <assert.h>


int main()
{
    yap::plainLogs(el::Level::Info);

    auto pdl = yap::read_pdl_file((std::string)::getenv("YAPDIR") + "/data/d4pi.pdl");
    const double D0_mass = pdl["D0"].mass();
    const double pi_mass = pdl["pi+"].mass();

    TLorentzVector D0(0., 0., 0., D0_mass);
    double fsp_masses[4] = {pi_mass, pi_mass, pi_mass, pi_mass};
    TGenPhaseSpace gen;
    assert( gen.SetDecay(D0, 4, fsp_masses) );

    TFile* f = TFile::Open("/home/ne53mad/YAPfork/buildRelease/output/D4pi_phsp_root.root", "RECREATE");

    EVENT* E = new EVENT;
    double BDT = 1;
    TTree* t = new TTree("t", "t");
    t->Branch("E", "EVENT", &E);
    t->Branch("BDT", &BDT, "BDT/D");

    TRandom3 rnd;

    // warm up
    LOG(INFO) << "warm up";
    for (unsigned i = 0; i < 1e7; ++i)
        gen.Generate();

    for (unsigned i = 0; i < 1e7; ++i) {
        if (i%100000 == 0)
            LOG(INFO) << i;

        E->reset();

        // take weight into account; accept generated point with a probability that equals its weight
        while (gen.Generate() < rnd.Uniform()) {
            ;
        }

        E->mom_uc_piPlus1  = *gen.GetDecay(0);
        E->mom_uc_piMinus1 = *gen.GetDecay(1);
        E->mom_uc_piPlus2  = *gen.GetDecay(2);
        E->mom_uc_piMinus2 = *gen.GetDecay(3);

        E->mom_piPlus1  = E->mom_uc_piPlus1 ;
        E->mom_piMinus1 = E->mom_uc_piMinus1;
        E->mom_piPlus2  = E->mom_uc_piPlus2 ;
        E->mom_piMinus2 = E->mom_uc_piMinus2;

        t->Fill();
    }

    t->Write();
    f->Close();

    return 0;
}

