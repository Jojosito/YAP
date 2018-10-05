// ***************************************************************
// This file was created using the bat-project script.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#ifndef __BAT__D4PI_FIT__H
#define __BAT__D4PI_FIT__H

#include "../bat_fit.h"
#include "../bat_fit_fractions.h"
#include "../fit_fitFraction.h"
#include "../tools.h"

#include "d4pi.h"

#include <DTo4piStructs.h>

//#include <BAT/BCAux.h>
#include <BAT/BCGaussianPrior.h>
#include <BAT/BCConstantPrior.h>
//#include <BAT/BCLog.h>
#include <BAT/BCParameterSet.h>

#include <TDirectory.h>
#include <TEventList.h>
#include <TFile.h>
#include <TTree.h>
#include <TKDTreeBinning.h>
#include <TH2D.h>

#include <assert.h>
#include <complex>
#include <future>
#include <memory>
#include <regex>

using namespace yap;

//-------------------------
inline size_t load_data_4pi(yap::DataSet& data, TTree& t, int N, const double BDT_cut,
                     const double K0_cut = 0, // cut away events with pi+pi- mass within +-K0_cut around K0 mass
                     const bool mc_phasespace_only = false) // only load phasespace MC events
{
    const auto old_size = data.size();
    const auto n_entries = t.GetEntries();

    // get event list (apply BDT cut), speeds up loading
    if (BDT_cut >= 0)
        t.Draw(">>eventList", ("BDT >= " + std::to_string(BDT_cut)).c_str(), "");
    else
        t.Draw(">>eventList", ("BDT < " + std::to_string(BDT_cut)).c_str(), "");
    auto eventList = (TEventList*)gDirectory->Get("eventList");
    LOG(INFO) << eventList->GetN() << " events after BDT cut";
    // apply phasespace cut
    if (mc_phasespace_only) {
        t.SetEventList(eventList);
        t.Draw(">>eventListPHSP", "mc_good_D0 && mc_good_DS && mc_direct_4pi", "");
        eventList = (TEventList*)gDirectory->Get("eventListPHSP");
        LOG(INFO) << eventList->GetN() << " events after phasespace cut";
    }
    auto eventListArray = eventList->GetList();

    // set branch addresses
    EVENT* E = nullptr;
    t.SetBranchAddress("E", &E);
    t.SetBranchStatus("*", 0);
    t.SetBranchStatus("mom_uc_pi*", 1);
    t.SetBranchStatus("mom_pi*", 1);

    if (N <= 0) // attempt to load all data
        N = eventList->GetN();
    else
        N = std::min(N, eventList->GetN());

    data.reserve(old_size + N);

    int n_loaded = 0;

    std::vector<yap::FourVector<double>> P; // pi+ pi- pi+ pi-
    std::vector<double> masses;

    // get K0 mass
    const double K0_mass = read_pdl_file((std::string)::getenv("YAPDIR") + "/data/d4pi.pdl")["K0"].mass();
    LOG(INFO) << "K0 mass = " << K0_mass;

    for (long i=0; i<eventList->GetN(); ++i) {
        unsigned long long n = eventListArray[i];

        if (t.GetEntry(n) < 1)
            continue;

        // K0 cut with UNconstrained 4 momenta
        if (K0_cut > 0) {
            masses.clear();
            masses.push_back((E->mom_uc_piPlus1 + E->mom_uc_piMinus1).M());
            masses.push_back((E->mom_uc_piPlus1 + E->mom_uc_piMinus2).M());
            masses.push_back((E->mom_uc_piPlus2 + E->mom_uc_piMinus1).M());
            masses.push_back((E->mom_uc_piPlus2 + E->mom_uc_piMinus2).M());

            bool cut(false);
            for (auto m : masses)
                if (fabs(m - K0_mass) < K0_cut) {
                    cut = true;
                    break;
                }

            if (cut)
                continue;
        }

        // Fill CONSTRAINED 4-momenta
        P.clear();
        // \todo ATTENTION: I had a bug in my analysis,
        // which means piPlus2 and piMinus1 (constrained) are swapped
        P.push_back(convert(E->mom_piPlus1));
        P.push_back(convert(E->mom_piPlus2)); // actually mom_piMinus1
        P.push_back(convert(E->mom_piMinus1));// actually mom_piPlus2
        P.push_back(convert(E->mom_piMinus2));

        data.push_back(P);

        if (++n_loaded >= N)
            break;

        if (n_loaded%10000 == 0)
            std::cout << "\r" << 100. * n_loaded/N << "%                           ";
    }
    std::cout << std::endl;

    if (data.size() == old_size)
        LOG(INFO) << "No data loaded.";
    else {
        LOG(INFO) << "Loaded " << data.size() - old_size << " data points (" << ((data.size() - old_size) * data[0].bytes() * 1.e-6) << " MB)"
                << " from a tree of size " << n_entries;
        if (old_size)
            LOG(INFO) << "Total data size now " << data.size() << " points (" << (data.bytes() * 1.e-6) << " MB)";
    }

    if (int(data.size() - old_size) < N)
        LOG(WARNING) << "could not load as many data points as requested.";

    return data.size() - old_size;
}

//-------------------------
inline void generate_fit_fraction_data(bat_fit& m, unsigned nPoints, unsigned nPartitions)
{
    auto T = yap::read_pdl_file((::getenv("YAPDIR") ? (std::string)::getenv("YAPDIR") + "/data" : ".") + "/evt.pdl");

    //assert(m.model()->initialStates()[0]->name() == "D0");
    auto isp_mass = T["D0"].mass();

    auto A = m.model()->massAxes();
    auto m2r = yap::squared(yap::mass_range(isp_mass, A, m.model()->finalStateParticles()));

    std::mt19937 g(0);
    // fill data set with nPoints points
    std::generate_n(std::back_inserter(m.fitFractionData()), nPoints,
            std::bind(yap::phsp<std::mt19937>, std::cref(*m.model()), isp_mass, A, m2r, g, std::numeric_limits<unsigned>::max()));
    m.fitFractionPartitions() = yap::DataPartitionBlock::create(m.fitFractionData(), nPartitions);

    yap::ImportanceSampler::calculate(m.fitFractionIntegral(), m.fitFractionPartitions(), true);
}

//-------------------------
inline bat_fit d4pi_fit(std::string name, std::vector<double> pars = {})
{
    bat_fit m(name, d4pi(), {});

    // check if pars need to be completed with free parameters
    bool addFreeParams = false;
    if ((int)pars.size() == m.firstParameter())
        addFreeParams = true;

    // free a_1 width
    //auto a_1   = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("a_1+")));
    //m.addParameter("width(a_1)", std::dynamic_pointer_cast<BreitWigner>(a_1->massShape())->width(), 0, 1.4);
    //m.GetParameters().Back().SetPriorConstant();

    if (d4pi_free_parameters) {
        /*if (a_sigma_pi or sigma_f_0_1370) {
            auto sigma  = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(500)")));
            auto width = std::dynamic_pointer_cast<BreitWigner>(sigma->massShape())->width();
            m.addParameter("width(f_0(500))", width, 0.5*width->value(), 1.5*width->value());
            m.GetParameters().Back().SetPriorConstant();
        }*/

        if (d4pi_f_0_1370_pipiS) {
            auto f_0_1370  = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(1370)")));
            auto width = std::dynamic_pointer_cast<BreitWigner>(f_0_1370->massShape())->width();
            m.addParameter("width(f_0(1370))", width, 0.5*width->value(), 1.5*width->value());
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.35, 0.15));
            if (addFreeParams)
                pars.push_back(width->value());
        }

        if (d4pi_a1_bowler) {
            auto a_1   = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("a_1+")));
            auto shape = std::dynamic_pointer_cast<BowlerMassShape>(a_1->massShape());

            double mass = shape->mass()->value();
            m.addParameter("mass(a_1)", shape->mass(), 0.9*mass, 1.1*mass);
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(mass, 0.04));
            if (addFreeParams)
                pars.push_back(mass);

            double width = shape->width()->value();
            m.addParameter("width(a_1)", shape->width(), 0.7*width, 1.5*width);
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(width, 0.175));
            if (addFreeParams)
                pars.push_back(width);

            double coupling = shape->coupling()->value();
            m.addParameter("K*K_coupling", shape->coupling(), 0.5*coupling, 2.*coupling);
            m.GetParameters().Back().SetPriorConstant();
            if (addFreeParams)
                pars.push_back(coupling);
        }

        /*if (d4pi_pi1300) {
            auto pi1300 = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("pi(1300)+")));
            auto shape  = std::dynamic_pointer_cast<ConstantWidthBreitWigner>(pi1300->massShape());

            m.addParameter("mass(pi(1300)+)", shape->mass(), 1.1, 1.5);
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(1.3, 0.1));

            m.addParameter("width(pi(1300)+)", shape->width(), 0.2, 0.6);
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(0.4, 0.2));
        }*/

        // f_0(980)
        if (not particles(*m.model(), yap::is_named("f_0")).empty()) {
            auto f_0 = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0")));
            auto f_0_980_flatte = std::dynamic_pointer_cast<Flatte>(f_0->massShape());

            double mass = f_0_980_flatte->mass()->value();
            m.addParameter("mass(f_0)", f_0_980_flatte->mass(), 0.9*mass, 1.1*mass);
            if (addFreeParams)
                pars.push_back(mass);

            /*for (auto ch : f_0_980_flatte->channels()) {
                double coupling = ch.Coupling->value();
                m.addParameter("f_0_coupling_" + ch.Particles[0]->name(), ch.Coupling, 0.5*coupling, 2*coupling);
            }*/
        }

        // rho
        if (not particles(*m.model(), yap::is_named("rho0")).empty()) {
            auto rho = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("rho0")));
            auto shape = std::dynamic_pointer_cast<GounarisSakurai>(rho->massShape());

            double mass = shape->mass()->value();
            m.addParameter("mass(rho0)", shape->mass(), 0.9*mass, 1.1*mass);
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(mass, 10.*0.00025));
            if (addFreeParams)
                pars.push_back(mass);

            double width = shape->width()->value();
            m.addParameter("width(rho0)", shape->width(), 0.7*width, 1.3*width);
            m.GetParameters().Back().SetPrior(new BCGaussianPrior(width, 5.*0.008));
            if (addFreeParams)
                pars.push_back(width);
        }
    }

    if (d4pi_bg_K) {
        auto K0   = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("K0")));
        auto shape = std::dynamic_pointer_cast<BreitWigner>(K0->massShape());

        m.addParameter("width(K0)", shape->width(), 0, 0.2);
        m.GetParameters().Back().SetPriorConstant();
        if (addFreeParams)
            pars.push_back(shape->width()->value());
    }

    if (d4pi_free_parameters and d4pi_bg_sigma /*and d4pi_bg_only*/) {
        auto sigma   = std::static_pointer_cast<DecayingParticle>(particle(*m.model(), is_named("f_0(500)")));
        auto shape = std::dynamic_pointer_cast<PoleMass>(sigma->massShape());

        // mass is complex
        m.addParameter("mass(f_0(500))", shape->mass(), std::complex<double>(0.2, -0.7), std::complex<double>(0.3, -0.16));
        if (addFreeParams) {
            pars.push_back(0);
            pars.push_back(0);
        }
    }

    if (not pars.empty()) {
        m.setParameters(pars);
        LOG(INFO) << "Amplitudes after setting parameters: ";
        printAmplitudes(*m.model());
    }

    if (d4pi_fix_amplitudes) {
        //fix_free_amplitudes(*m.model());
        m.fixAmplitudes();
        m.fixAdmixtures();
    }

    // set priors: range around expected value
    double range = 2.;
    for (const auto& fa : m.freeAmplitudes()) {
        //double re = real(fa->value());
        //double im = imag(fa->value());

        double ab = abs(fa->value());
        double ar = deg(arg(fa->value()));
        m.setPriors(fa, new BCConstantPrior(0, range*ab),
                new BCConstantPrior(ar - 180, ar + 180));
        m.setRealImagRanges(fa, -range*ab, range*ab, -range*ab, range*ab);
        m.setAbsArgRanges(fa, 0, range*ab,
                          ar - 180, ar + 180);
    }
    //m.setRealImagRanges(range);

    LOG(INFO) << "Amplitudes: ";
    printAmplitudes(*m.model());

    return m;
}

//-------------------------
inline bat_fit_fractions d4pi_fit_fractions(std::string name, std::vector<double> pars = {})
{
    bat_fit_fractions m(name, d4pi(), {});

    // check if pars need to be completed with free parameters
    bool addFreeParams = false;
    if ((int)pars.size() == m.firstParameter())
        addFreeParams = true;



    if (not pars.empty()) {
        m.setParameters(pars);
        LOG(INFO) << "Amplitudes after setting parameters: ";
        printAmplitudes(*m.model());
    }

    if (d4pi_fix_amplitudes) {
        //fix_free_amplitudes(*m.model());
        m.fixAmplitudes();
        m.fixAdmixtures();
    }

    // set priors: range around expected value
    double range = 2.;
    for (const auto& fa : m.freeAmplitudes()) {
        //double re = real(fa->value());
        //double im = imag(fa->value());

        double ab = abs(fa->value());
        double ar = deg(arg(fa->value()));
        m.setPriors(fa, new BCConstantPrior(0, range*ab),
                new BCConstantPrior(ar - 180, ar + 180));
        m.setRealImagRanges(fa, -range*ab, range*ab, -range*ab, range*ab);
        m.setAbsArgRanges(fa, 0, range*ab,
                          ar - 180, ar + 180);
    }
    //m.setRealImagRanges(range);

    LOG(INFO) << "Amplitudes: ";
    printAmplitudes(*m.model());

    return m;
}


//-------------------------
inline std::vector<double> get_bat_parameter_values(bat_fit_fractions& m)
{
    std::vector<double> parameterValues;
    auto& batParameters = m.GetParameters();
    for (unsigned i = 0; i < batParameters.Size(); ++i)
    {
        auto& batPar = batParameters.At(i);
        double val(0.);

        if (batPar.Fixed())
        {
            val = batPar.GetFixedValue();
        }
        else
        {
            val = batPar.GetPriorMean();
        }

        parameterValues.push_back(val);
    }

    return parameterValues;
}


//-------------------------
inline bat_fit_fractions d4pi_fit_fractions(std::string name, std::string pars)
{
    bat_fit_fractions m(name, d4pi(), {});

    LOG(INFO) << "BAT parameters:";
    auto& batParameters = m.GetParameters();
    for (unsigned i = 0; i < batParameters.Size(); ++i)
    {
        auto& batPar = batParameters.At(i);
        LOG(INFO) << batPar.GetName();

        // default value: fixed to 0
        batPar.Fix(0);
        batPar.SetLimits(0, 0);
    }


    //
    // set BAT parameters
    //
    LOG(INFO) << "trying to match parameters:";
    std::istringstream parsSs(pars);
    for (std::string line; std::getline(parsSs, line); )
    {
        // remove parameter id from line
        std::regex regexReal("real\\(\\d+");
        std::regex regexImag("imag\\(\\d+");
        std::regex regexOther("\"\\(\\d+");

        line = std::regex_replace(line, regexReal, "real(");
        line = std::regex_replace(line, regexImag, "imag(");
        line = std::regex_replace(line, regexOther, "\"");

        bool found = false;

        // loop over parameters
        for (unsigned i = 0; i < batParameters.Size(); ++i)
        {
            auto& batPar = batParameters.At(i);

            // remove parameter id from batParName
            std::string batParName = batPar.GetName();
            batParName = std::regex_replace(batParName, regexReal, "real(");
            batParName = std::regex_replace(batParName, regexImag, "imag(");
            batParName = std::regex_replace(batParName, std::regex("^\\d+"), "");

            if (line.find(batParName) != std::string::npos) {
                LOG(INFO) << "found param " << batParName;

                // parse string
                std::regex regexValueAndError("[0-9\\ +-\\.]+$");
                std::smatch valueAndErrorMatch;

                std::regex_search(line, valueAndErrorMatch, regexValueAndError);
                std::string valueAndError = valueAndErrorMatch.str(0);

                LOG(INFO) << "line = " << line;
                //LOG(INFO) << "valueAndErrorMatch = " << valueAndErrorMatch;
                LOG(INFO) << "valueAndError = " << valueAndError;

                std::string plusMinus(" +- ");
                std::size_t plusMinusFound = valueAndError.find(plusMinus);

                double value = std::stod(valueAndError);
                double error = std::stod(valueAndError.substr(plusMinusFound + plusMinus.size()));

                LOG(INFO) << "parsed value = " << value << " +- " << error;

                // set priors
                batPar.Unfix();
                batPar.SetPrior(new BCGaussianPrior(value, error));
                batPar.SetLimits(value - 5.*error, value + 5.*error);

                found = true;
                break;
            }
        }
        if (!found) {
            LOG(INFO) << "NOT found param for " << line;
            exit(0);
        }
    }

    // set YAP parameters
    m.setParameters(get_bat_parameter_values(m));


    // print
    LOG(INFO) << "Amplitudes: ";
    printAmplitudes(*m.model());
    m.PrintSummary();

    return m;
}


//-------------------------
std::string to_short_string(const DecayTree& dt)
{
    std::string s = yap::to_string(dt);

    s = std::regex_replace(s, std::regex("D0 \\[m = 0\\] --> "), "");
    s = std::regex_replace(s, std::regex(".*pi\\+.*pi\\-.*"), "");
    s = std::regex_replace(s, std::regex("\\[m = ..\\]"), "");
    s = std::regex_replace(s, std::regex("\\[m = .\\]"), "");
    s = std::regex_replace(s, std::regex(".*-->"), "_");
    s = std::regex_replace(s, std::regex(", S = .."), "");
    s = std::regex_replace(s, std::regex(", S = ."), "");
    s = std::regex_replace(s, std::regex("\n"), ",");
    s = std::regex_replace(s, std::regex(" "), "");
    s = std::regex_replace(s, std::regex(" "), "");
    s = std::regex_replace(s, std::regex("pi[\\+\\-],"), "");
    s = std::regex_replace(s, std::regex(",,"), ",");
    s = std::regex_replace(s, std::regex(",$"), "");
    s = std::regex_replace(s, std::regex(",$"), "");
    s = std::regex_replace(s, std::regex("\\("), "_");
    s = std::regex_replace(s, std::regex("\\)"), "");
    s = std::regex_replace(s, std::regex("\\+"), "_plus");
    s = std::regex_replace(s, std::regex("\\-"), "_minus");
    s = std::regex_replace(s, std::regex(",L=0"), ",S");
    s = std::regex_replace(s, std::regex(",L=1"), ",P");
    s = std::regex_replace(s, std::regex(",L=2"), ",D");
    s = std::regex_replace(s, std::regex(",L=3"), ",F");

    return s;
}

//-------------------------
void draw_overlap_integrals(const DecayTreeVectorIntegral& dtvi)
{
    auto& diagonals = dtvi.diagonals();
    auto offDiagonals = dtvi.offDiagonals();

    TFile file("output/overlapIntegrals.root", "RECREATE");

    // for all decay trees
    {
        TH2D histo("overlapIntegrals", "overlapIntegrals", diagonals.size(), 0, diagonals.size(), diagonals.size(), 0, diagonals.size());

        for (unsigned i = 0; i < diagonals.size(); ++i) {
            auto str = to_string(*dtvi.decayTrees().at(i)->freeAmplitude()->decayChannel());
            str.erase(0, 7);

            histo.GetXaxis()->SetBinLabel(i+1, str.c_str());
            histo.GetYaxis()->SetBinLabel(diagonals.size()-i, str.c_str());
        }

        for (unsigned i = 0; i < diagonals.size(); ++i) {
            for (unsigned j = i; j < diagonals.size(); ++j) {
                double val = 0.;
                if (i == j) {
                    val = 1.;
                }
                else {
                    double norm = 1. / sqrt(diagonals.at(i).value() * diagonals.at(j).value());
                    val = norm * abs(offDiagonals.at(i).at(j-i-1).value());
                }
                histo.Fill(j+0.5, diagonals.size()-i-0.5, val);
            }
        }

        histo.Write();
    }

    // for grouped decay trees
    {
        auto groupedDecayTrees = group_by_free_amplitudes(dtvi.decayTrees());

        // new grouped matrix
        RealIntegralElementVector groupedDiagonals(groupedDecayTrees.size());
        ComplexIntegralElementMatrix groupedOffDiagonals(groupedDecayTrees.size() - 1);
        for (size_t i = 0; i < groupedOffDiagonals.size(); ++i)
            groupedOffDiagonals[i] = ComplexIntegralElementMatrix::value_type(groupedDecayTrees.size() - i - 1);

        std::vector<unsigned> limits({0});
        for (auto dts : groupedDecayTrees)
            limits.push_back(limits.back() + dts.size());

        // fill grouped matrix
        for (unsigned iGroup = 0; iGroup < groupedDecayTrees.size(); ++iGroup) {
            for (unsigned jGroup = iGroup; jGroup < groupedDecayTrees.size(); ++jGroup) {
                ComplexIntegralElement element(0);
                for (unsigned i = limits.at(iGroup); i < limits.at(iGroup + 1); ++ i)
                    for (unsigned j = limits.at(jGroup); j < limits.at(jGroup + 1); ++ j)
                        element += dtvi.component(i, j);

                if (iGroup == jGroup) {
                    assert(fabs(imag(element.value())) < 1.e-7);
                    groupedDiagonals.at(iGroup).value() = real(element.value());
                }
                else
                    groupedOffDiagonals.at(iGroup).at(jGroup - iGroup - 1) = element;
            }
        }


        TH2D histo("groupedOverlapIntegrals", "groupedOverlapIntegrals", groupedDiagonals.size(), 0, groupedDiagonals.size(), groupedDiagonals.size(), 0, groupedDiagonals.size());

        for (unsigned i = 0; i < groupedDiagonals.size(); ++i) {
            auto str = to_short_string(*groupedDecayTrees.at(i).at(0));

            histo.GetXaxis()->SetBinLabel(i+1, str.c_str());
            histo.GetYaxis()->SetBinLabel(groupedDiagonals.size()-i, str.c_str());
        }

        for (unsigned i = 0; i < groupedDiagonals.size(); ++i) {
            for (unsigned j = i; j < groupedDiagonals.size(); ++j) {
                double val = 0.;
                if (i == j) {
                    val = 1.;
                }
                else {
                    double norm = 1. / sqrt(groupedDiagonals.at(i).value() * groupedDiagonals.at(j).value());
                    val = norm * abs(groupedOffDiagonals.at(i).at(j-i-1).value());
                }
                histo.Fill(j+0.5, groupedDiagonals.size()-i-0.5, val);
            }
        }

        histo.Write();

    }

    file.Close();
}

//-------------------------
inline void d4pi_normalizeFitFractions(bat_fit& m)
{
    LOG(INFO) << "\nd4pi_normalizeFitFractions";
    auto groupedDecayTrees = group_by_free_amplitudes(m.fitFractionIntegral().integrals().at(0).Integral.decayTrees());

    // aim so that fixed a1 wave will have ~40%
    double sum(0);
    auto ff = fit_fractions(m.fitFractionIntegral().integrals().at(0).Integral, groupedDecayTrees);
    for (size_t i = 0; i < ff.size(); ++i) {

        if (groupedDecayTrees[i].empty())
            continue;
        if (!groupedDecayTrees[i][0])
            continue;

        ff = fit_fractions(m.fitFractionIntegral().integrals().at(0).Integral, groupedDecayTrees);

        sum += ff[i].value();
    }

    const double aim = 0.6 * sum / (groupedDecayTrees.size() - 1);

    for (size_t i = 0; i < ff.size(); ++i) {

        if (groupedDecayTrees[i].empty())
            continue;
        if (!groupedDecayTrees[i][0])
            continue;

        ff = fit_fractions(m.fitFractionIntegral().integrals().at(0).Integral, groupedDecayTrees);

        double frac = ff[i].value();

        std::shared_ptr<FreeAmplitude> fa;
        for (auto ddtv : groupedDecayTrees[i][0]->daughterDecayTreeVector()) {
            fa = ddtv->freeAmplitude();
            if (fa->variableStatus() != yap::VariableStatus::fixed)
                break;
        }

        if (not fa or fa->variableStatus() == yap::VariableStatus::fixed)
            fa =  groupedDecayTrees[i][0]->freeAmplitude();

        if (fa and fa->variableStatus() != yap::VariableStatus::fixed)
            fa->setValue(fa->value() / sqrt(frac) * sqrt(aim));


        LOG(INFO) << to_string(*fa) << " = " << fa->value() << "; \t ff =" << frac*100. << " %";

    }
}

//-------------------------
inline void d4pi_printFitFractions(bat_fit& m)
{
    // sum of integrals
    double sumIntegrals(0);
    // loop over admixtures
    for (const auto& mci : m.fitFractionIntegral().integrals()) {
        sumIntegrals += mci.Admixture->value() * integral(mci.Integral).value();
    }
    LOG(INFO) << "\nsum of integrals = " << sumIntegrals;

    double sum(0);

    LOG(INFO) << "\nFit fractions for single decay trees:";
    for (const auto& mci : m.fitFractionIntegral().integrals()) {
        double sum_admixture(0);
        auto ff = fit_fractions(mci.Integral);
        for (size_t i = 0; i < ff.size(); ++i) {
            double fit_frac = mci.Admixture->value()  * integral(mci.Integral).value() / sumIntegrals * ff[i].value();
            LOG(INFO) << to_string(*mci.Integral.decayTrees()[i]) << "\t" << fit_frac*100. << " %; fit fraction in admixture = " << ff[i].value()*100. << " %";
            sum += fit_frac;
            sum_admixture += ff[i].value();
        }
        LOG(INFO) << "Sum for admixture = " << sum_admixture*100 << " %";
        LOG(INFO) << "Deviation from 100% = " << (sum_admixture-1.)*100 << " %";
    }
    LOG(INFO) << "Sum = " << sum*100 << " %";
    LOG(INFO) << "Deviation from 100% = " << (sum-1.)*100 << " %";

    draw_overlap_integrals(m.fitFractionIntegral().integrals()[0].Integral);



    LOG(INFO) << "\nFit fractions grouped by free amplitudes:";
    sum = 0;
    std::multimap<double, unsigned> ff_map;
    auto groupedDecayTrees = group_by_free_amplitudes(m.fitFractionIntegral().integrals().at(0).Integral.decayTrees());
    auto ff = fit_fractions(m.fitFractionIntegral().integrals().at(0).Integral, groupedDecayTrees);
    for (size_t i = 0; i < ff.size(); ++i) {

        if (groupedDecayTrees[i].empty())
            continue;
        if (!groupedDecayTrees[i][0])
            continue;

        LOG(INFO) << "" ;
        LOG(INFO) << yap::daughter_name()(*(groupedDecayTrees[i][0]))
                  << "   L = " << yap::orbital_angular_momentum()(*(groupedDecayTrees[i][0]))
                  << " :: " << std::to_string(groupedDecayTrees[i].size()) << " amplitudes:";
        for (const auto& dt : groupedDecayTrees[i])
            LOG(INFO) << yap::to_string(*dt);

        double frac = ff[i].value();

        LOG(INFO) << "\t" << frac*100. << " %";
        sum += frac;

        ff_map.insert( std::pair<double, unsigned>(-frac, i) );
    }
    LOG(INFO) << "Sum = " << sum*100 << " %\n";


    LOG(INFO) << "\nFit fractions grouped by free amplitudes, sorted by decreasing fit fraction:";
    for (auto kv : ff_map) {
        unsigned i = kv.second;

        LOG(INFO) << "" ;
        LOG(INFO) << yap::daughter_name()(*(groupedDecayTrees[i][0]))
                  << "   L = " << yap::orbital_angular_momentum()(*(groupedDecayTrees[i][0]))
                  << " :: " << std::to_string(groupedDecayTrees[i].size()) << " amplitudes:";
        for (const auto& dt : groupedDecayTrees[i])
            LOG(INFO) << yap::to_string(*dt);
        LOG(INFO) << "\t" << -kv.first*100. << " %";
    }
}

//-------------------------
inline double chi2_adaptiveBinning(bat_fit& batFit, bool sigma = false, unsigned EvtsPerBin = 30)
{
    auto massAxes = batFit.model()->massAxes({{0, 1}, {0, 3}, {1, 2}, {2, 3}, {0, 2}});

    unsigned dataDim = 5;
    unsigned dataSize_meas = batFit.fitData().size();
    std::vector<double> data_meas; // invariant mass squares; stride = 1 for the points and = N (the dataSize) for the dimension
    data_meas.reserve(dataDim * dataSize_meas);

    // fill data_meas
    for (auto& indices : massAxes) { // must first loop over mass axes
        for (auto& dp : batFit.fitData()) { // then over data points
            FourVector<double> s;
            for (auto fsp : indices->indices())
                s += batFit.model()->fourMomenta()->finalStateMomenta(dp)[fsp];

            data_meas.push_back(s*s);
        }
    }

    // adaptive binning
    TKDTreeBinning kdBins(dataSize_meas, dataDim, &data_meas[0], dataSize_meas/EvtsPerBin);

    /*std::cout << "bin contents (measured):\n";
    for (unsigned i = 0; i < kdBins.GetNBins(); ++i) {
        std::cout << kdBins.GetBinContent(i) << "  ";
    }
    std::cout << "\n";*/


    // bin integralData
    std::vector<double> histo_exp(kdBins.GetNBins(), 0.);
    std::vector<double> histo_exp_sigma2(kdBins.GetNBins(), 0.);
    std::vector<unsigned> histo_exp_N(kdBins.GetNBins(), 0.);
    std::vector<double> ss(5, 0.); // invariant mass squares

    for (auto& dp : batFit.integralData()) {
        for (unsigned i = 0; i < massAxes.size(); ++i) {
            FourVector<double> s;
            for (auto fsp : massAxes[i]->indices())
                s += batFit.model()->fourMomenta()->finalStateMomenta(dp)[fsp];

            ss[i] = s*s;
        }

        // fill bin with intensity
        histo_exp.at(kdBins.FindBin(&ss[0])) += intensity(*batFit.model(), dp);
        histo_exp_sigma2.at(kdBins.FindBin(&ss[0])) += pow(intensity(*batFit.model(), dp), 2); // compare arXiv:1712.08609
        histo_exp_N.at(kdBins.FindBin(&ss[0])) += 1;
    }

    /*std::cout << "bin contents (expected):\n";
    for (unsigned i = 0; i < histo_exp_N.size(); ++i) {
        std::cout << histo_exp_N.at(i) << "  ";
    }
    std::cout << "\n";

    std::cout << "bin contents (expected) * intensity:\n";
    for (unsigned i = 0; i < histo_exp.size(); ++i) {
        std::cout << histo_exp.at(i) << "  ";
    }
    std::cout << "\n";*/


    // calculate chi^2
    double n_exp = std::accumulate(histo_exp.begin(), histo_exp.end(), 0.0);

    double norm_exp = dataSize_meas/n_exp;

    double chi2(0.);
    double NDF = -double(batFit.GetNFreeParameters());

    for (unsigned i = 0; i < histo_exp.size(); ++i) {
        // at least this much events
        if (histo_exp_N.at(i) < 5.) {
            //std::cout<< "x  ";
            continue;
        }

        double NbExp = histo_exp.at(i) * norm_exp;
        double Nb = kdBins.GetBinContent(i);
        double sigma2 = 0.;
        if (sigma)
            sigma2 = histo_exp_sigma2.at(i) * norm_exp * norm_exp;

        double chi2inc = pow(Nb - NbExp, 2) / (NbExp + sigma2);
        //std::cout<< chi2inc << "  ";
        chi2 += chi2inc;
        NDF += 1.;
    }
    std::cout << "\n";

    if (sigma)
        std::cout << "total chi2 with adaptive binning, taking uncertainty into account = " << chi2 << "\n";
    else
        std::cout << "total chi2 with adaptive binning = " << chi2 << "\n";
    std::cout << "NDF = " << NDF << "\n";
    std::cout << "reduced chi2 = " << chi2/NDF << "\n";

    return chi2;

}

//-------------------------
inline double chi2(bat_fit& batFit, std::vector<double> sbinUpperBounds = {0.237108, 0.549527, 0.734061})
{
    unsigned nBins = sbinUpperBounds.size() + 1;

    //batFit.model()->calculate(batFit.integralData());
    //batFit.model()->calculate(batFit.fitData());

    auto massAxes = batFit.model()->massAxes({{0, 1}, {0, 3}, {1, 2}, {2, 3}, {0, 2}});
    //auto massAxes = batFit.model()->massAxes();

    std::vector<double> ss(5, 0.); // invariant mass squares
    std::vector<unsigned> pointCoords(5, 0); // in 5D space

    std::vector<unsigned> histo_exp_N(pow(nBins, 5), 0.); // predicted, number of events in bins
    std::vector<double> histo_exp(pow(nBins, 5), 0.); // predicted
    std::vector<double> histo_meas(pow(nBins, 5), 0.); // measured

    // phsp data -> predicted events
    for (auto& dp : batFit.integralData()) {
        for (unsigned i = 0; i < massAxes.size(); ++i) {
            FourVector<double> s;
            for (auto fsp : massAxes[i]->indices())
                s += batFit.model()->fourMomenta()->finalStateMomenta(dp)[fsp];

            ss[i] = s*s;
        }

        // determine bin number
        unsigned iCoord(0);
        for (auto& c : ss) {
          for (unsigned i = 0; i < sbinUpperBounds.size(); ++i) {
            if (c < sbinUpperBounds.at(i)) {
              pointCoords.at(iCoord) = i;
              break;
            }
            pointCoords.at(iCoord) = sbinUpperBounds.size();
          }
          ++iCoord;
        }

        unsigned iBin(0);
        for (unsigned i=0; i<pointCoords.size(); ++i) {
          iBin += pointCoords.at(i) * pow(nBins, i);
        }

        // fill bin with intensity
        histo_exp.at(iBin) += intensity(*batFit.model(), dp);
        histo_exp_N.at(iBin) += 1;
    }


    /*std::cout << "histo_exp\n";
    for (auto c : histo_exp) {
      std::cout<< c << "  ";
    }
    std::cout << "\n";*/


    // data -> measured events
    for (auto& dp : batFit.fitData()) {
        for (unsigned i = 0; i < massAxes.size(); ++i) {
            FourVector<double> s;
            for (auto fsp : massAxes[i]->indices())
                s += batFit.model()->fourMomenta()->finalStateMomenta(dp)[fsp];

            ss[i] = s*s;
        }

        // determine bin number
        unsigned iCoord(0);
        for (auto& c : ss) {
          for (unsigned i = 0; i < sbinUpperBounds.size(); ++i) {
            if (c < sbinUpperBounds.at(i)) {
              pointCoords.at(iCoord) = i;
              break;
            }
            pointCoords.at(iCoord) = sbinUpperBounds.size();
          }
          ++iCoord;
        }

        unsigned iBin(0);
        for (unsigned i=0; i<pointCoords.size(); ++i) {
          iBin += pointCoords.at(i) * pow(nBins, i);
        }

        // fill bin with 1
        histo_meas.at(iBin) += 1;
    }


    /*std::cout << "histo_meas\n";
    for (auto c : histo_meas) {
      std::cout<< c << "  ";
    }
    std::cout << "\n";


    std::cout << "chi2 increments\n";*/

    // calculate chi^2
    double n_exp = std::accumulate(histo_exp.begin(), histo_exp.end(), 0.0);
    double n_meas = std::accumulate(histo_meas.begin(), histo_meas.end(), 0.0);

    double norm_exp = n_meas/n_exp;

    double chi2(0.);
    double NDF = -double(batFit.GetNFreeParameters());

    // see arXiv:1611.09253v2 [hep-ex], Eq. 6
    for (unsigned i = 0; i < histo_exp.size(); ++i) {
        double NbExp = histo_exp.at(i) * norm_exp;

        // at least this much events
        if (histo_exp_N.at(i) < 5.) {
            //std::cout<< "x  ";
            continue;
        }

        double Nb = histo_meas.at(i);

        double chi2inc = pow(Nb - NbExp, 2) / NbExp;
        //std::cout<< chi2inc << "  ";
        chi2 += chi2inc;
        NDF += 1.;
    }
    std::cout << "\n";

    std::cout << "total chi2 = " << chi2 << "\n";
    std::cout << "NDF = " << NDF << "\n";
    std::cout << "reduced chi2 = " << chi2/NDF << "\n";

    return chi2;
}

//-------------------------
std::string forRootAsTxt(bat_fit& m, std::vector<double> params) {
    std::string descriptor;
    std::string data;

    // L
    double L = m.LogLikelihood(params);

    descriptor += "LogLikelihood/F:";
    data += std::to_string(L) + " ";

    // fit fractions
    double sum = 0;
    auto groupedDecayTrees = group_by_free_amplitudes(m.fitFractionIntegral().integrals().at(0).Integral.decayTrees());
    auto ff = fit_fractions(m.fitFractionIntegral().integrals().at(0).Integral, groupedDecayTrees);
    for (size_t i = 0; i < ff.size(); ++i) {

        if (groupedDecayTrees[i].empty())
            continue;
        if (!groupedDecayTrees[i][0])
            continue;

        descriptor += "FF_" + to_short_string(*groupedDecayTrees[i][0]) + "/F:";
        data += std::to_string(ff[i].value()) + " ";

        sum += ff[i].value();
    }

    descriptor += "FF_total/F:";
    data += std::to_string(sum) + " ";

    //
    // background fractions
    //
    // sum of integrals
    double sumIntegrals(0);
    // loop over admixtures
    for (const auto& mci : m.fitFractionIntegral().integrals()) {
        sumIntegrals += mci.Admixture->value() * integral(mci.Integral).value();
    }

    unsigned i_adm(0);
    for (const auto& mci : m.fitFractionIntegral().integrals()) {
        auto ff = fit_fractions(mci.Integral);
        if (ff.size() != 1)
            continue;

        double fit_frac = mci.Admixture->value()  * integral(mci.Integral).value() / sumIntegrals * ff[0].value();

        descriptor += "ADM_" + std::to_string(i_adm++) + "/F:";
        data += std::to_string(fit_frac) + " ";
    }

    //
    // complex amplitudes
    //


    //
    // pars, within limits?
    //
    auto batPars = m.GetParameters();
    assert( params.size() >= batPars.Size());
    bool limit(false);

    for (unsigned i = 0; i < batPars.Size(); ++i) {
        if (batPars[i].IsAtLimit(params[i])) {
            LOG(INFO) << "parameter " << i << " " << batPars[i].GetName() << " at limit:";
            LOG(INFO) << "value = " << params[i] << "; lower limit = " << batPars[i].GetLowerLimit()
                    << "; upper limit = " << batPars[i].GetUpperLimit();
            limit = true;
        }
        descriptor += "PAR_" + batPars[i].GetSafeName() + "/F:";
        data += std::to_string(params[i]) + " ";
    }

    descriptor += "at_limit/B:";
    data += std::to_string(limit) + " ";


    descriptor.pop_back(); // remove last :

    return descriptor + "\n" + data;
}

//-------------------------
void write_weighted_integral_data(bat_fit& m, std::string filename)
{
    TFile* f = new TFile(filename.c_str(), "RECREATE");
    TTree t_data("data", "data");
    TTree t_integralData("integralData", "integral data");
    TTree t_fitFractionData("fitFractionData", "fitFraction data");
    TTree t_fit("fit",  "weighted integration data, aka fit result");

    // prepare storing of helicity angles
    ParticleCombinationSet pcs;
    auto D0s = particles(*m.model(), is_named("D0"));
    if (D0s.size() == 1)
        pcs = (*(D0s.begin()))->particleCombinations();

    // add all daughter non-fsp pcs
    // inefficient, but does the job
    for (unsigned i=0; i<4; ++i)
        for (auto pc : pcs)
            for (auto d : pc->daughters())
                if (not is_final_state_particle_combination(*d))
                    pcs.insert(d);

    std::vector<float> angles(pcs.size() * 2, 0);
    std::vector<std::string> branchnames_angles;

    for (auto pc : pcs) {
        branchnames_angles.push_back("theta_" + to_string(*pc));
        branchnames_angles.push_back("phi_" + to_string(*pc));
    }

    // make name ROOT-branchable
    for (auto& s : branchnames_angles) {
        s = std::regex_replace(s, std::regex("[\\(\\) ]"), "");
        s = std::regex_replace(s, std::regex("->"), "_");
        s = std::regex_replace(s, std::regex("\\+"), "_");
        s = std::regex_replace(s, std::regex(";"), "__");
    }

    float m_12, m_23, m_14, m_34, m_24, m_13,  m_124, m_123, m_234, m_134;

    //
    // data
    //
    {
        t_data.Branch("m_12 ", &m_12,  "m_12/F");
        t_data.Branch("m_23 ", &m_23,  "m_23/F");
        t_data.Branch("m_14 ", &m_14,  "m_14/F");
        t_data.Branch("m_34 ", &m_34,  "m_34/F");
        t_data.Branch("m_24 ", &m_24,  "m_24/F");
        t_data.Branch("m_13 ", &m_13,  "m_13/F");
        t_data.Branch("m_124", &m_124, "m_124/F");
        t_data.Branch("m_123", &m_123, "m_123/F");
        t_data.Branch("m_234", &m_234, "m_234/F");
        t_data.Branch("m_134", &m_134, "m_134/F");

        unsigned i(0);
        for (auto name : branchnames_angles)
            t_data.Branch(name.c_str(), &angles[i++], (name + "/F").c_str());

        for (auto& dp : m.fitData()) {
            auto fsm = m.model()->fourMomenta()->finalStateMomenta(dp);

            m_12  = abs(FourVector<double>(fsm[0] + fsm[1]));
            m_23  = abs(FourVector<double>(fsm[1] + fsm[2]));
            m_14  = abs(FourVector<double>(fsm[0] + fsm[3]));
            m_34  = abs(FourVector<double>(fsm[2] + fsm[3]));
            m_24  = abs(FourVector<double>(fsm[1] + fsm[3]));
            m_13  = abs(FourVector<double>(fsm[0] + fsm[2]));
            m_124 = abs(FourVector<double>(fsm[0] + fsm[1] + fsm[3]));
            m_123 = abs(FourVector<double>(fsm[0] + fsm[1] + fsm[2]));
            m_234 = abs(FourVector<double>(fsm[1] + fsm[2] + fsm[3]));
            m_134 = abs(FourVector<double>(fsm[0] + fsm[2] + fsm[3]));

            i = 0;
            for (auto pc : pcs) {
                auto hel_angles = m.model()->helicityAngles()(dp, m.fitData(), pc);
                //LOG(INFO) << "helicity angles for pc " << to_string(*pc) << ": " << angles.phi << ", " << angles.theta;
                angles[i++] = hel_angles.theta;
                angles[i++] = hel_angles.phi;
            }

            t_data.Fill();
        }

        t_data.Write();
    }

    //
    // integralData
    //
    {
        t_integralData.Branch("m_12 ", &m_12,  "m_12/F");
        t_integralData.Branch("m_23 ", &m_23,  "m_23/F");
        t_integralData.Branch("m_14 ", &m_14,  "m_14/F");
        t_integralData.Branch("m_34 ", &m_34,  "m_34/F");
        t_integralData.Branch("m_24 ", &m_24,  "m_24/F");
        t_integralData.Branch("m_13 ", &m_13,  "m_13/F");
        t_integralData.Branch("m_124", &m_124, "m_124/F");
        t_integralData.Branch("m_123", &m_123, "m_123/F");
        t_integralData.Branch("m_234", &m_234, "m_234/F");
        t_integralData.Branch("m_134", &m_134, "m_134/F");

        unsigned i(0);
        for (auto name : branchnames_angles)
            t_integralData.Branch(name.c_str(), &angles[i++], (name + "/F").c_str());

        for (auto& dp : m.integralData()) {
            auto fsm = m.model()->fourMomenta()->finalStateMomenta(dp);

            m_12  = abs(FourVector<double>(fsm[0] + fsm[1]));
            m_23  = abs(FourVector<double>(fsm[1] + fsm[2]));
            m_14  = abs(FourVector<double>(fsm[0] + fsm[3]));
            m_34  = abs(FourVector<double>(fsm[2] + fsm[3]));
            m_24  = abs(FourVector<double>(fsm[1] + fsm[3]));
            m_13  = abs(FourVector<double>(fsm[0] + fsm[2]));
            m_124 = abs(FourVector<double>(fsm[0] + fsm[1] + fsm[3]));
            m_123 = abs(FourVector<double>(fsm[0] + fsm[1] + fsm[2]));
            m_234 = abs(FourVector<double>(fsm[1] + fsm[2] + fsm[3]));
            m_134 = abs(FourVector<double>(fsm[0] + fsm[2] + fsm[3]));

            i = 0;
            for (auto pc : pcs) {
                auto hel_angles = m.model()->helicityAngles()(dp, m.fitData(), pc);
                //LOG(INFO) << "helicity angles for pc " << to_string(*pc) << ": " << angles.phi << ", " << angles.theta;
                angles[i++] = hel_angles.theta;
                angles[i++] = hel_angles.phi;
            }

            t_integralData.Fill();
        }

        t_integralData.Write();
    }


    //
    // fitFractionData
    //
    {
        t_fitFractionData.Branch("m_12 ", &m_12,  "m_12/F");
        t_fitFractionData.Branch("m_23 ", &m_23,  "m_23/F");
        t_fitFractionData.Branch("m_14 ", &m_14,  "m_14/F");
        t_fitFractionData.Branch("m_34 ", &m_34,  "m_34/F");
        t_fitFractionData.Branch("m_24 ", &m_24,  "m_24/F");
        t_fitFractionData.Branch("m_13 ", &m_13,  "m_13/F");
        t_fitFractionData.Branch("m_124", &m_124, "m_124/F");
        t_fitFractionData.Branch("m_123", &m_123, "m_123/F");
        t_fitFractionData.Branch("m_234", &m_234, "m_234/F");
        t_fitFractionData.Branch("m_134", &m_134, "m_134/F");

        unsigned i(0);
        for (auto name : branchnames_angles)
            t_fitFractionData.Branch(name.c_str(), &angles[i++], (name + "/F").c_str());

        for (auto& dp : m.fitFractionData()) {
            auto fsm = m.model()->fourMomenta()->finalStateMomenta(dp);

            m_12  = abs(FourVector<double>(fsm[0] + fsm[1]));
            m_23  = abs(FourVector<double>(fsm[1] + fsm[2]));
            m_14  = abs(FourVector<double>(fsm[0] + fsm[3]));
            m_34  = abs(FourVector<double>(fsm[2] + fsm[3]));
            m_24  = abs(FourVector<double>(fsm[1] + fsm[3]));
            m_13  = abs(FourVector<double>(fsm[0] + fsm[2]));
            m_124 = abs(FourVector<double>(fsm[0] + fsm[1] + fsm[3]));
            m_123 = abs(FourVector<double>(fsm[0] + fsm[1] + fsm[2]));
            m_234 = abs(FourVector<double>(fsm[1] + fsm[2] + fsm[3]));
            m_134 = abs(FourVector<double>(fsm[0] + fsm[2] + fsm[3]));

            i = 0;
            for (auto pc : pcs) {
                auto hel_angles = m.model()->helicityAngles()(dp, m.fitData(), pc);
                //LOG(INFO) << "helicity angles for pc " << to_string(*pc) << ": " << angles.phi << ", " << angles.theta;
                angles[i++] = hel_angles.theta;
                angles[i++] = hel_angles.phi;
            }

            t_fitFractionData.Fill();
        }

        t_fitFractionData.Write();
    }


    //
    // fit
    //
    {
        float intens;

        // prepare writing of individual admixtures' and waves' intensities
        std::vector<std::string> branchnames_admixtures;
        unsigned i_adm(0);
        for (const auto& mc : m.model()->components())
            branchnames_admixtures.push_back("intens_adm_" + std::to_string(i_adm++));
        std::vector<float> intens_admixtures(branchnames_admixtures.size());

        std::vector<std::string> branchnames_waves;
        auto groupedDecayTrees = group_by_free_amplitudes(m.fitFractionIntegral().integrals().at(0).Integral.decayTrees());

        unsigned i = 0;
        for (; i < groupedDecayTrees.size(); ++i)
            branchnames_waves.push_back("intens_" + to_short_string(*groupedDecayTrees[i][0]));

        for (auto dt : m.fitFractionIntegral().integrals().at(0).Integral.decayTrees())
            groupedDecayTrees.push_back({dt});

        for (; i < groupedDecayTrees.size(); ++i)
            branchnames_waves.push_back("intens_" + to_short_string(*groupedDecayTrees[i][0]) + "_" + std::to_string(i));


        std::vector<float> intens_waves(branchnames_waves.size());

        t_fit.Branch("m_12 ", &m_12,  "m_12/F");
        t_fit.Branch("m_23 ", &m_23,  "m_23/F");
        t_fit.Branch("m_14 ", &m_14,  "m_14/F");
        t_fit.Branch("m_34 ", &m_34,  "m_34/F");
        t_fit.Branch("m_24 ", &m_24,  "m_24/F");
        t_fit.Branch("m_13 ", &m_13,  "m_13/F");
        t_fit.Branch("m_124", &m_124, "m_124/F");
        t_fit.Branch("m_123", &m_123, "m_123/F");
        t_fit.Branch("m_234", &m_234, "m_234/F");
        t_fit.Branch("m_134", &m_134, "m_134/F");
        t_fit.Branch("intensity", &intens, "intens/F");

        for (unsigned i = 0; i < branchnames_admixtures.size(); ++i)
            t_fit.Branch(branchnames_admixtures[i].c_str(), &intens_admixtures[i], (branchnames_admixtures[i] + "/F").c_str());

        for (unsigned i = 0; i < branchnames_waves.size(); ++i)
            t_fit.Branch(branchnames_waves[i].c_str(), &intens_waves[i], (branchnames_waves[i] + "/F").c_str());

        for (unsigned i = 0; i < branchnames_angles.size(); ++i)
            t_fit.Branch(branchnames_angles[i].c_str(), &angles[i], (branchnames_angles[i] + "/F").c_str());

        // force calculate all integrals
        yap::ImportanceSampler::calculate(m.modelIntegral(), m.integralPartitions(), true);

        for (auto& dp : m.integralData()) {
            auto fsm = m.model()->fourMomenta()->finalStateMomenta(dp);

            m_12  = abs(FourVector<double>(fsm[0] + fsm[1]));
            m_23  = abs(FourVector<double>(fsm[1] + fsm[2]));
            m_14  = abs(FourVector<double>(fsm[0] + fsm[3]));
            m_34  = abs(FourVector<double>(fsm[2] + fsm[3]));
            m_24  = abs(FourVector<double>(fsm[1] + fsm[3]));
            m_13  = abs(FourVector<double>(fsm[0] + fsm[2]));
            m_124 = abs(FourVector<double>(fsm[0] + fsm[1] + fsm[3]));
            m_123 = abs(FourVector<double>(fsm[0] + fsm[1] + fsm[2]));
            m_234 = abs(FourVector<double>(fsm[1] + fsm[2] + fsm[3]));
            m_134 = abs(FourVector<double>(fsm[0] + fsm[2] + fsm[3]));

            intens = intensity(*m.model(), dp);

            for (unsigned i = 0; i < intens_admixtures.size(); ++i)
                intens_admixtures[i] = intensity(m.model()->components()[i], dp);

            for (unsigned i = 0; i < intens_waves.size(); ++i)
                intens_waves[i] = intensity(groupedDecayTrees[i], dp);

            unsigned i(0);
            for (auto pc : pcs) {
                auto hel_angles = m.model()->helicityAngles()(dp, m.fitData(), pc);
                angles[i++] = hel_angles.theta;
                angles[i++] = hel_angles.phi;
            }

            t_fit.Fill();
        }

        t_fit.Write();
    }

    f->Close();
}


#endif
