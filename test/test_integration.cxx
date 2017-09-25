#include <catch.hpp>

#include "helperFunctions.h"

#include <CompensatedSum.h>
#include <DataPartition.h>
#include <DataPoint.h>
#include <DecayTree.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <Model.h>
#include <ModelIntegral.h>

#include <future>
#include <memory>
#include <vector>

/**
 *  Test the integration
 */

TEST_CASE("integration")
{
    // disable debug logs in test
    //yap::disableLogs(el::Level::Debug);
    yap::plainLogs(el::Level::Debug);

    const unsigned nPoints = 1000;

    auto M = d4pi();
    auto data = generate_data(*M, nPoints);
    auto partitions = yap::DataPartitionBlock::create(data, 4);

    // integrating lambda
    auto sum_of_intensities = [M](yap::DataPartition& D)
        {
            M->calculate(D);
            yap::CompensatedSum<double> I(0.);
            for (auto& d : D)
                I += intensity(*M, d);
            return I.sum;
        };
    
    // create threads for calculation on each partition
    std::vector<std::future<double> > partial_sums;
    partial_sums.reserve(partitions.size());
    for (auto& P : partitions)
        partial_sums.push_back(std::async(std::launch::async, sum_of_intensities, std::ref(*P)));

    // wait for each partition to finish calculating
    double bruteForceIntegral = std::accumulate(partial_sums.begin(), partial_sums.end(), 0.,
                                                [](double & l, std::future<double>& s)
                                                {return l += s.get();});

    DEBUG("bruteForceIntegral = " << bruteForceIntegral);

    yap::ModelIntegral mi(*M);

    yap::ImportanceSampler::calculate(mi, partitions);
    double smartIntegral = integral(mi).value();
    DEBUG("smartIntegral = " << smartIntegral);

    REQUIRE(bruteForceIntegral / nPoints == Approx(smartIntegral));


    //
    // fit fractions
    //
    double sum(0);
    auto ff = fit_fractions(mi.integrals()[0].Integral);
    for (size_t i = 0; i < ff.size(); ++i) {
        //LOG(INFO) << yap::to_string(*mi.integrals()[0].Integral.decayTrees()[i]) << "\t" << ff[i].value()*100. << "%" ;
        sum += ff[i].value();
    }
    //LOG(INFO) << "Sum = " << sum*100 << " %";

    //
    // check components
    //
    REQUIRE(mi.integrals().size() == 1);
    REQUIRE(M->components().size() == 1);
    auto& dtvi = mi.integrals()[0].Integral;
    auto& mc = M->components()[0];
    std::vector<double> smartIntegrals;
    std::vector<double> bruteForceIntegrals;

    // original free amplitude values
    std::vector<std::complex<double>> freeAmps;

    yap::DecayTreeVector dts;
    for (unsigned i = 0; i < mc.decayTrees().size(); ++i) {
        if (i == 0 or (i > 0 and mc.decayTrees()[i]->freeAmplitude() == mc.decayTrees()[i-1]->freeAmplitude()))
            dts.push_back(mc.decayTrees()[i]);
        else {
            auto I_dt = integral(dtvi, dts);
            smartIntegrals.push_back(I_dt.value());
            //LOG(INFO) << "integral = " << smartIntegrals.back();
            freeAmps.push_back(mc.decayTrees()[i]->freeAmplitude()->value());

            dts.clear();
            dts.push_back(mc.decayTrees()[i]);
        }
    }
    auto I_dt = integral(dtvi, dts);
    smartIntegrals.push_back(I_dt.value());
    //LOG(INFO) << "integral = " << smartIntegrals.back();
    freeAmps.push_back(mc.decayTrees().back()->freeAmplitude()->value());

    // set free amps to 0
    for (auto& dt : mc.decayTrees()) {
        dt->freeAmplitude()->setValue(0);
    }

    unsigned iAmp = 0;
    for (unsigned i = 0; i < mc.decayTrees().size(); ++i) {
        if (i > 0 and mc.decayTrees()[i]->freeAmplitude() == mc.decayTrees()[i-1]->freeAmplitude())
            continue;

        mc.decayTrees()[i]->freeAmplitude()->setValue(freeAmps[iAmp++]);

        /*LOG(INFO) << "free amplitudes";
        for (auto& dt : mc.decayTrees()) {
            LOG(INFO) << "  " << dt->freeAmplitude()->value();
        }*/

        // create threads for calculation on each partition
        std::vector<std::future<double> > partial_sums;
        partial_sums.reserve(partitions.size());
        for (auto& P : partitions)
            partial_sums.push_back(std::async(std::launch::async, sum_of_intensities, std::ref(*P)));

        // wait for each partition to finish calculating
        double bruteForceIntegral = std::accumulate(partial_sums.begin(), partial_sums.end(), 0.,
                                                    [](double & l, std::future<double>& s)
                                                    {return l += s.get();});

        bruteForceIntegrals.push_back(bruteForceIntegral / nPoints);
        //LOG(INFO) << "bruteForceIntegral = " << bruteForceIntegrals.back();
        mc.decayTrees()[i]->freeAmplitude()->setValue(0.);
    }

    REQUIRE(smartIntegrals.size() == bruteForceIntegrals.size());

    for (unsigned i = 0; i < smartIntegrals.size(); ++i) {
        REQUIRE(smartIntegrals[i] == Approx(bruteForceIntegrals[i]));
    }


}

