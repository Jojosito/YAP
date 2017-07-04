// ***************************************************************
// This file was created using the bat-project script
// for project yap_test.
// bat-project is part of Bayesian Analysis Toolkit (BAT).
// BAT can be downloaded from http://mpp.mpg.de/bat
// ***************************************************************

#include "bat_fit.h"
#include "models/d3pi.h"
#include "models/d4pi_fit.h"
#include "models/dkkpi.h"
#include "tools.h"

#include <Group.h>
#include <HelicityFormalism.h>
#include <ImportanceSampler.h>
#include <logging.h>
#include <make_unique.h>
#include <MassRange.h>
#include <PHSP.h>
#include <Sort.h>
#include <ZemachFormalism.h>

#include <algorithm>
#include <random>

int main()
{
    yap::plainLogs(el::Level::Info);

    // create model
    bat_fit m(d4pi_fit("D4pi_fit"));
    double D_mass = 1.8648400; // D0

    // get FSP mass ranges
    auto m2r = yap::squared(mass_range(D_mass, m.axes(), m.model()->finalStateParticles()));

    // generate integration data
    generate_fit_fraction_data(m, 100000, 4);

    d4pi_printFitFractions(m);

    return 0;
}
