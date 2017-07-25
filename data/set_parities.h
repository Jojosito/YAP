/*  YAP - Yet another PWA toolkit
    Copyright 2015, Technische Universitaet Muenchen,
    Authors: Daniel Greenwald, Johannes Rauch

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/// \file

#ifndef yap_set_parities_h
#define yap_set_parities_h

#include <Exceptions.h>
#include <logging.h>
#include <ParticleTable.h>
#include <QuantumNumbers.h>

inline void set_parity(yap::ParticleTable& pdl, const std::string& name, int parity)
{
    try {
        pdl[name].quantumNumbers().setP(parity);
    }
    catch (yap::exceptions::Exception& e) {
        LOG(INFO) << e.what();
    }
}

/// \function set_parities
/// set parities for selected mesons.
/// THIS IS NOT COMPLETE
inline void set_parities(yap::ParticleTable& pdl) {

    // light unflavored mesons
    set_parity(pdl, "pi+", -1);
    set_parity(pdl, "pi-", -1);
    set_parity(pdl, "pi0", -1);

    set_parity(pdl, "pi(1300)+", -1);
    set_parity(pdl, "pi(1300)-", -1);
    set_parity(pdl, "pi(1300)0", -1);

    set_parity(pdl, "pi(1800)+", -1);
    set_parity(pdl, "pi(1800)-", -1);
    set_parity(pdl, "pi(1800)0", -1);

    set_parity(pdl, "eta", -1);

    set_parity(pdl, "f_0(500)", +1);

    set_parity(pdl, "rho+", -1);
    set_parity(pdl, "rho-", -1);
    set_parity(pdl, "rho0", -1);

    set_parity(pdl, "omega", -1);

    set_parity(pdl, "eta'", -1);

    set_parity(pdl, "f_0", +1);

    set_parity(pdl, "a_0+", +1);
    set_parity(pdl, "a_0-", +1);
    set_parity(pdl, "a_00", +1);

    set_parity(pdl, "phi", -1);

    set_parity(pdl, "h_1", +1);

    set_parity(pdl, "b_1+", +1);
    set_parity(pdl, "b_1-", +1);
    set_parity(pdl, "b_10", +1);

    set_parity(pdl, "a_1+", +1);
    set_parity(pdl, "a_1-", +1);
    set_parity(pdl, "a_10", +1);

    set_parity(pdl, "a_1(1420)+", +1);
    set_parity(pdl, "a_1(1420)-", +1);
    set_parity(pdl, "a_1(1420)0", +1);

    set_parity(pdl, "a_1(1640)+", +1);
    set_parity(pdl, "a_1(1640)-", +1);
    set_parity(pdl, "a_1(1640)0", +1);

    set_parity(pdl, "f_2", +1);

    set_parity(pdl, "f_1", +1);

    set_parity(pdl, "eta(2S)", -1);

    set_parity(pdl, "pi(2S)+", -1);
    set_parity(pdl, "pi(2S)-", -1);
    set_parity(pdl, "pi(2S)0", -1);

    set_parity(pdl, "a_2+", +1);
    set_parity(pdl, "a_2-", +1);
    set_parity(pdl, "a_20", +1);

    set_parity(pdl, "eta(1405)", -1);

    set_parity(pdl, "omega(2S)", -1);

    set_parity(pdl, "rho(2S)+", -1);
    set_parity(pdl, "rho(2S)-", -1);
    set_parity(pdl, "rho(2S)0", -1);

    set_parity(pdl, "eta(1475)", -1);

    set_parity(pdl, "f_0(1500)", +1);

    set_parity(pdl, "omega(1650)", -1);

    set_parity(pdl, "omega(3)(1670)", -1);


    // strange mesons
    set_parity(pdl, "K+", -1);
    set_parity(pdl, "K-", -1);

    set_parity(pdl, "K0", -1);
    set_parity(pdl, "anti-K0", -1);

    set_parity(pdl, "K_S0", -1);
    set_parity(pdl, "K_L0", -1);


    // charmed mesons
    set_parity(pdl, "D+", -1);
    set_parity(pdl, "D-", -1);

    set_parity(pdl, "D0", -1);
    set_parity(pdl, "anti-D0", -1);

}

#endif
