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

#ifndef yap_ParticleCombinationCache_h
#define yap_ParticleCombinationCache_h

#include "ParticleCombination.h"
#include "ParticleIndex.h"
#include "WeakPtrCache.h"

#include <memory>
#include <set>
#include <string>

namespace yap {

/// \class ParticleCombinationCache
/// \brief Caches list of ParticleCombination's
/// \author Johannes Rauch, Daniel Greenwald

class ParticleCombinationCache : public WeakPtrCache<const ParticleCombination>
{
public:

    /// implements equivalence checking
    bool equiv(const shared_ptr_type& A, const shared_ptr_type& B) const override
    { return ParticleCombination::equivUpAndDown(A, B); }

    /// Default constructor
    ParticleCombinationCache() = default;

    /// Construct cache from vector of ISP's ParticleCombination's
    ParticleCombinationCache(std::vector<shared_ptr_type> V);

    /// retrieve or create final-state particle ParticleCombination
    /// \param index Index of particle
    /// \param two_lambda Spin projection of particle
    shared_ptr_type fsp(ParticleIndex index, int two_lambda = 0)
    { return operator[](create_fsp(index, two_lambda)); }

    /// retrieve or create copy of ParticleCombination with new spin projection
    /// \param other ParticleCombination to copy
    /// \param two_lambda new spin projection of particle
    shared_ptr_type copy(const ParticleCombination& other, int two_lambda)
    { return operator[](create_copy(other, two_lambda)); }
    
    /// retrieve or create composite particle from daughters.
    /// copies daughters into composite, setting copies' parents = shared_from_this()
    /// \param D ParticleCombinationVector of daughters to create composite from
    /// \param two_lambda spin projection of particle
    shared_ptr_type composite(const ParticleCombinationVector& D, int two_lambda = 0)
    { return operator[](create_composite(D, two_lambda)); }

    using WeakPtrCache::find;

    /// retrieve final-state particle ParticleCombination
    /// Does not add to the cache if ParticleCombination is not found.
    /// \param index Index of particle
    /// \param two_lambda Spin projection of particle
    weak_ptr_type find(ParticleIndex index, int two_lambda = 0)
    { return find(create_fsp(index, two_lambda)); }

    /// retrieve copy of ParticleCombination with new spin projection
    /// Does not add to the cache if ParticleCombination is not found.
    /// \param other ParticleCombination to copy
    /// \param two_lambda new spin projection of particle
    weak_ptr_type find(const ParticleCombination& other, int two_lambda)
    { return find(create_copy(other, two_lambda)); }

    /// retrieve composite particle ParticleCombination from cache.
    /// Does not add to the cache if ParticleCombination is not found.
    /// \param D vector of daughters to construct ParticleCombination from
    /// \param two_lambda Spin projection of ParticleCombinatin to create
    /// \return weak_ptr to ParticleCombination; is empty if not found.
    weak_ptr_type find(const ParticleCombinationVector& D, int two_lambda = 0)
    { return find(create_composite(D, two_lambda)); }

    /// check if cache contains element equating to pc
    /// \param I vector of ParticleIndex's to build ParticleCombination from for checking equivalence
    weak_ptr_type find(const std::vector<ParticleIndex>& I) const;

    /// Check consistency of cache.
    bool consistent() const;

protected:

    /// set lineage: copy each daughter, add pc as parent to copy,
    /// swap copy for daughter, and call setLineage on each daughter.
    void setLineage(shared_ptr_type pc);

private:

    /// add to cache
    void addToCache(shared_ptr_type pc) override;

    /// create final-state particle
    shared_ptr_type create_fsp(ParticleIndex index, int two_lambda = 0) const
    { return shared_ptr_type(new ParticleCombination(index, two_lambda)); }

    /// create copy
    shared_ptr_type create_copy(const ParticleCombination& other, int two_lambda) const;

    /// create composite ParticleCombination
    shared_ptr_type create_composite(const ParticleCombinationVector& D, int two_lambda) const;

};

/// convert to string
std::string to_string(const ParticleCombinationCache& C);

/// streamer
inline std::ostream& operator<<(std::ostream& os, const ParticleCombinationCache& C)
{ os << to_string(C); return os; }

}

#endif
