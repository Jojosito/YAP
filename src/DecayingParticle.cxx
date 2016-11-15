#include "DecayingParticle.h"

#include "BlattWeisskopf.h"
#include "container_utils.h"
#include "DecayChannel.h"
#include "DecayTree.h"
#include "FinalStateParticle.h"
#include "FreeAmplitude.h"
#include "logging.h"
#include "MassShape.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleTable.h"
#include "SpinAmplitude.h"
#include "VariableStatus.h"

#include <functional>
#include <memory>

namespace yap {

//-------------------------
const is_of_type<DecayingParticle> is_decaying_particle{};

//-------------------------
DecayingParticle::DecayingParticle(const std::string& name, const QuantumNumbers& q,
                                   double radial_size, std::shared_ptr<MassShape> mass_shape) :
    Particle(name, q),
    MassShape_(mass_shape),
    RadialSize_(std::make_shared<PositiveRealParameter>(radial_size))
{
    if (MassShape_)
        MassShape_->setOwner(this);
}

//-------------------------
DecayingParticle::DecayingParticle(const ParticleTableEntry& pde, double radial_size,
                                   std::shared_ptr<MassShape> mass_shape) :
    DecayingParticle(pde.name(), pde.quantumNumbers(), radial_size, mass_shape)
{
}

//-------------------------
bool DecayingParticle::consistent() const
{
    bool C = Particle::consistent();

    if (Channels_.empty()) {
        FLOG(ERROR) << "no channels specified.";
        return false; // further checks require at least one channel
    }

    // check no channel is empty
    if (std::any_of(Channels_.begin(), Channels_.end(), std::logical_not<DecayChannelVector::value_type>())) {
        FLOG(ERROR) << "DecayChannel vector contains nullptr";
        C &= false;
    }

    // check consistency of all channels
    std::for_each(Channels_.begin(), Channels_.end(), [&](const std::shared_ptr<DecayChannel>& dc) {if (dc) C &= dc->consistent();});

    if (MassShape_) {
        C &= MassShape_->consistent();

        if (MassShape_->owner() != this) {
            FLOG(ERROR) << "MassShape's owner is not this";
            C &= false;
        }
    }

    return C;
}

//-------------------------
// helper function
bool pcs_has_pc(const ParticleCombinationSet& pcs, const ParticleCombination& pc, size_t d)
{ return pcs.find(pc.daughters()[d]) != pcs.end(); }

//-------------------------
// helper function, checks that decay channel of free amplitude of decay tree has particle combination
bool dt_valid_for_pc(const DecayTreeVector::value_type& dt, const ParticleCombination& pc, size_t d)
{ return pcs_has_pc(dt->freeAmplitude()->decayChannel()->particleCombinations(), pc, d); }

//-------------------------
void DecayingParticle::addDecayChannel(std::shared_ptr<DecayChannel> c)
{
    if (!c)
        throw exceptions::Exception("DecayChannel empty", "DecayingParticle::addDecayChannel");

    if (c->particleCombinations().empty())
        throw exceptions::Exception(std::string("DecayChannel has no ParticleCombinations - ") + to_string(*c),
                                    "DecayingParticle::addDecayChannel");

    // check ISP
    if (!Channels_.empty() and c->model() != model())
        throw exceptions::Exception("Model mismatch", "DecayingParticle::addDecayChannel");

    // check charge
    if (charge(c->daughters()) != quantumNumbers().Q())
        throw exceptions::Exception("Charge of channel not equal to decaying particle ("
                                    + std::to_string(charge(c->daughters())) + " != " + std::to_string(quantumNumbers().Q()) + ")",
                                    "DecayChannel::addDecayChannel");

    // check that DecayChannel has SpinAmplitudes
    if (c->spinAmplitudes().empty())
        throw exceptions::Exception("DecayChannel has no SpinAmplitudes", "DecayingState::addDecayChannel");

    // check with mass shape
    if (MassShape_)
        MassShape_->checkDecayChannel(*c);
    
    Channels_.push_back(c);

    // create necessary BlattWeisskopf objects
    for (const auto& sa : Channels_.back()->spinAmplitudes()) {
        // if BW is not already stored for L, add it
        if (BlattWeisskopfs_.find(sa->L()) == BlattWeisskopfs_.end())
            BlattWeisskopfs_.emplace(sa->L(), std::make_shared<BlattWeisskopf>(sa->L(), this));
    }

    // add particle combinations
    for (auto pc : Channels_.back()->particleCombinations())
        addParticleCombination(*pc);

    /////////////////////////
    /// create decay trees for channel:

    /// loop over spin amplitudes of channel
    for (auto& sa : Channels_.back()->spinAmplitudes()) {

        // create new FreeAmplitude
        auto fa = std::make_shared<FreeAmplitude>(Channels_.back(), sa);

        // loop over possible parent spin projections
        for (const auto& two_M : sa->twoM()) {

            // loop over possible daughter spin projections
            for (const auto& two_m : sa->twoM(two_M)) {

                // loop over particle combinations of this decaying particle
                for (const auto& pc : particleCombinations()) {

                    // initialize vector of possible decay trees from initial DecayTree created above
                    DecayTreeVector DTV(1, std::make_shared<DecayTree>(fa, two_M, two_m));
                    modifyDecayTree(*DTV[0]);

                    // loop over daughters in channel
                    for (size_t d = 0; d < Channels_.back()->daughters().size() and !DTV.empty(); ++d) {

                        if (is_decaying_particle(Channels_.back()->daughters()[d])) {

                            // get daughter's DecayTreeVectorMap
                            const auto& m_dtv_map = static_cast<const DecayingParticle&>(*Channels_.back()->daughters()[d]).decayTrees();

                            // get DecayTreeVector for desired spin projection of daughter
                            auto it  = m_dtv_map.find(two_m[d]);
                            if (it == m_dtv_map.end()) {
                                DTV.clear();
                                break;
                            }
                            const auto& dtv = it->second;
                            
                            // create temp tree vector to store new copies into
                            DecayTreeVector DTV_temp;

                            // loop over daughter DecayTree's valid for the the particle combination
                            for (const auto& dt : dtv) {
                                if (pcs_has_pc(dt->freeAmplitude()->decayChannel()->particleCombinations(), *pc, d)) {

                                    // loop over parent DecayTree's
                                    for (const auto& DT : DTV) {
                                        // add copy of DT to DTV_temp
                                        DTV_temp.push_back(std::make_shared<DecayTree>(*DT));
                                        // add daughter DecayTree to it
                                        DTV_temp.back()->setDaughterDecayTree(d, dt);
                                    }

                                }
                            }

                            // repalce parent's DecayTreeVector with the new one with daughter d set into them
                            DTV = DTV_temp;
                        }
                        
                        // else if FinalStateParticle
                        else if (is_final_state_particle(Channels_.back()->daughters()[d])) {

                            // if channels daughter cannot accommodate daughter
                            if (!pcs_has_pc(Channels_.back()->daughters()[d]->particleCombinations(), *pc, d))
                                DTV.clear();

                        
                        }

                        // else was not DecayingParticle or FinalStateParticle
                        else
                            DTV.clear();
                        
                    } // ends loop over daughters

                    // if decay trees were created, add them into DecayTrees_
                    // if they aren't already present in it
                    if (!DTV.empty()) {
                        // this creates the empty DecayTreeVector if not yet created (hence the above conditional)
                        auto& dtv_M = DecayTrees_[two_M];
                        for (const auto& DT : DTV) {
                            // check that DT isn't already equal to one in dtv_M
                            if (std::none_of(dtv_M.begin(), dtv_M.end(), [&DT](const std::shared_ptr<DecayTree>& dt) {return *dt == *DT;}))
                                dtv_M.push_back(DT);
                        }
                    }

                } // ends loop over particle combinations of this decaying particle
            } // ends loop over spin projections of daughters
        } // ends loop over spin projection of parent
    } // ends loop over spin amplitude
    
    if (MassShape_)
        MassShape_->addDecayChannel(c);
}

//-------------------------
void DecayingParticle::addAllPossibleSpinAmplitudes(DecayChannel& dc, bool conserve_parity) const
{
    auto two_j = spins(dc.daughters());
    auto p = (conserve_parity) ? quantumNumbers().P() * parity(dc.daughters()) : 0;

    // create spin amplitudes
    // loop over possible S: |j1-j2| <= S <= (j1+j2)
    for (unsigned two_S = std::abs<int>(two_j[0] - two_j[1]); two_S <= two_j[0] + two_j[1]; two_S += 2)
        // loop over possible L: |J-s| <= L <= (J+s)
        for (unsigned L = std::abs<int>(quantumNumbers().twoJ() - two_S) / 2; L <= (quantumNumbers().twoJ() + two_S) / 2; ++L)
            // check parity conservation (also fulfilled if parity = 0)
            if (p * pow_negative_one(L) >= 0)
                // add SpinAmplitude retrieved from cache
                dc.addSpinAmplitude(const_cast<Model*>(dc.model())->spinAmplitudeCache()->spinAmplitude(quantumNumbers().twoJ(), two_j, L, two_S));
}

//-------------------------
std::shared_ptr<DecayChannel> DecayingParticle::addDecay(const ParticleVector& daughters, bool conserve_parity)
{
    auto dc = std::make_shared<DecayChannel>(daughters);
    addAllPossibleSpinAmplitudes(*dc, conserve_parity);
    addDecayChannel(dc);
    return dc;
}

//-------------------------
const Model* DecayingParticle::model() const
{
    return Channels_.empty() ? nullptr : Channels_[0]->model();
}

//-------------------------
void DecayingParticle::registerWithModel()
{
    for (auto& l_bw : BlattWeisskopfs_)
        l_bw.second->registerWithModel();

    for (auto& c : Channels_)
        c->registerWithModel();

    if (MassShape_)
        MassShape_->registerWithModel();
}

//-------------------------
void DecayingParticle::addParticleCombination(const ParticleCombination& pc)
{
    Particle::addParticleCombination(pc);

    // add also to all BlattWeiskopf barrier factors
    for (auto& l_bw : BlattWeisskopfs_)
        if (pc.daughters().size() == 2)
            l_bw.second->addParticleCombination(pc);

    // add to DecayChannels,
    // if DecayChannel contains particle combination with same content (without checking parent)
    // this is for the setting of ParticleCombination's with parents
    for (auto& dc : Channels_)
        if (std::any_of(dc->particleCombinations().begin(), dc->particleCombinations().end(), std::bind(&equal_down, pc.shared_from_this(), std::placeholders::_1)))
            dc->addParticleCombination(pc);

    if (MassShape_)
        MassShape_->addParticleCombination(pc);
}

//-------------------------
void DecayingParticle::pruneParticleCombinations()
{
    Particle::pruneParticleCombinations();

    for (auto& c : Channels_)
        c->pruneParticleCombinations();
}

//-------------------------
// helper function
void get_paddings(const DecayingParticle& dp, size_t& name_padding, size_t& sa_padding)
{
    name_padding = std::max(name_padding, dp.name().length());

    for (const auto& c : dp.channels()) {
        if (!c) continue;

        sa_padding = std::max(sa_padding, to_string(c->spinAmplitudes()).length());

        for (const auto& d : c->daughters()) {
            if (!d) continue;

            name_padding = std::max(name_padding, d->name().length());

            if (is_decaying_particle(*d))
                get_paddings(static_cast<const DecayingParticle&>(*d), name_padding, sa_padding);
        }
    }
}

//-------------------------
// helper function
std::string pad_right(const std::string& s, size_t len, char c = ' ')
{ return s + std::string(s.length() < len ? len - s.length() : 0, c); }

//-------------------------
std::string to_decay_string(const DecayingParticle& dp, unsigned level)
{
    static size_t name_padding = 0;
    static size_t sa_padding = 0;

    if (level == 0)
        get_paddings(dp, name_padding, sa_padding);

    std::string s;
    unsigned i = 0;
    for (const auto& c : dp.channels()) {
        if (!c) continue;

        if (i++ > 0)
            s += pad_right("\n", level * (3 * name_padding + 8 + sa_padding));

        s += pad_right(dp.name(), name_padding) + " -> ";
        for (const auto d : c->daughters())
            if (d)
                s += pad_right(d->name(), name_padding + 1);

        s += pad_right(to_string(c->spinAmplitudes()), sa_padding);

        for (const auto& d : filter(c->daughters(), [](const std::shared_ptr<Particle>& p){return p and is_decaying_particle(*p);}))
            s += ", " + to_decay_string(static_cast<const DecayingParticle&>(*d), level + 1);
    }

    return s;
}

//-------------------------
void DecayingParticle::modifyDecayTree(DecayTree& dt)
{
    if (!dt.freeAmplitude())
        throw exceptions::Exception("DecayTree has nullptr free amplitude", "DecayingParticle::modifyDecayTree");

    if (!dt.freeAmplitude()->spinAmplitude())
        throw exceptions::Exception("FreeAmplitude's SpinAmplitude is nullptr", "DecayingParticle::modifyDecayTree");

    // find BlattWeisskopf object
    if (dt.freeAmplitude()->spinAmplitude()->L() > 0) {
        auto l_bw = BlattWeisskopfs_.find(dt.freeAmplitude()->spinAmplitude()->L());
        if (l_bw == BlattWeisskopfs_.end())
            throw exceptions::Exception("No Blatt-Weisskopf factor found for L = "
                                        + std::to_string(dt.freeAmplitude()->spinAmplitude()->L()),
                                        "DecayingParticle::modifyDecayTree");

        if (!l_bw->second)
            throw exceptions::Exception("BlattWeisskopf is nullptr", "DecayingParticle::modifyDecayTree");

        // Add BlattWeisskopf object
        dt.addAmplitudeComponent(*l_bw->second);
    }

    if (MassShape_)
        dt.addAmplitudeComponent(*MassShape_);
}

//-------------------------
ParticleSet particles(DecayingParticle& dp)
{
    ParticleSet S = {dp.shared_from_this()};
    for (const auto& dc : dp.channels()) {
        auto s = particles(*dc);
        S.insert(s.begin(), s.end());
    }
    return S;
}

//-------------------------
DecayTreeSet decay_trees(const DecayingParticle& dp)
{
    DecayTreeSet S;
    for (const auto& m_dtv : dp.decayTrees())
        S.insert(m_dtv.second.begin(), m_dtv.second.end());
    return S;
}

//-------------------------
FreeAmplitudeSet free_amplitudes(const DecayingParticle& dp)
{
    FreeAmplitudeSet S;
    for (const auto& m_dtv : dp.decayTrees()) {
        auto s = free_amplitudes(m_dtv.second);
        S.insert(s.begin(), s.end());
    }
    return S;
}

}
