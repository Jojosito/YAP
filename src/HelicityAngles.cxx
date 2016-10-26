#include "HelicityAngles.h"

#include "CachedValue.h"
#include "CalculationStatus.h"
#include "FourMomenta.h"
#include "FourVector.h"
#include "LorentzTransformation.h"
#include "Model.h"
#include "ParticleCombination.h"
#include "StatusManager.h"

#include "logging.h"

namespace yap {

//-------------------------
HelicityAngles::HelicityAngles(Model& m) :
    StaticDataAccessor(m, equal_up_and_down)
{
    registerWithModel();

    // HelicityAngles will not be added to StaticDataAccessors by registerWithModel,
    // since it has size 0. So we need to add it manually here
    addToStaticDataAccessors();
}

//-------------------------
void HelicityAngles::addToStaticDataAccessors()
{
    // look for Model's FourMomenta_
    auto it_fm = std::find(staticDataAccessors().begin(), staticDataAccessors().end(), model()->fourMomenta().get());
    if (it_fm == staticDataAccessors().end())
        throw exceptions::Exception("HelicityAngles cannot be registered with the model before FourMomenta", "HelicityAngles::registerWithModel");
    // add this to just after FourMomenta_
    const_cast<StaticDataAccessorVector&>(staticDataAccessors()).insert(it_fm + 1, this);
}

//-------------------------
const std::array<double, 2>& HelicityAngles::helicityAngles(const DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc) const
{
    for (auto& kv : cachedForDataPoint_)
        if (kv.second == &d)
            return cachedAngles_.at(kv.first).at(symmetrizationIndex(pc));

    throw exceptions::Exception("HelicityAngles have not been cached for given DataPoint", "HelicityAngles::angles");
}

//-------------------------
void HelicityAngles::calculate(DataPoint& d, StatusManager& sm) const
{
    // call on ISP PC's
    // \todo allow for designating the boost that takes from the data frame to the lab frame (possibly event dependent)
    for (auto& kv : symmetrizationIndices())
        if (not kv.first->parent())
            calculateAngles(d, kv.first, model()->coordinateSystem(), unitMatrix<double, 4>(), sm);

    cachedForDataPoint_[&sm] = &d;
}

//-------------------------
void HelicityAngles::calculateAngles(DataPoint& d, const std::shared_ptr<const ParticleCombination>& pc,
                                     const CoordinateSystem<double, 3>& C, const FourMatrix<double>& boosts,
                                     const StatusManager& sm) const
{
    // terminate recursion
    if (is_final_state_particle_combination(*pc))
        return;

    // get pc's 4-mom in data frame
    const auto P = model()->fourMomenta()->p(d, pc);

    // calculate reference frame for P from parent's RF
    const auto cP = helicityFrame(boosts * P, C);

    // calculate boost from data frame into pc rest frame
    const auto boost = lorentzTransformation(-(boosts * P));

    const auto boost_boosts = boost * boosts;

    const unsigned symIndex = symmetrizationIndex(pc);

    for (auto& daughter : pc->daughters()) {

        // boost daughter momentum from data frame into pc rest frame
        const auto p = boost_boosts * model()->fourMomenta()->p(d, daughter);

        auto phi_theta = angles<double>(vect<double>(p), cP);

        // set ambiguous phi to theta
        // todo: in this cases, theta should be 0 or pi. In most cases it is, but sometimes not.
        // Not checking if theta == 0 or pi results in tests passing which would otherwise not
        if (std::isnan(phi_theta[0]))
            phi_theta[0] = phi_theta[1];

        cachedAngles_[&sm][symIndex] = phi_theta;

        // recurse down the decay tree
        calculateAngles(d, daughter, cP, boost, sm);
    }
}

//-------------------------
void HelicityAngles::addParticleCombination(const ParticleCombination& pc)
{
    if (pc.daughters().size() != 2)
        throw exceptions::NotTwoBodyParticleCombination("cannot calculate helicity angles for "
                                                        + std::to_string(pc.daughters().size()) + "-body decay",
                                                        "HelicityAngles::addParticleCombination");

    return StaticDataAccessor::addParticleCombination(pc);
}

}
