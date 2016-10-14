#include "MassShapeWithNominalMass.h"

#include "CachedValue.h"
#include "DataPartition.h"
#include "Filters.h"
#include "FourMomenta.h"
#include "logging.h"
#include "Model.h"
#include "Parameter.h"
#include "ParticleFactory.h"
#include "Resonance.h"

namespace yap {

//-------------------------
MassShapeWithNominalMass::MassShapeWithNominalMass(double m) :
    MassShape(),
    Mass_(std::make_shared<RealParameter>(m))
{
    addParameter(Mass_);
}

//-------------------------
void MassShapeWithNominalMass::setParameters(const ParticleTableEntry& entry)
{
    if (Mass_->value() < 0)
        *Mass_ = entry.mass();
}

//-------------------------
RealParameter& mass_parameter(Particle& p)
{
    if (!is_resonance(p) or !has_mass()(p))
        throw exceptions::Exception("Particle has no mass parameter", "mass_parameter");
    auto m = std::static_pointer_cast<MassShapeWithNominalMass>(static_cast<Resonance*>(&p)->massShape())->mass();
    if (!m)
        throw exceptions::Exception("Particle's accessible mass parameter is nullptr", "mass_parameter");
    return *m;
}

}




