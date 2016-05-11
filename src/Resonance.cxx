#include "Resonance.h"

#include "Exceptions.h"
#include "logging.h"
#include "MassShape.h"
#include "StatusManager.h"

namespace yap {

//-------------------------
Resonance::Resonance(const QuantumNumbers& q, double mass, std::string name, double radialSize, std::shared_ptr<MassShape> massShape) :
    DecayingParticle(q, mass, name, radialSize),
    MassShape_(massShape)
{
    if (!MassShape_)
        throw exceptions::Exception("MassShape unset", "Resonance::Resonance");

    MassShape_->setResonance(this);
}

//-------------------------
std::complex<double> Resonance::amplitude(DataPoint& d, const std::shared_ptr<ParticleCombination>& pc, int two_m, StatusManager& sm) const
{
    return DecayingParticle::amplitude(d, pc, two_m, sm) * MassShape_->amplitude(d, pc, sm);
}

//-------------------------
bool Resonance::consistent() const
{
    bool C = DecayingParticle::consistent();

    if (!MassShape_->consistent()) {
        FLOG(ERROR) << "MassShape not consistent";
        C &= false;
    }
    if (MassShape_->resonance() != this) {
        FLOG(ERROR) << "MassShape's resonance is not this";
        C &= false;
    }

    return C;
}

//-------------------------
void Resonance::addToModel()
{
    DecayingParticle::addToModel();
    MassShape_->addToModel();
}

//-------------------------
unsigned Resonance::addParticleCombination(std::shared_ptr<ParticleCombination> c)
{
    unsigned index = DecayingParticle::addParticleCombination(c);
    MassShape_->addParticleCombination(c);
    return index;
}

}
