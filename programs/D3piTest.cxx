#include "ParticleFactory.h"

#include <assert.h>
#include <iostream>

#include "logging.h"
INITIALIZE_EASYLOGGINGPP

int main( int argc, char** argv)
{
  yap::ParticleFactory factory;

  // final state particles
  yap::FinalStateParticle* piPlus = factory.createFinalStateParticle(211);
  yap::FinalStateParticle* piMinus = factory.createFinalStateParticle(-211);

  yap::InitialState = InitialStateParticle(QuantumNumbers(0, -1, -1, 1, 0)

  double radialSize = 1.;
  unsigned L = 0;
  yap::Resonance* rho = factory.createResonanceBreitWigner(113, radialSize);
  factory.createChannel(rho, piPlus, piMinus, L);



  assert(rho->consistent());
}