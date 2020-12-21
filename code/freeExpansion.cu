# include <iostream>
# include "uammd.cuh"
# include "Integrator/VerletNVE.cuh"

using namespace uammd;
using std::make_shared;
using std::endl;

int main(int argc, char *argv[]){

  auto sys = make_shared<System>(argc, argv);

  int numberOfParticles = 100000;
  auto particles
    = make_shared<ParticleData>(numberOfParticles, sys);

  real L = 128;
  Box box(make_real3(L, L, std::numeric_limits<real>::infinity()));
  bool periodicityX = true, periodicityY = true,
       periodicityZ = false;
  box.setPeriodicity(periodicityX, periodicityY,
                     periodicityZ);

  {
    auto position
      = particles->getPos(access::location::cpu,
                          access::mode::write);

    for(int i = 0; i < numberOfParticles; ++i)
      position[i]
        = make_real4(sys->rng().uniform3(-0.5, 0.5), 0)*L;
  }

  using Verlet = VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = 0.01;
  VerletParams.initVelocities=true;
  VerletParams.energy = 1.0;

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);

  std::string outputFile = "freeExpansion.dat";
  std::ofstream out(outputFile);

  int numberOfSteps = 1000;
  int printEverynSteps = 100;

  for(int step = 0; step < numberOfSteps; ++step) {
    integrator->forwardTime();

    if(printEverynSteps > 0
       and step % printEverynSteps == 1) {
      /* ... Output particle positions ... */
      auto position
        = particles->getPos(access::location::cpu,
                            access::mode::read);
      const int * index = particles->getIdOrderedIndices(access::location::cpu);

      out<<endl;
      for(int id = 0; id < numberOfParticles; ++id)
        out<<box.apply_pbc(make_real3(position[index[id]]))<<endl;
    }
  }

  sys->finish();

  return 0;
}
