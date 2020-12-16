# include "uammd.cuh"
# include "utils/InitialConditions.cuh"
# include "Interactor/Potential/Potential.cuh"
# include "Interactor/NeighbourList/CellList.cuh"
# include "Interactor/PairForces.cuh"
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

    Box box(make_real3(L, L, L));
    bool periodicityX = true, periodicityY = true,
         periodicityZ = true;
    box.setPeriodicity(periodicityX, periodicityY,
                       periodicityZ);
  {
    auto position
      = particles->getPos(access::location::cpu,
                          access::mode::write);

    auto initial =  initLattice(box.boxSize,
                                numberOfParticles, sc);

    std::copy(initial.begin(), initial.end(), position.begin());
  }

  using Verlet = VerletNVE::VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = 0.01;
  VerletParams.initVelocities=true;
  VerletParams.energy = 1.0;

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);

  auto LJPotential = make_shared<Potential::LJ>(sys);
  {
    Potential::LJ::InputPairParameters LJParams;
    LJParams.epsilon = 1.0;
    LJParams.sigma = 1.0;
    LJParams.cutOff = 2.5*LJParams.sigma;
    LJParams.shift = false;
    LJPotential->setPotParameters(0, 0, LJParams);
  }

  {
    using LJForces = PairForces<Potential::LJ>;
    LJForces::Parameters interactionParams;
    interactionParams.box = box;

    auto allParticles
      = make_shared<ParticleGroup>(particles, sys, "All");

    auto interaction
      = make_shared<LJForces>(particles, allParticles,
                              sys, interactionParams,
                                LJPotential);

    integrator->addInteractor(interaction);
  }

  std::string outputFile = "Lennard-Jones.dat";
  std::ofstream out(outputFile);

  int numberOfSteps = 1000;
  int printEverynSteps = 100;

  for(int step = 0; step < numberOfSteps; ++step) {
    integrator->forwardTime();

    if(printEverynSteps > 0
       && step % printEverynSteps == 1) {
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
