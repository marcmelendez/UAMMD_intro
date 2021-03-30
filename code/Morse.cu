# include "uammd.cuh"
# include "utils/InitialConditions.cuh"
# include "Interactor/Potential/Potential.cuh"
# include "Interactor/NeighbourList/CellList.cuh"
# include "Interactor/PairForces.cuh"
# include "Integrator/VerletNVE.cuh"

using namespace uammd;
using std::make_shared;
using std::endl;

struct Morse {
  real De, a, r0, rc;

  Morse(real i_De, real i_a, real i_r0, real i_rc):
       De(i_De), a(i_a), r0(i_r0), rc(i_rc){} //!
  real getCutOff() { return rc; } //!
  struct ForceEnergy {
    real4 * force;
    real * energy;
    Box box;
    real De, a, r0, rc;

    ForceEnergy(Box i_box, real i_rc,
                real4 * i_force, real * i_energy,
                real i_De, real i_a, real i_r0):
                box(i_box), rc(i_rc),
                force(i_force), energy(i_energy),
                De(i_De), a(i_a), r0(i_r0){} //!
    __device__ real4 compute(real4 ri, real4 rj){
      const real3 rij = box.apply_pbc(make_real3(rj)-make_real3(ri));
      const real r2 = dot(rij, rij);
      const real r = sqrtf(r2);
      if(r2 > 0 and r < rc){
        const real expar = exp(-a*(r - r0));
        const real oneminusexpar = real(1.0) - expar;
        return make_real4(2.0*De*a*(oneminusexpar*expar/r)*rij,
                          De*(oneminusexpar*oneminusexpar - real(1.0)));
      }
      return real4();
    }

    __device__ void set(int id, real4 total){
      force[id] += make_real4(total.x, total.y, total.z, 0);
      energy[id] += real(0.5)*total.w;
    }
  }; //!
  ForceEnergy getForceTransverser(Box box, std::shared_ptr<ParticleData> sys){
    auto force = sys->getForce(access::location::gpu,
                               access::mode::readwrite).raw();
    auto energy = sys->getEnergy(access::location::gpu,
                                 access::mode::readwrite).raw();
    return ForceEnergy(box, rc, force, energy, De, a, r0);
  }
}; //!

int main(int argc, char *argv[]){

  auto sys = make_shared<System>(argc, argv);

  int numberOfParticles = 100000;
  auto particles
    = make_shared<ParticleData>(numberOfParticles, sys);//!

  real L = 128;

  Box box(make_real3(L, L, L));
  bool periodicityX = true, periodicityY = true,
       periodicityZ = true;
  box.setPeriodicity(periodicityX, periodicityY,
                     periodicityZ); //!
  {
    auto position
      = particles->getPos(access::location::cpu,
                          access::mode::write);

    auto initial =  initLattice(box.boxSize,
                                numberOfParticles, sc);

    std::copy(initial.begin(), initial.end(), position.begin());
  } //!

  using Verlet = VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = 0.01;
  VerletParams.initVelocities = true;
  VerletParams.energy = 1.0;//!

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);//!

  real De = 1.0;
  real a = 2.0;
  real r0 = 1.0;
  real rc = 6.5*r0;

  auto MorsePotential = make_shared<Morse>(De, a, r0, rc);
  {
    using MorseForces = PairForces<Morse>;
    MorseForces::Parameters interactionParams;
    interactionParams.box = box;

    auto interaction
      = make_shared<MorseForces>(particles, sys,
                              interactionParams,
                              MorsePotential);

    integrator->addInteractor(interaction);
  } //!

  std::string outputFile = "Morse.dat";
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
        out<<box.apply_pbc(make_real3(position[index[id]]))<<endl; //!
    }
  } //!

  sys->finish();

  return 0;
}
