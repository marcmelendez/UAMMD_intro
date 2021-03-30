# include "uammd.cuh"
# include "utils/InputFile.h"
# include "utils/InitialConditions.cuh"
# include "Interactor/Potential/Potential.cuh"
# include "Interactor/NeighbourList/CellList.cuh"
# include "Interactor/PairForces.cuh"
# include "Integrator/VerletNVE.cuh"

using namespace uammd;
using std::make_shared;
using std::endl;

struct Buckingham {
  real A, B, C;
  real rc;
  real shift;

  Buckingham(real i_A, real i_B, real i_C, real i_rc):
             A(i_A), B(i_B), C(i_C), rc(i_rc) {
    real invr2 = real(1.0)/(rc*rc);
    real invr6 = invr2*invr2*invr2;
    shift = A*expf(-B*rc) - C*invr6;
  }

  real getCutOff() { return rc; } //!

  struct ForceEnergy {
    real4 * force;
    real * energy;
    Box box;

    real A, B, C, rc, shift;

    ForceEnergy(Box i_box, real i_rc,
                real4 * i_force, real * i_energy,
                real i_A, real i_B, real i_C, real i_shift):
                box(i_box), rc(i_rc),
                force(i_force), energy(i_energy),
                A(i_A), B(i_B), C(i_C), shift(i_shift) {}

    __device__ real4 compute(real4 ri, real4 rj){
      const real3 rij = box.apply_pbc(make_real3(rj)-make_real3(ri));
      const real r2 = dot(rij, rij);
      if(r2 > 0 and r2 < rc*rc) {
        const real r = sqrtf(r2);
        const real AexpmBr = A*expf(-B*r);
        const real invr2 = real(1.0)/r2;
        const real invr6 = invr2*invr2*invr2;
        return make_real4((-B*AexpmBr/r + 6*C*invr6*invr2)*rij,
                          AexpmBr - C*invr6 + shift);
      }
      return real4();
    }

    __device__ void set(int id, real4 total){
      force[id] += make_real4(total.x, total.y, total.z, 0);
      energy[id] += real(0.5)*total.w;
    }
  }; //!

  struct Virial {
    real4 * position;
    real * virial;
    Box box;

    real A, B, C, rc;

    Virial(Box i_box, real i_rc,
           real4 * i_position, real * i_virial,
           real i_A, real i_B, real i_C):
           box(i_box), rc(i_rc),
           position(i_position), virial(i_virial),
           A(i_A), B(i_B), C(i_C){}

    __device__ real compute(real4 ri, real4 rj){
      const real3 rij = box.apply_pbc(make_real3(rj)-make_real3(ri));
      const real r2 = dot(rij, rij);
      if(r2 > 0 and r2 < rc*rc) {
        const real r = sqrtf(r2);
        const real AexpmBr = A*expf(-B*r);
        const real invr2 = real(1.0)/r2;
        const real invr6 = invr2*invr2*invr2;
        return B*AexpmBr*r - 6*C*invr6;
      } else {
        return real(0.0);
      }
    }

    __device__ void set(int id, real total){
      virial[id] += total;
    }
  }; //!

  ForceEnergy getForceEnergyTransverser(Box box, std::shared_ptr<ParticleData> sys){
    auto force = sys->getForce(access::location::gpu,
                               access::mode::readwrite).raw();
    auto energy = sys->getEnergy(access::location::gpu,
                                 access::mode::readwrite).raw();
    return ForceEnergy(box, rc, force, energy, A, B, C, shift);
  }

  Virial getComputeTransverser(Box box, std::shared_ptr<ParticleData> sys) {
    auto position = sys->getPos(access::location::gpu,
                                access::mode::read).raw();
    auto virial = sys->getVirial(access::location::gpu,
                                 access::mode::readwrite).raw();
    return Virial(box, rc, position, virial, A, B, C);
  }
}; //!

double getTotalEnergy(std::shared_ptr<Integrator> integrator,
                      std::shared_ptr<ParticleData> particles){
  {
    auto energy
      = particles->getEnergy(access::location::cpu,
                             access::mode::write);
    std::fill(energy.begin(), energy.end(), real(0.0));
  }

  double totalEnergy = 0; //!
  integrator->sumEnergy(); //!
  for(auto interactor: integrator->getInteractors()){
    interactor->sumEnergy();
  } //!
  {
    auto energy
      = particles->getEnergy(access::location::cpu,
                             access::mode::read);
    for(int i = 0; i < particles->getNumParticles(); ++i) {
      totalEnergy += energy[i];
    }
  }
  return totalEnergy;
} //!

real3 getTotalMomentum(std::shared_ptr<ParticleData> particles){
    auto velocity
      = particles->getVel(access::location::cpu,
                          access::mode::read);
    auto mass
      = particles->getMass(access::location::cpu,
                           access::mode::read);

    real3 totalMomentum = make_real3(0.0, 0.0, 0.0);

    for(int i = 0; i < particles->getNumParticles(); ++i) {
      totalMomentum += mass[i]*velocity[i];
    }

  return totalMomentum;
} //!

double getThermalEnergy(std::shared_ptr<ParticleData> particles){
  int N = particles->getNumParticles();
  auto velocity
    = particles->getVel(access::location::cpu,
                          access::mode::read);
  auto mass
    = particles->getMass(access::location::cpu,
                           access::mode::read);

  real3 Vcm = make_real3(0.0, 0.0, 0.0);
  double M = real(0.0);

  for(int i = 0; i < N; ++i) {
    Vcm += mass[i]*velocity[i];
    M += mass[i];
  }
  Vcm /= M; //!
  double kineticEnergy = real(0.0);
  for(int i = 0; i < N; ++i) {
    kineticEnergy
     += real(0.5)*mass[i]*dot(velocity[i] - Vcm, velocity[i] - Vcm);
  }

  return real(2.0/(3.0*N))*kineticEnergy;
}//!

double getTotalVirial(std::shared_ptr<ParticleData> particles){
  auto virial
    = particles->getVirial(access::location::cpu,
                           access::mode::read);

  double totalVirial = 0.0;

  for(int i = 0; i < particles->getNumParticles(); ++i) {
    totalVirial += virial[i];
  }

  return totalVirial;
}

int main(int argc, char * argv[])
{

  auto sys = make_shared<System>(argc, argv);

  int numberOfParticles = 10000;
  auto particles
    = make_shared<ParticleData>(numberOfParticles, sys);

  real L = 64;

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

  {
    auto mass
      = particles->getMass(access::location::cpu,
                           access::mode::write);

    std::fill(mass.begin(), mass.end(), real(1.0));
  }

  using Verlet = VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = 0.001;
  VerletParams.initVelocities = true;
  VerletParams.energy = 2.0;

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);

  real A = 37101;
  real B = 4.4148;
  real C = 113.09;
  real rc = 6.0;

  auto BuckinghamPotential = make_shared<Buckingham>(A, B, C, rc);

  using BuckinghamForces = PairForces<Buckingham>;
  BuckinghamForces::Parameters interactionParams;
  interactionParams.box = box;

  auto interaction
    = make_shared<BuckinghamForces>(particles, sys,
                                      interactionParams,
                                      BuckinghamPotential);

  integrator->addInteractor(interaction);

  std::string outputFile = "helium.dat";
  std::ofstream out(outputFile);

  std::string macroFile = "heliumMacro.dat";
  std::ofstream macro(macroFile);

  int numberOfSteps = 10000;
  int printEverynSteps = 100;

  for(int step = 0; step < numberOfSteps; ++step) {
    integrator->forwardTime();

    if(printEverynSteps > 0
       and step % printEverynSteps == 1) {
      {
        auto virial
          = particles->getVirial(access::location::cpu,
                                 access::mode::write);

        std::fill(virial.begin(), virial.end(), real(0.0));
      }

      interaction->compute(0);

      real thermalEnergy = getThermalEnergy(particles);
      real totalVirial = getTotalVirial(particles);
      real volume = box.boxSize.x*box.boxSize.y*box.boxSize.z;
      real pressure = particles->getNumParticles()*thermalEnergy/volume
                      + totalVirial/(6*volume);

      auto position
        = particles->getPos(access::location::cpu,
                            access::mode::read);
      const int * index = particles->getIdOrderedIndices(access::location::cpu);

      out<<endl;
      for(int id = 0; id < numberOfParticles; ++id)
        out<<box.apply_pbc(make_real3(position[index[id]]))<<endl;

      macro<<step*VerletParams.dt<<" ";
      macro<<getTotalEnergy(integrator, particles)<<" ";
      macro<<getTotalMomentum(particles)<<" ";
      macro<<thermalEnergy<<" "<<pressure<<endl;
    }
  }

  sys->finish();

  return 0;
}


