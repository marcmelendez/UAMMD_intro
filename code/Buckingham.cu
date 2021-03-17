# include "uammd.cuh"
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

  Buckingham(real i_A, real i_B, real i_C, real i_rc):
             A(i_A), B(i_B), C(i_C), rc(i_rc){}

  real getCutOff() { return rc; }

  struct ForceEnergy {
    real4 * force;
    real * energy;
    Box box;

    real A, B, C, rc;

    ForceEnergy(Box i_box, real i_rc,
                real4 * i_force, real * i_energy,
                real i_A, real i_B, real i_C):
                box(i_box), rc(i_rc),
                force(i_force), energy(i_energy),
                A(i_A), B(i_B), C(i_C){}

    __device__ real4 compute(real4 ri, real4 rj){
      const real3 rij = box.apply_pbc(make_real3(rj)-make_real3(ri));
      const real r2 = dot(rij, rij);
      const real r = sqrtf(r2);
      if(r2 > 0 and r < rc) {
        const real expBr = exp(-B*r);
        const real invr2 = real(1.0)/r2;
        const real invr6 = invr2*invr2*invr2;
        return make_real4((-A*B*expBr + 6*C*invr6*invr2)*rij,
                          A*expBr - C*invr6);
      }
      return real4();
    }

    __device__ void set(int id, real4 total){
      force[id] += make_real4(total.x, total.y, total.z, 0);
      energy[id] += total.w;
    }
  };

  struct Virial {
    real4 * position;
    real4 * virial;
    Box box;

    real A, B, C, rc;

    Virial(Box i_box, real i_rc,
           real4 * i_position, real4 * i_virial,
           real i_A, real i_B, real i_C):
           box(i_box), rc(i_rc),
           position(i_position), virial(i_virial),
           A(i_A), B(i_B), C(i_C){}

    __device__ real compute(real4 ri, real4 rj){
      const real3 rij = box.apply_pbc(make_real3(rj)-make_real3(ri));
      const real r2 = dot(rij, rij);
      real3 Fij = make_real3(ForceEnergy.compute(ri, rj));
      return dot(Fij, rij);
    }

    __device__ void set(int id, real total){
      virial[id].x += total;
    }
  };

  ForceEnergy getForceTransverser(Box box, std::shared_ptr<ParticleData> sys){
    auto force = sys->getForce(access::location::gpu,
                               access::mode::readwrite).raw();
    auto energy = sys->getEnergy(access::location::gpu,
                                 access::mode::readwrite).raw();
    return ForceEnergy(box, rc, force, energy, A, B, C);
  }

  Virial getComputeTransverser(Box box, std::shared_ptr<ParticleData> sys) {
    auto position = sys->getPos(access::location::gpu,
                                access::mode::read).raw();
    auto virial = sys->getVirial(access:location::gpu,
                                 access::mode::readwrite).raw();
    return Virial(box, rc, position, virial, A, B, C);
  }
};

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

  using Verlet = VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = 0.01;
  VerletParams.initVelocities = true;
  VerletParams.energy = 1.0;

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);

  real A = 37101;
  real B = 4.4148;
  real C = 113.09;
  real rc = 12;
  auto BuckinghamPotential = make_shared<Buckingham>(A, B, C, rc);
  {
    using BuckinghamForces = PairForces<Buckingham>;
    BuckinghamForces::Parameters interactionParams;
    interactionParams.box = box;

    auto interaction
      = make_shared<BuckinghamForces>(particles, sys,
                                      interactionParams,
                                      BuckinghamPotential);

    integrator->addInteractor(interaction);
  }

  std::string outputFile = "helium.dat";
  std::ofstream out(outputFile);

  int numberOfSteps = 1000;
  int printEverynSteps = 100;

  for(int step = 0; step < numberOfSteps; ++step) {
    integrator->forwardTime();

    if(printEverynSteps > 0
       and step % printEverynSteps == 1) {
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
