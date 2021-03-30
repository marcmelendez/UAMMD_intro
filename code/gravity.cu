# include "uammd.cuh"
# include "utils/InitialConditions.cuh"
# include "Interactor/Potential/Potential.cuh"
# include "Interactor/NeighbourList/CellList.cuh"
# include "Interactor/PairForces.cuh"
# include "Integrator/VerletNVE.cuh"

using namespace uammd;
using std::make_shared;
using std::endl;

struct gravity {
  real G;

  gravity(real i_G) : G(i_G) {}

  real getCutOff() { return std::numeric_limits<real>::infinity(); }

  struct ForceEnergy {
    real4 * force;
    real * energy;
    real * mass;
    real G;

    ForceEnergy(real4 * i_force, real * i_energy,
                real * i_mass, real i_G):
                force(i_force), energy(i_energy),
                mass(i_mass), G(i_G){}

    __device__ real getInfo(int index) {
      return mass[index];
    }

    __device__ real4 compute(real4 ri, real4 rj, real mi, real mj){
      const real3 rij = make_real3(rj)-make_real3(ri);
      const real r2 = dot(rij, rij);
      const real r = sqrtf(r2);
      if(r2 > 0){
        return make_real4(G*(mi*mj/(r2*r))*rij,
                          -G*mi*mj/r);
      }
      return real4();
    }

    __device__ void set(int id, real4 total){
      force[id] += make_real4(total.x, total.y, total.z, 0);
      energy[id] += 0.5*total.w;
    }
  };

  ForceEnergy getForceTransverser(Box box, std::shared_ptr<ParticleData> sys){
    auto force = sys->getForce(access::location::gpu,
                               access::mode::readwrite).raw();
    auto energy = sys->getEnergy(access::location::gpu,
                                 access::mode::readwrite).raw();
    auto mass = sys->getMass(access::location::gpu,
                             access::mode::read).raw();
    return ForceEnergy(force, energy, mass, G);
  }
};//!

int main(int argc, char *argv[]){

  auto sys = make_shared<System>(argc, argv);
  int numberOfParticles;

  std::ifstream in("solarSystem.dat");
  if(!std::ifstream("solarSystem.dat").good()) {
    sys->log<System::CRITICAL>("File solarSystem.dat not found.");
  }

  in>>numberOfParticles;

  auto particles
    = make_shared<ParticleData>(numberOfParticles, sys);

  {
    auto position
      = particles->getPos(access::location::cpu,
                          access::mode::write);
    auto velocity
      = particles->getVel(access::location::cpu,
                          access::mode::write);
    auto mass
      = particles->getMass(access::location::cpu,
                           access::mode::write);

    sys->log<System::MESSAGE>("Reading data for astronomical objects:");
    std::string objectName;

    for(int i = 0; i < numberOfParticles; ++i) {
      in>>position[i].x>>position[i].y>>position[i].z
        >>velocity[i].x>>velocity[i].y>>velocity[i].z
        >>mass[i]>>objectName;
      position[i].w = real(0.0);
      sys->log<System::MESSAGE>(objectName);
    }
  }//!

  using Verlet = VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = 0.01;
  VerletParams.initVelocities = false;

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);

  real G = 8.888e-10;
  auto GPotential = make_shared<gravity>(G);
  {
    using GForces = PairForces<gravity>;

    auto interaction
      = make_shared<GForces>(particles, sys, GPotential);

    integrator->addInteractor(interaction);
  }//!

  std::string outputFile = "gravity.dat";
  std::ofstream out(outputFile);

  int numberOfSteps = 1000000;
  int printEverynSteps = 1000;

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
        out<<make_real3(position[index[id]])<<endl;
    }
  }

  sys->finish();

  return 0;
}//!
