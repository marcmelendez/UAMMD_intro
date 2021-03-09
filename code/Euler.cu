# include "uammd.cuh"
# include "utils/InitialConditions.cuh"
# include "Interactor/Potential/Potential.cuh"
# include "Interactor/NeighbourList/CellList.cuh"
# include "Interactor/PairForces.cuh"
# include "Integrator/Integrator.cuh"

using namespace uammd;
using std::make_shared;
using std::endl;

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

__global__ void EulerStep(int N, real dt,
                          real4 * position,
                          real3 * velocity,
                          real * mass,
                          real4 * force){
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i > N) return;

  position[i] += make_real4(velocity[i]*dt, real(0.0));
  velocity[i] += make_real3(force[i].x*dt/mass[i],
                            force[i].y*dt/mass[i],
                            force[i].z*dt/mass[i]);

  return;
} //!

__global__ void kineticEnergy(int N,
                              real3 * velocity,
                              real * mass,
                              real * energy){
  const int i = blockIdx.x*blockDim.x + threadIdx.x;
  if(i > N) return;

  energy[i] += real(0.5)*mass[i]*dot(velocity[i], velocity[i]);

  return;
} //!

class Euler : public Integrator {
  real dt;

  public:
    struct Parameters {
      real dt;
    };

    Euler(std::shared_ptr<ParticleData> particles,
          std::shared_ptr<System> sys,
          Parameters params) :
          Integrator(particles,
                     std::make_shared<ParticleGroup>(particles, sys, "ALL"),
                     sys, "Euler"),
          dt(params.dt){} //!
    virtual void forwardTime() override{
      auto particles = pd;

      {
        auto force
          = particles->getForce(access::location::cpu,
                                access::mode::write);

        std::fill(force.begin(), force.end(),
                  make_real4(0.0, 0.0, 0.0, 0.0));
      } //!
      for(auto interaction: interactors) {
        interaction->sumForce(0);
      } //!
      auto position
        = particles->getPos(access::location::gpu,
                            access::mode::readwrite);

      auto velocity
        = particles->getVel(access::location::gpu,
                            access::mode::readwrite);

      auto mass
        = particles->getMass(access::location::gpu,
                             access::mode::read);

      auto force
        = particles->getForce(access::location::gpu,
                              access::mode::read);

      int numberOfParticles = particles->getNumParticles();
      int Nthreads = 128;
      int Nblocks = numberOfParticles/Nthreads
                    + ((numberOfParticles%Nthreads)?1:0);

      EulerStep<<<Nblocks,Nthreads,0,0>>>(numberOfParticles, dt,
                                          position.raw(),
                                          velocity.raw(),
                                          mass.raw(),
                                          force.raw());
    } //!
    virtual real sumEnergy() override{
      auto particles = pd;

      auto velocity
        = particles->getVel(access::location::gpu,
                            access::mode::read);

      auto mass
        = particles->getMass(access::location::gpu,
                             access::mode::read);

      auto energy
        = particles->getEnergy(access::location::gpu,
                               access::mode::write);

      int numberOfParticles = particles->getNumParticles();
      int Nthreads = 128;
      int Nblocks = numberOfParticles/Nthreads
                    + ((numberOfParticles%Nthreads)?1:0);

      kineticEnergy<<<Nblocks,Nthreads,0,0>>>(numberOfParticles,
                                              velocity.raw(),
                                              mass.raw(),
                                              energy.raw());
      return 0;
    }
};

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

  {
    auto velocity
      = particles->getVel(access::location::cpu,
                          access::mode::write);

    for(int i = 0; i < numberOfParticles; ++i) {
      velocity[i] = make_real3(sys->rng().gaussian(0,1),
                               sys->rng().gaussian(0,1),
                               sys->rng().gaussian(0,1));
    }
  }

  {
    auto mass
      = particles->getMass(access::location::cpu,
                           access::mode::write);

    std::fill(mass.begin(), mass.end(), real(1.0));
  } //!

  Euler::Parameters simParams;
  simParams.dt = 0.001;

  auto integrator
    = make_shared<Euler>(particles, sys, simParams); //!

  auto LJPotential = make_shared<Potential::LJ>(sys);
  {
    Potential::LJ::InputPairParameters LJParams;
    LJParams.epsilon = 1.0;
    LJParams.sigma = 1.0;
    LJParams.cutOff = 2.5*LJParams.sigma;
    LJParams.shift = true;
    LJPotential->setPotParameters(0, 0, LJParams);
  } //!

  {
    using LJForces = PairForces<Potential::LJ>;
    LJForces::Parameters interactionParams;
    interactionParams.box = box;

    auto interaction
      = make_shared<LJForces>(particles, sys,
                              interactionParams,
                              LJPotential);

    integrator->addInteractor(interaction);
  } //!

  std::string outputFile = "Lennard-Jones.dat";
  std::ofstream out(outputFile);
  std::string macroFile = "LJmacro.dat";
  std::ofstream macro(macroFile);

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
        out<<box.apply_pbc(make_real3(position[index[id]]))<<endl; //!

      macro<<step*simParams.dt<<" ";
      macro<<getTotalEnergy(integrator, particles)<<" ";
      macro<<getTotalMomentum(particles)<<" ";
      macro<<getThermalEnergy(particles)<<endl;
    }
  }

  sys->finish();

  return 0;
} //!
