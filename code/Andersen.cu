# include "uammd.cuh"
# include "utils/InputFile.h"
# include "utils/InitialConditions.cuh"
# include "Interactor/Potential/Potential.cuh"
# include "Interactor/NeighbourList/CellList.cuh"
# include "Interactor/PairForces.cuh"
# include "Integrator/VerletNVE.cuh" //!

using namespace uammd;
using std::make_shared;
using std::endl;

struct InputParameters {
  int numberOfParticles;
  real L;
  real dt;
  std::string outputFile;
  std::string macroFile;
  int numberOfSteps;
  int printEverynSteps;
  real epsilon;
  real sigma;
  real cutOff;
  real particleEnergy;
  real thermalEnergy;
  real meanFreeTime; //!
  real mass;
  int checkpointEverynSteps;
  std::string inputFile; //!
};

InputParameters readParameterFile(std::shared_ptr<System> sys)
{
  if(!std::ifstream("data.main").good()) {
    sys->log<System::WARNING>("File data.main not found. Creating file with default values.");
    std::ofstream defaultParameters("data.main");
    if(not defaultParameters.is_open()) {
      sys->log<System::CRITICAL>("Unable to create data.main file. Halting program.");
      exit(-1);
    }
    defaultParameters<<"numberOfParticles 100000"<<endl;
    defaultParameters<<"boxSize 128"<<endl;
    defaultParameters<<"timeStep 0.01"<<endl;
    defaultParameters<<"outputFile Lennard-Jones.dat"<<endl;
    defaultParameters<<"measurementsFile LJmacro.dat"<<endl;
    defaultParameters<<"numberOfSteps 10000"<<endl;
    defaultParameters<<"printEverynSteps 1000"<<endl;
    defaultParameters<<"epsilon 1.0"<<endl;
    defaultParameters<<"sigma 1.0"<<endl;
    defaultParameters<<"cutOff 2.5"<<endl;
    defaultParameters<<"particleEnergy 1.0"<<endl;
    defaultParameters<<"thermalEnergy 1.0"<<endl;
    defaultParameters<<"meanFreeTime 1.0"<<endl; //!
    defaultParameters<<"mass 1.0"<<endl;
  }
  InputFile parameterFile("data.main", sys);
  InputParameters params;

  parameterFile.getOption("numberOfParticles",
    InputFile::Required)>>params.numberOfParticles;
  parameterFile.getOption("boxSize",
    InputFile::Required)>>params.L;
  parameterFile.getOption("timeStep",
    InputFile::Required)>>params.dt;
  parameterFile.getOption("outputFile",
    InputFile::Required)>>params.outputFile;
  parameterFile.getOption("measurementsFile",
    InputFile::Required)>>params.macroFile;
  parameterFile.getOption("numberOfSteps",
    InputFile::Required)>>params.numberOfSteps;
  parameterFile.getOption("printEverynSteps",
    InputFile::Required)>>params.printEverynSteps;
  parameterFile.getOption("epsilon",
    InputFile::Required)>>params.epsilon;
  parameterFile.getOption("sigma",
    InputFile::Required)>>params.sigma;
  parameterFile.getOption("cutOff",
    InputFile::Required)>>params.cutOff;
  parameterFile.getOption("particleEnergy",
    InputFile::Required)>>params.particleEnergy;
  parameterFile.getOption("thermalEnergy",
    InputFile::Required)>>params.thermalEnergy;
  parameterFile.getOption("meanFreeTime",
    InputFile::Required)>>params.meanFreeTime; //!
  params.checkpointEverynSteps = 0;
  parameterFile.getOption("checkpointEverynSteps",
    InputFile::Optional)>>params.checkpointEverynSteps;
  parameterFile.getOption("inputFile",
    InputFile::Optional)>>params.inputFile; //!
  parameterFile.getOption("mass",
    InputFile::Required)>>params.mass;

  return params;
}

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

__global__ void thermalise(int N,
                           real3 * velocity,
                           real thermalVelocity,
                           real probability){
  const int id = blockIdx.x*blockDim.x + threadIdx.x;
  if(id > N) return;

  Saru rng(id);

  if(rng.f() < probability) {
    velocity[id] = make_real3(rng.gf(0, thermalVelocity),
                              rng.gf(0, thermalVelocity).x);
  }

  return;
}

class Andersen: public VerletNVE {
  public:
    struct Parameters: VerletNVE::Parameters {
      real thermalEnergy;
      real meanFreeTime;
      real mass;
      std::shared_ptr<ParticleData> particles;
    } AndersenParameters;

  Andersen(std::shared_ptr<ParticleData> particles,
           std::shared_ptr<System> sys,
           Parameters params) :
           VerletNVE(particles, sys, params){
    AndersenParameters.particles = particles;
    AndersenParameters.thermalEnergy = params.thermalEnergy;
    AndersenParameters.meanFreeTime = params.meanFreeTime;
  }

  virtual void forwardTime() override{
    VerletNVE::forwardTime();

    real collisionProbability
      = AndersenParameters.dt/AndersenParameters.meanFreeTime;

    real thermalVelocity = sqrtf(AndersenParameters.thermalEnergy
                                 /AndersenParameters.mass);

    auto velocity
      = AndersenParameters.particles->getVel(access::location::gpu,
                                             access::mode::readwrite);

    int numberOfParticles = AndersenParameters.particles->getNumParticles();
    int Nthreads = 128;
    int Nblocks = numberOfParticles/Nthreads
                  + ((numberOfParticles%Nthreads)?1:0);

    thermalise<<<Nblocks,Nthreads,0,0>>>(numberOfParticles,
                                         velocity.raw(),
                                         thermalVelocity,
                                         collisionProbability);
  }
};

int main(int argc, char *argv[]){

  auto sys = make_shared<System>(argc, argv);

  sys->log<System::MESSAGE>("Reading parameters from data.main."); //!
  InputParameters simParams = readParameterFile(sys); //!

  int numberOfParticles = simParams.numberOfParticles;
  auto particles
    = make_shared<ParticleData>(numberOfParticles, sys);

  real L = simParams.L;

  Box box(make_real3(L, L, L));
  bool periodicityX = true, periodicityY = true,
       periodicityZ = true;
  box.setPeriodicity(periodicityX, periodicityY,
                     periodicityZ); //!

  if(simParams.inputFile.empty()) {
    sys->log<System::MESSAGE>("Creating new initial positions.");

    auto position
      = particles->getPos(access::location::cpu,
                          access::mode::write);

    auto initial =  initLattice(box.boxSize,
                                numberOfParticles, sc);

    std::copy(initial.begin(), initial.end(), position.begin());
  } else {
    sys->log<System::MESSAGE>("Reading initial positions and velocities from file.");
    auto position
      = particles->getPos(access::location::cpu,
                          access::mode::write);
    auto velocity
      = particles->getVel(access::location::cpu,
                          access::mode::write);

    std::string inputFile = simParams.inputFile;
    std::ifstream in(inputFile);

    for(int i = 0; i < numberOfParticles; ++i) {
      in>>position[i].x>>position[i].y>>position[i].z
        >>velocity[i].x>>velocity[i].y>>velocity[i].z;
      position[i].w = 0;
    }
  } //!

  {
    auto mass
      = particles->getMass(access::location::cpu,
                           access::mode::write);
    std::fill(mass.begin(), mass.end(), simParams.mass);
  }

  using Verlet = Andersen;
  Verlet::Parameters VerletParams;
  VerletParams.dt = simParams.dt;
  VerletParams.initVelocities = true;
  VerletParams.energy = simParams.particleEnergy;
  VerletParams.thermalEnergy = simParams.thermalEnergy;
  VerletParams.mass = simParams.mass;
  VerletParams.meanFreeTime = simParams.meanFreeTime; //!

  if(simParams.inputFile.empty()) {
    sys->log<System::MESSAGE>("UAMMD will generate new velocities.");
    VerletParams.initVelocities = true;
  } else {
    VerletParams.initVelocities = false;
  }

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);//!

  auto LJPotential = make_shared<Potential::LJ>(sys);
  {
    Potential::LJ::InputPairParameters LJParams;
    LJParams.epsilon = simParams.epsilon;
    LJParams.sigma = simParams.sigma;
    LJParams.cutOff = simParams.cutOff;
    LJParams.shift = true;
    LJPotential->setPotParameters(0, 0, LJParams);
  }

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

  std::string outputFile = simParams.outputFile;
  std::ofstream out(outputFile);

  std::string macroFile = simParams.macroFile;
  std::ofstream macro(macroFile);

  int numberOfSteps = simParams.numberOfSteps;
  int printEverynSteps = simParams.printEverynSteps;

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

    if(simParams.checkpointEverynSteps > 0
       and step % simParams.checkpointEverynSteps == 1) {
      auto position
        = particles->getPos(access::location::cpu,
                            access::mode::read);
      auto velocity
        = particles->getVel(access::location::cpu,
                            access::mode::read);


      std::string checkpointFile
        = "checkpoint."
          + std::to_string(step/simParams.checkpointEverynSteps)
          + ".dat";
      std::ofstream checkpoint(checkpointFile);

      for(int i = 0; i < numberOfParticles; ++i) {
        checkpoint<<position[i].x<<" "<<position[i].y<<" "
                  <<position[i].z<<" "<<velocity[i]<<endl;
      }
    } //!
  }

  sys->finish();

  return 0;
}
