# include "uammd.cuh"
# include "utils/InputFile.h"
# include "utils/InitialConditions.cuh"
# include "Interactor/Potential/Potential.cuh"
# include "Interactor/NeighbourList/CellList.cuh"
# include "Interactor/PairForces.cuh"
# include "Integrator/BrownianDynamics.cuh" //!

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
  real thermalEnergy;
  real viscosity;
  real hydrodynamicRadius; //!
  std::string inputFile;
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
    defaultParameters<<"timeStep 0.001"<<endl;
    defaultParameters<<"outputFile Lennard-Jones.dat"<<endl;
    defaultParameters<<"measurementsFile LJmacro.dat"<<endl;
    defaultParameters<<"numberOfSteps 10000"<<endl;
    defaultParameters<<"printEverynSteps 1000"<<endl;
    defaultParameters<<"epsilon 1.0"<<endl;
    defaultParameters<<"sigma 1.0"<<endl;
    defaultParameters<<"cutOff 2.5"<<endl;
    defaultParameters<<"thermalEnergy 1.0"<<endl;
    defaultParameters<<"viscosity 1.0"<<endl;
    defaultParameters<<"hydrodynamicRadius 1.0"<<endl; //!
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
  parameterFile.getOption("thermalEnergy",
    InputFile::Required)>>params.thermalEnergy;
  parameterFile.getOption("viscosity",
    InputFile::Required)>>params.viscosity;
  parameterFile.getOption("hydrodynamicRadius",
    InputFile::Required)>>params.hydrodynamicRadius; //!
  parameterFile.getOption("inputFile",
    InputFile::Optional)>>params.inputFile;

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

    std::string inputFile = simParams.inputFile;
    std::ifstream in(inputFile);

    for(int i = 0; i < numberOfParticles; ++i) {
      in>>position[i].x>>position[i].y>>position[i].z;
      position[i].w = 0;
    }
  }

  using EulerMaruyama = BD::EulerMaruyama;
  EulerMaruyama::Parameters EMParams;
  EMParams.dt = simParams.dt;
  EMParams.temperature = simParams.thermalEnergy;
  EMParams.viscosity = simParams.viscosity;
  EMParams.hydrodynamicRadius = simParams.hydrodynamicRadius;

  auto integrator
    = make_shared<EulerMaruyama>(particles, sys, EMParams);//!

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
        out<<box.apply_pbc(make_real3(position[index[id]]))<<endl;

      macro<<step*simParams.dt<<" ";
      macro<<getTotalEnergy(integrator, particles)<<" ";
      macro<<" N/A N/A N/A "; /* Undefined total momentum */
      macro<<simParams.thermalEnergy<<endl; //!
    }
  }

  sys->finish();

  return 0;
}
