# include "uammd.cuh"
# include "utils/InitialConditions.cuh"
# include "Interactor/BondedForces.cuh"
# include "Integrator/VerletNVE.cuh"

using namespace uammd;
using std::make_shared;
using std::endl;

int main(int argc, char *argv[]){

  auto sys = make_shared<System>(argc, argv);

  int numberOfParticles = 101;
  auto particles
    = make_shared<ParticleData>(numberOfParticles, sys);

  {
    auto position
      = particles->getPos(access::location::cpu,
                          access::mode::write);
    auto velocity
      = particles->getVel(access::location::cpu,
                          access::mode::write);

    real amplitude = 0.1;
    real stringlength = 1.0;
    int modenum = 1;
    for(int i = 0; i < numberOfParticles; ++i) {
      position[i].x = i*(stringlength/(numberOfParticles - 1));
      position[i].y = amplitude*sin(modenum*M_PI*position[i].x);
      position[i].z = position[i].w = 0;
      velocity[i].x = velocity[i].y = velocity[i].z = 0;
    }
  }

  using Verlet = VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = 0.001;
  VerletParams.initVelocities=false;

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);

  {
    std::ofstream bondInfo("data.bonds");
    if(not bondInfo.is_open()) {
      sys->log<System::CRITICAL>("Unable to create data.bonds file. Halting program.");
      exit(-1);
    }

    bondInfo<<(numberOfParticles - 1)<<endl;
    for(int i = 0; i < numberOfParticles - 1; ++i) {
      bondInfo<<i<<" "<<(i + 1)<<" 1000.0 0.01"<<endl;
    }
    bondInfo<<"2"<<endl;
    bondInfo<<"0 0 0 0 1000.0 0.0"<<endl;
    bondInfo<<"100 1 0 0 1000.0 0.0"<<endl;
  }

  {
    using HarmonicBonds = BondedForces<BondedType::Harmonic>;
    HarmonicBonds::Parameters bondParameters;
    bondParameters.file = "data.bonds";
    auto bonds = make_shared<HarmonicBonds>(particles, sys, bondParameters);

    integrator->addInteractor(bonds);
  }

  std::string outputFile = "vibratingString.dat";
  std::ofstream out(outputFile);

  int numberOfSteps = 20000;
  int printEverynSteps = 2000;

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
        out<<position[index[id]]<<endl;
    }
  }

  sys->finish();

  return 0;
}
