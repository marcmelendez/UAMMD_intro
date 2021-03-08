# include <iostream>
# include "uammd.cuh"
# include "Integrator/VerletNVE.cuh"

using namespace uammd;
using std::make_shared;
using std::endl;

void density(std::string filename, real crossSection,
             int numberOfParticles, real4 * position,
             real min, real max, int numberOfBins){ //!

  std::ofstream densityFile;
  densityFile.open(filename, std::ios_base::app);

  real binSize = (max - min)/numberOfBins;
  static std::vector<real> bin;
  bin.resize(numberOfParticles);
  std::fill(bin.begin(), bin.end(), real(0.0)); //!

  for(int i = 0; i < numberOfParticles; ++i) {
    int binNumber = (int) floorf((position[i].z - min)/binSize);

    if(binNumber >= 0 and binNumber < numberOfBins) {
      bin[binNumber]++;
    }
  } //!

  for(int i = 0; i < numberOfBins; ++i) {
    densityFile<<min + (i + real(0.5))*binSize<<" "
               <<bin[i]/(numberOfParticles*binSize*crossSection)
               <<endl;
  }
  densityFile<<endl;
} //!

int main(int argc, char *argv[]){

  auto sys = make_shared<System>(argc, argv);

  int numberOfParticles = 100000;
  auto particles
    = make_shared<ParticleData>(numberOfParticles, sys);//!

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
  }//!

  using Verlet = VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = 0.01;
  VerletParams.initVelocities = true;
  VerletParams.energy = 1.0;//!

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);//!

  std::string outputFile = "density.dat";

  int numberOfSteps = 10000;
  int printEverynSteps = 2000;

  for(int step = 0; step < numberOfSteps; ++step) {
    integrator->forwardTime();

    if(printEverynSteps > 0
       and step % printEverynSteps == 1) {
      auto position
        = particles->getPos(access::location::cpu,
                            access::mode::read);

      density(outputFile, L*L,
              numberOfParticles, position.raw(),
              -real(1.5)*L, real(1.5)*L, 100);
    }
  }

  sys->finish();

  return 0;
} //!
