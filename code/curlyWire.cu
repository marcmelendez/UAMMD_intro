# include "uammd.cuh"
# include "utils/InitialConditions.cuh"
# include "Interactor/BondedForces.cuh"
# include "Interactor/AngularBondedForces.cuh"
# include "Interactor/TorsionalBondedForces.cuh"
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

    real wireLength = 1.0;
    for(int i = 0; i < numberOfParticles; ++i) {
      position[i].x
        = wireLength/(2*M_PI)*cos(2*M_PI*i/(numberOfParticles - 1));
      position[i].y
        = wireLength/(2*M_PI)*sin(2*M_PI*i/(numberOfParticles - 1));
      position[i].z = 0.0001*i;
      position[i].w = 0;
      velocity[i].x = velocity[i].y = velocity[i].z = 0;
    }
  }

  real L = std::numeric_limits<real>::infinity();
  Box box(make_real3(L, L, L));

  using Verlet = VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = real(0.00001);
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
  }

  {
    std::ofstream angularInfo("data.angularForces");
    if(not angularInfo.is_open()) {
      sys->log<System::CRITICAL>("Unable to create data.angularForces file. Halting program.");
      exit(-1);
    }

    real K = 1.0, theta0 = 4*M_PI/(numberOfParticles - 1);
    angularInfo<<(numberOfParticles - 2)<<endl;
    for(int i = 0; i < numberOfParticles - 2; ++i) {
      angularInfo<<i<<" "<<(i + 1)<<" "<<(i + 2)<<" "<<K<<" "<<theta0<<endl;
    }
  }

  {
    std::ofstream torsionalInfo("data.torsionalForces");
    if(not torsionalInfo.is_open()) {
      sys->log<System::CRITICAL>("Unable to create data.torsionalForces file. Halting program.");
      exit(-1);
    }

    real K = 0.1, phi0 = 0.01;
    torsionalInfo<<(numberOfParticles - 3)<<endl;
    for(int i = 0; i < numberOfParticles - 2; ++i) {
      torsionalInfo<<i<<" "<<(i + 1)<<" "<<(i + 2)<<" "<<(i + 3)
                   <<" "<<K<<" "<<phi0<<endl;
    }
  }

  {
    using HarmonicBonds = BondedForces<BondedType::Harmonic>;
    HarmonicBonds::Parameters bondParameters;
    bondParameters.file = "data.bonds";
    auto bonds = make_shared<HarmonicBonds>(particles, sys, bondParameters);

    integrator->addInteractor(bonds);
  }

  {
    using angularPotentials
      = AngularBondedForces<AngularBondedForces_ns::AngularBond>;
    angularPotentials::Parameters angularParameters;
    angularParameters.readFile = "data.angularForces";
    auto angularForces
      = make_shared<angularPotentials>(particles, sys,
                                       angularParameters,
                                       std::make_shared<AngularBondedForces_ns::AngularBond>(box));
    integrator->addInteractor(angularForces);
  }

  {
    using torsionalPotentials
     = TorsionalBondedForces<TorsionalBondedForces_ns::TorsionalBond>;
    torsionalPotentials::Parameters torsionalParameters;
    torsionalParameters.readFile = "data.torsionalForces";
    auto torsionalForces
     = make_shared<torsionalPotentials>(particles, sys,
                                        torsionalParameters,
                                        std::make_shared<TorsionalBondedForces_ns::TorsionalBond>(box));
    integrator->addInteractor(torsionalForces);
  }

  std::string outputFile = "curlyWire.dat";
  std::ofstream out(outputFile);

  int numberOfSteps = 600000;
  int printEverynSteps = 10000;

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
