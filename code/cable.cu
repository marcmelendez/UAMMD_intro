# include "uammd.cuh"
# include "utils/InitialConditions.cuh"
# include "Interactor/BondedForces.cuh"
# include "Interactor/AngularBondedForces.cuh"
# include "Interactor/ExternalForces.cuh"
# include "Integrator/VerletNVE.cuh"

using namespace uammd;
using std::make_shared;
using std::endl;

struct gravitationalForce{
  real g;
  gravitationalForce(real numericalValueOfg):g(numericalValueOfg){} //!
  __device__ __forceinline__ real3 force(const real4 &position,
                                         const real &mass){
    return make_real3(0.0f, -mass*g, 0.0f);
  } //!
  __device__ __forceinline__ real energy(const real4 &position,
                                         const real &mass){
    return mass*g*position.y;
  } //!
  std::tuple<const real4 *, const real *> getArrays(ParticleData *particles){
    auto position = particles->getPos(access::location::gpu, access::mode::read);
    auto mass = particles->getMass(access::location::gpu, access::mode::read);
    return std::make_tuple(position.raw(), mass.raw());
  } //!
}; //!

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
    auto mass
      = particles->getMass(access::location::cpu,
                          access::mode::write);

    real ropeLength = 1.0;
    for(int i = 0; i < numberOfParticles; ++i) {
      position[i].x = i*(ropeLength/(numberOfParticles - 1));
      position[i].y = position[i].z = position[i].w = real(0.0);
      velocity[i].x = velocity[i].y = velocity[i].z = real(0.0);
      mass[i] = real(0.001);
    }
  } //!

  real L = std::numeric_limits<real>::infinity();
  Box box(make_real3(L, L, L)); //!

  using Verlet = VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = real(0.00001); //!
  VerletParams.initVelocities=false;

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);//!

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
    bondInfo<<"10"<<endl;

    real cableLength = 1.0;
    for(int i = 0; i < 10; ++i)
      bondInfo<<i<<" "<<i*(cableLength/(numberOfParticles - 1))<<" 0 0 1000.0 0.0"<<endl; //!
  }

  {
    std::ofstream angularInfo("data.angularForces");
    if(not angularInfo.is_open()) {
      sys->log<System::CRITICAL>("Unable to create data.angularForces file. Halting program.");
      exit(-1);
    }

    real K = 2.0, theta0 = 0.0;
    angularInfo<<(numberOfParticles - 2)<<endl;
    for(int i = 0; i < numberOfParticles - 2; ++i) {
      angularInfo<<i<<" "<<(i + 1)<<" "<<(i + 2)<<" "<<K<<" "<<theta0<<endl;
    }
  } //!

  {
    using HarmonicBonds = BondedForces<BondedType::Harmonic>;
    HarmonicBonds::Parameters bondParameters;
    bondParameters.file = "data.bonds";
    auto bonds = make_shared<HarmonicBonds>(particles, sys, bondParameters); //!

    integrator->addInteractor(bonds);
  } //!

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
  } //!

  {
    auto gravity
      = make_shared<ExternalForces<gravitationalForce>>(particles, sys, make_shared<gravitationalForce>(real(9.8)));

    integrator->addInteractor(gravity);
  } //!

  std::string outputFile = "cable.dat";
  std::ofstream out(outputFile);

  int numberOfSteps = 150000;
  int printEverynSteps = 5000; //!

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
