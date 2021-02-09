# include "uammd.cuh"
# include "utils/InitialConditions.cuh"
# include "Interactor/BondedForces.cuh"
# include "Interactor/ExternalForces.cuh"
# include "Integrator/VerletNVE.cuh"

using namespace uammd;
using std::make_shared;
using std::endl;

struct Morse { //!
  struct BondInfo {
    real De, a, r0;
  }; //!
  static __host__ BondInfo readBond(std::istream &in) {
    BondInfo MorseParams;
    in>>MorseParams.De>>MorseParams.a>>MorseParams.r0;
    return MorseParams;
  } //!
  inline __device__ real energy(int i, int j,
                                 const real3 &rij,
                                 const BondInfo &MorseParams) {
    real r = sqrtf(dot(rij, rij));
    real oneminusexpar
      = real(1.0) - exp(-MorseParams.a*(r - MorseParams.r0));
    return MorseParams.De*(oneminusexpar*oneminusexpar - real(1.0));
  } //!
  inline __device__ real3 force(int i, int j,
                                const real3 &rij,
                                const BondInfo &MorseParams) {
    real r = sqrtf(dot(rij, rij));
    real expar = exp(-MorseParams.a*(r - MorseParams.r0));
    return real(2.0)*MorseParams.De*MorseParams.a
                    *((expar - 1)*expar/r)*rij;
  }
}; //!

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
      position[i].x = (1 + i)*(ropeLength/(numberOfParticles - 1)); //!
      position[i].y = position[i].z = position[i].w = real(0.0);
      velocity[i].x = velocity[i].y = velocity[i].z = real(0.0);
      mass[i] = real(0.001);
    }
  }

  using Verlet = VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = real(0.00001);
  VerletParams.initVelocities=false;

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);//!

  {
    std::ofstream bondInfo("data.bonds");
    if(not bondInfo.is_open()) {
      sys->log<System::CRITICAL>("Unable to create data.bonds file. Halting program.");
      exit(-1);
    }

    real De = 0.025, a = 200.0, r0 = 0.01;
    bondInfo<<(numberOfParticles - 1)<<endl;
    for(int i = 0; i < numberOfParticles - 1; ++i) {
      bondInfo<<i<<" "<<(i + 1)<<" "
              <<De<<" "<<a<<" "<<r0<<endl;
    }
    bondInfo<<"1"<<endl;
    bondInfo<<"0 0 0 0 "<<De<<" "<<a<<" "<<r0<<endl; //!
  }

  {
    using MorseBonds = BondedForces<Morse>;
    MorseBonds::Parameters bondParameters;
    bondParameters.file = "data.bonds";
    auto bonds = make_shared<MorseBonds>(particles, sys,
                                         bondParameters);

    integrator->addInteractor(bonds);
  } //!

  {
    auto gravity
      = make_shared<ExternalForces<gravitationalForce>>(particles, sys, make_shared<gravitationalForce>(real(9.8)));

    integrator->addInteractor(gravity);
  } //!

  std::string outputFile = "MorseChain.dat";
  std::ofstream out(outputFile);

  int numberOfSteps = 150000;
  int printEverynSteps = 5000;

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
