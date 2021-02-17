# include "uammd.cuh"
# include "utils/InitialConditions.cuh"
# include "Interactor/Potential/Potential.cuh"
# include "Interactor/NeighbourList/CellList.cuh"
# include "Interactor/PairForces.cuh"
# include "Integrator/VerletNVE.cuh"

using namespace uammd;
using std::make_shared;
using std::endl;

struct diatomic {
  real De, a, r0, rc;
  real sigma, epsilon;
  Box box;
  thrust::device_vector<int> bonds; //!

  real getCutOff() { return rc; }

  diatomic(real i_De, real i_a, real i_r0,
           real i_sigma, real i_epsilon,
           real i_rc, Box i_box):
       De(i_De), a(i_a), r0(i_r0),
       sigma(i_sigma), epsilon(i_epsilon),
       rc(i_rc), box(i_box) {} //!

  void createBonds(std::shared_ptr<ParticleData> particles,
                   std::shared_ptr<System> sys){

    int numberOfParticles = particles->getNumParticles();

    bonds.resize(numberOfParticles, -1);
    int * bondlist = thrust::raw_pointer_cast(bonds.data()); //!
  	{
      auto position
        = particles->getPos(access::location::cpu,
                            access::mode::read);
      auto velocity
        = particles->getVel(access::location::cpu,
                            access::mode::readwrite);
      auto mass
        = particles->getMass(access::location::cpu,
                             access::mode::read);
      auto id
        = particles->getId(access::location::cpu,
                           access::mode::read);
      auto index
        = particles->getIdOrderedIndices(access::location::cpu); //!
      for(int i = 0; i < numberOfParticles; ++i) {
        int j = bonds[id[i]];
        if(j >= 0) {
          const real3 rij
            = box.apply_pbc(make_real3(position[index[j]])
                            - make_real3(position[i]));
          const real r2 = dot(rij, rij);
          if(r2 > rc*rc) {
            bonds[j] = bonds[id[i]] = -1;
          }
        }
      } //!
      for(int i = 0; i < numberOfParticles - 1; ++i) {
        if(bonds[id[i]] >= 0) continue;
        const real3 ri = make_real3(position[i]);
        for(int j = i + 1; j < numberOfParticles; ++j) {
          if(bonds[id[j]] >= 0) continue;
          const real3 rij
            = box.apply_pbc(make_real3(position[j]) - ri);
          const real r2 = dot(rij, rij); //!
          if(r2 <= rc*rc) {
            bonds[id[j]] = id[i];
            bonds[id[i]] = id[j]; //!
            const real r = sqrtf(r2);
            const real expar = exp(-a*(r - r0));
            const real oneminusexpar = real(1.0) - expar;
            const real invr2 = (sigma*sigma)/r2;
            const real invr6 = invr2*invr2*invr2;
            const real deltaE = 4*epsilon*invr6*(invr6 - real(1.0))
                                - De*(oneminusexpar*oneminusexpar - real(1.0));
            if(deltaE >= 0) {
              velocity[i] += sqrtf(deltaE/mass[i])*rij/r;
              velocity[j] -= sqrtf(deltaE/mass[j])*rij/r;
            } else {
              velocity[i] -= sqrtf(-deltaE/mass[i])*rij/r;
              velocity[j] += sqrtf(-deltaE/mass[j])*rij/r;
            }
            break;
          }
        }
      }
    }
  } //!

  struct ForceEnergy{
    real4 * force;
    real * energy;
    Box box;
    real De, a, r0, rc;
    real sigma, epsilon;
    int * bondlist;
    int * id;

    ForceEnergy(Box i_box, real i_rc,
                real4 * i_force, real * i_energy,
                real i_De, real i_a, real i_r0,
                real i_sigma, real i_epsilon,
                int * i_bondlist, int * i_id):
                box(i_box), rc(i_rc),
                force(i_force), energy(i_energy),
                De(i_De), a(i_a), r0(i_r0),
                sigma(i_sigma), epsilon(i_epsilon),
                bondlist(i_bondlist), id(i_id){} //!
    __device__ int getInfo(int index){
      return id[index];
    } //!
    __device__ real4 compute(real4 ri, real4 rj, int id_i, int id_j){
      const real3 rij = box.apply_pbc(make_real3(rj)-make_real3(ri));
      const real r2 = dot(rij, rij);
      if(r2 > 0 and r2 < rc*rc) {
        if(bondlist[id_i] == id_j) {
          const real r = sqrtf(r2);
          const real expar = exp(-a*(r - r0));
          const real oneminusexpar = real(1.0) - expar;
          return make_real4(-2.0*De*a*((expar - 1)*expar/r)*rij,
                            De*(oneminusexpar*oneminusexpar
                                - real(1.0)));
        } else {
          const real invr2 = (sigma*sigma)/r2;
          const real invr6 = invr2*invr2*invr2;
          return make_real4(epsilon*invr6*invr2
                            *(real(24.0) - real(48.0)*invr6)*rij,
                            4*epsilon*invr6*(invr6 - real(1.0)));
        }
      }
      else return real4();
    } //!

    __device__ void set(int id, real4 total){
      force[id] += make_real4(total.x, total.y, total.z, 0);
      energy[id] += total.w;
    }
  };

  ForceEnergy getForceTransverser(Box box, std::shared_ptr<ParticleData> particles){
    auto force
      = particles->getForce(access::location::gpu,
                            access::mode::readwrite).raw();
    auto energy
      = particles->getEnergy(access::location::gpu,
                             access::mode::readwrite).raw();
    auto id
      = particles->getId(access::location::gpu,
                         access::mode::read);
    int * bondlist = thrust::raw_pointer_cast(bonds.data());
    return ForceEnergy(box, rc, force, energy, De, a, r0,
                       sigma, epsilon, bondlist, id.raw());
  }
}; //!

int main(int argc, char *argv[]){

  auto sys = make_shared<System>(argc, argv);

  int numberOfParticles = 1000;
  auto particles
    = make_shared<ParticleData>(numberOfParticles, sys); //!

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
    auto mass
      = particles->getMass(access::location::cpu,
                           access::mode::write);
    real atomicMass = 1.0;
    std::fill(mass.begin(), mass.end(), atomicMass);
} //!

  using Verlet = VerletNVE;
  Verlet::Parameters VerletParams;
  VerletParams.dt = 0.01;
  VerletParams.initVelocities = true;
  VerletParams.energy = 0.001; //!

  auto integrator
    = make_shared<Verlet>(particles, sys, VerletParams);//!

  real De = 1.0;
  real a = 2.0;
  real r0 = 1.0;
  real sigma = 1.0;
  real epsilon = 1.0;
  real rc = 6.5*r0;

  auto diatomicPotential
   = make_shared<diatomic>(De, a, r0, sigma, epsilon, rc, box);
  {
    using diatomicForces = PairForces<diatomic>;
    diatomicForces::Parameters interactionParams;
    interactionParams.box = box;

    auto interaction
      = make_shared<diatomicForces>(particles, sys,
                                    interactionParams,
                                    diatomicPotential);

    integrator->addInteractor(interaction);
  }

  std::string outputFile = "diatomic.dat";
  std::ofstream out(outputFile);

  int numberOfSteps = 100000;
  int printEverynSteps = 100; //!

  for(int step = 0; step < numberOfSteps; ++step) {
    diatomicPotential->createBonds(particles, sys); //!
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
    }
  }

  sys->finish();

  return 0;
}
