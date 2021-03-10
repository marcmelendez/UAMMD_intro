# Example programs

The comment ``//!`` marks the ends of code snippets in the text.

## Introduction

* `helloWorld.cpp`: The iconic "Hello, World!" example written in C++
   (*section*: Beginner C++).
* `rubberBall.cpp`: Simulation of a bouncing rubber ball (*section*: simulating
   classical particle physics).

## Our first UAMMD simulation

* `minimal.cu`: A minimal UAMMD program that only prints a Hello message and
  then exits.
* `freeExpansion.cu`: Ideal gas expanding freely in the vertical direction with
  periodic boundary conditions in the other two directions (*sections*:
  1.1&ndash;1.4).
* `Lennard-Jones.cu`: Particles in a periodic box interacting through a
  Lennard-Jones potential (*section*: 1.5).
* `superLennard-Jones.cu`: Improved Lennard-Jones system with logging and
  external parameter file (*sections*: 1.6&ndash;1.7).

## Unleash your own potential

* `vibratingString.cu`: A chain of particles linked by harmonic springs and held
   tight at the ends behaving like a vibrating string (*section*: 2.1).
* `swingingRope.cu`: A chain of particles simulating a rope hanging from one of
   its ends, swinging back and forth in a gravitational field (*section*: 2.2).
* `cable.cu`: A chain of particles simulating a flexible cable bending under its
   own weight (*section*: 2.4).
* `curlyWire.cu`: A chain of particles behaving like a curly wire bending due to
   angular and torsional forces (*section*: 2.4).
* `MorseChain.cu`: A chain of particles connected with Morse bond potentials and
   attached to the origin in a uniform gravitational field (*section*: 2.5).
* `Morse.cu`: Particles in a periodic box interacting through a Morse potential
   (*section*: 2.6).
* `diatomic.cu`: A modification of `Morse.cu` in which bonds can only link
   particles together in pairs. All other interactions are treated with
   Lennard-Jones potentials (*section*: 2.7).

## Measuring

* `rubberBallE.cpp`: A version of the `rubberBall.cpp` from the introduction
   which also outputs the mechanical energy of the ball (*section*: 3.1).
* `measurements.cu`: An improvement of `superLennard-Jones.cu` that outputs the
   total energy and linear momentum, and the thermal energy, which is
   proportional to the absolute kinetic temperature (*sections*:
   3.1&ndash;3.2).
* `checkpoints.cu`: An extension of `measurements.cu` with the ability restore
   simulations from saved checkpoints (*section*: 3.3).
* `Langevin.cu`: Lennard-Jones molecular dynamics simulation with a Langevin
   thermostat (*section*: 3.4).
* `Brownian_dynamics.cu:`: Brownian dynamics simulation with Lennard-Jones
   interactions (*section*: 3.5).
* `shear.cu:`: Brownian dynamics simulation with Lennard-Jones particles
   suspended in a shear flow (*section*: 3.5).
* `Andersen.cu`: Lennard-Jones molecular dynamics simulation with an Andersen
   thermostat (*section*: 3.6).
* `Euler.cu`: Lennard-Jones molecular dynamics integrated with the Euler scheme
  (*section*: 3.6).
* `density.cu`: Density profile of the ideal free gas expansion from sections
   1.1&ndash;1.4 (*section*: 3.7).
