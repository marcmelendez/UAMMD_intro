# Example programs

The comment ``//!`` marks the ends of code snippets in the text.

## Introduction

* `helloWorld.cpp`: The iconic "Hello, World!" example written in C++.
* `rubberBall.cpp`: Simulation of a bouncing rubber ball.

## Our first UAMMD simulation

* `minimal.cu`: A minimal UAMMD program that only prints a Hello message and
  then exits.
* `freeExpansion.cu`: Ideal gas expanding freely in the vertical direction with
  periodic boundary conditions in the other two directions.
* `Lennard-Jones.cu`: Particles in a periodic box interacting through a
  Lennard-Jones potential.
* `superLennard-Jones.cu`: Improved Lennard-Jones system with logging and
  external parameter file.

## Unleash your own potential

* `vibratingString.cu`: A chain of particles linked by harmonic springs and held
   tight at the ends behaving like a vibrating string.
* `swingingRope.cu`: A chain of particles simulating a rope hanging from one of
   its ends, swinging back and forth in a gravitational field.
* `cable.cu`: A chain of particles simulating a flexible cable bending under its
   own weight.
* `curlyWire.cu`: A chain of particles behaving like a curly wire bending due to
   angular and torsional forces.
* `MorseChain.cu`: A chain of particles connected with Morse bond potentials and
   attached to the origin in a uniform gravitational field.
* `Morse.cu`: Particles in a periodic box interacting through a Morse potential.
* `diatomic.cu`: A modification of `Morse.cu` in which bonds can only link
   particles together in pairs. All other interactions are treated with
   Lennard-Jones potentials.

## Measuring

* `rubberBallE.cpp`: A version of the `rubberBall.cpp` from the introduction 
   which also outputs the mechanical energy of the ball.
