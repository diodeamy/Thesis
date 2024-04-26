# Group Meeting 2 (26.04)

## What was done since last time, highlights

- only got to attend the Kapteyn CIT system on Monday this week, so I have not figured out a satisfying workflow with norma and the public network available for computations
- familiarised myself with Illustris API and the data available: a managed to get a rudimentary 3D particle distribution (I graphed my first Universe!!)
- read and understood (to some extent) the two relevant ages from [paper by M. Cautun et al.](https://ui.adsabs.harvard.edu/abs/2014MNRAS.441.2923C/abstract) on Cosmic Web evolution and started to familiarise myself with NEXUS through

## Questions that came up

- what are merger trees and how relevant are they actually to what I'm personally doing
- "for voids, around half of their [...] mass has streamed out into sheets and filaments": what are sheets? are there any other morphological structures that I should know? if so which ones? (ANSWERED: look at the lecture website from Riens webpage!!!)
- do you have actual books where I can learn a little bit more about cosmic structure formation? (like for the master course, are they most relevant?)
- has NEXUS/NEXUS+ been related to the structure found by phase-space-based methods?
- where can I actually find NEXUS?
- for DTFE, is the density at the boundary supposed to be 0 like for an isolated universe or similar to the density inside the box?

## What I found out

### Important

- periodic boundary conditions on the data ( all represented in fourier space on a grid, that's why each box repeats infinitely)
- on a grid waves run from the fundamental wave : smallest wave is 2pi/L: Niquest frequency: you need at least 2 cells to represent a sine
- **you can always shift a box central point because the box repeats infinitely** so you still have data, even on a cluster that is near the edges
- "fingers of god": some clusters have huge velocities which means when you observe them, there is a large random redshift just because of the velocity
- do 128^3 datapoints
- what grid? should be in the order of the # of points because of shot noise (cube root of N)
- use 128^3, 256^3, 512^3 powers of 2

## APs

- do 3 timesteps: get positions, try to make paths of the clumps: do it for 1000 particles: get different colours for the position and plot them together to see movement
- for the love of god make your graphs look nicer
