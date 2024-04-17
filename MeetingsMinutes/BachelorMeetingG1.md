# Group Meeting 1 (17.04)

## What I found out

### Important

We are the officially the "COSMIC WEB CLUB"! :D

"PROJECT = PLAY"

"density is a convolution"

We will likely be helped by Master student *Tomas Bruyn* and PHD student Rob Scholtens <3

Throughout the project we will be using the [Illustris](https://www.illustris-project.org/) simulation collection. This uses the [Arepo](https://wwwmpa.mpa-garching.mpg.de/~volker/arepo/) moving-mesh code (Galilean-invariant hydrodynamicl simulations on a moving mesh). We will be working on darkmatter-only simulations of $455^3$ particles because this is enough.

#### Reading Data

We'll be working with Illustris-3-Dark. The column names indicate the following:

- $L_{box}[Mpc]$ indicates the length of the box, so that $L_{box}^3$ is the volume of the box
- $N_{snap}$ is the number of timesteps taken
- $m_{DM}$ is not the mass of a particle in the classical sense at all. *It's just some mass for some area of the universe that is simulated?*
- the snapshots are the actual particle distributions at each point in spacetime
- groupcat (Group catalog) are dark matter halos. Dark matter halos are centers where gas falls in, stars and galaxies form, etc.

#### Transforming particles to density ditributions

- to transform particle distributions to density distributions we can find 3 different methods on github

  - CIC looks only at each bin content so it has a very sharp cutoff and you lose some quality
  - TSC also looks at the bins immediately nearby and assigns a lower weight to the particle density in them, which means there is less sharp cutoff
  - DTFE is a lot more involved but it produces some pretty crisp results. It uses Voronoi and tessellation and Delauney triangulation (graph theory dual of Voronoi tessellations)

-

### Nice to know

- look into Carlos Frank and the big four: they were the first (westerners) to actually simulate the formation of cosmic sstructure and they did so when Rien was a *PHD* student on VAX machines at only $32^3$ particles
- simulation starts at redshift of z = 127 because around that time the first density formations start to behave nonlinearly. Before that the simulation just sortf of analytically but linearly hops in time?

## APs

Look at a redshift of z = 0 and do the following:

- play with the shit
- turn particle distributions to density distributions
- make visualisations (3D, xy slices, etc.)
- (can do) make a redhsift map by looking also at the velocities of particles
- (can do) get the group catalog
- (can do) make contour plotss from potentials
