# Individual Meeting 3 (31.05)

## Code

- what threshold values should I keep in mind? the programs may have bad default values

- Field interpolation: **unaveraged** or **volume averaged**? (recommended: density-> vol averaged, velocity-> unaveraged)

- should the number of grid cells be the same as the number of particles? (or should they be gid cells = sqrt(number of particles)). For illustris Dark3 should the grid size fed to DTFE be 455 or something like 512/256? (choose 256 or 128)

- NEXUS thresholds:

  - threshold_1: *'--minSize 5.e13' gives the smallest mass (in Msolar/h units) that an object must have to be considered a valid cluster.'* (10^(12.5) group) 10^13- cluster

  - threshold_2: *'--minSize 10' specify the minimum volume in (Mpc/h)^3 that a distinct filament must have to be considered a valid feature'* here the features are filaments

  - threshold_3: -||- same as for threshold_2 but related to walls

  - **look** 2013 paper

- should I also get velocity fields from DTFE?

## Writing Skeleton

- what is the purpose of this research really? how can i structure it in terms of a research question and subquestions?
(what are the original ; what determines whether a mass element ends up in a filament, cluster or a void; identification of proto-filamentary -cluster -void structure in the primordial universe) (parts form a current structure that moved from far away? how far) (coherence of the structure: how irregular)(identify which mass elements are gravitationally bound to a given filament)

- Illustris Simulation

    Introduce the simulation, the data, (any important characteristics that must be mentioned in the report?) and add some of the initial histograms of the data

    Data, particles

- NEXUS

    Talk about the algorithm a bit (mention the important parameters); there are a lot of graphs that are in the pythonNEXUS.py file. Any of note? If not, then make my own to show what Cautun et.al. 2014 showed but to a bigger redshift.

- lagrangian description of evolving mass distribution

- Gravitational Collapse bit

    Here I can put some snapshots from the gifs to highlight visually how objects do tend to go void->wall->filament->node

- DTFE

    Find a way to put intermediary results here (and again mention important characteristics and numbers of the algorithm)
