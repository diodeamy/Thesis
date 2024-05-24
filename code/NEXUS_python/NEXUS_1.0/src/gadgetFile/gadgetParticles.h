#ifndef GADGET_PARTICLE_HEADER
#define GADGET_PARTICLE_HEADER
#include <vector>


//! data structure to store only particle positions
struct Particle_p
{
    float pos[3];	//particle position
};

//! data structure to store particle positions and masses
struct Particle_pm
{
    float pos[3];	//particle position
    float mass;         //particle mass
};

//! data structure to store both particle positions, velocities and masses
struct Particle_pvm
{
    float pos[3];	//particle position
    float vel[3];	//particle velocity
    float mass;         //particle mass
};

//! data structure to keep the full details of each particle
struct Particle	
{
    float pos[3];	//particle position
    float vel[3];	//particle velocity
    float mass;         //particle mass
    int id;		//particle id
};



/* Rescale the particle positions to the [0,1] interval. */
template< typename particle> inline void rescalePosition(particle *p,
                                                         int const noParticles,
                                                         float const boxLength)
{
    for (int i=0; i<noParticles; ++i)
        for (int j=0; j<3; ++j)
            p[i].pos[j] /= boxLength;
}

/* Return true if the particle is within the box boundaries. */
template< typename particle> inline bool isParticleInBox(particle &p,
                                                         float const box[])
{
    if ( p.pos[0]>=box[0] and p.pos[0]<=box[1] and p.pos[1]>=box[2] and p.pos[1]<=box[3] and p.pos[2]>=box[4] and p.pos[2]<=box[5] )
        return true;
    return false;
}



#endif
