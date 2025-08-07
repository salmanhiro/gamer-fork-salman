#include "GAMER.h"
#include <cstdlib>
#include <cmath>

#ifdef PARTICLE
#ifdef SUPPORT_GSL
   #include <gsl/gsl_rng.h>
   #include <gsl/gsl_randist.h>
#endif

extern bool   ParStream_Use_Massive;
extern double  ParStream_SigmaX;
extern double  ParStream_SigmaY;
extern double  ParStream_SigmaZ;
extern double  ParStream_BulkSigmaX;
extern double  ParStream_BulkSigmaY;
extern double  ParStream_BulkSigmaZ;
extern double  ParStream_Width;
extern double ParStream_Mass;

// Simple Gaussian random number generator using Box-Muller transform
double rand_normal_parpos(double mean, double stddev) {
    const double u1 = rand() / (RAND_MAX + 1.0);
    const double u2 = rand() / (RAND_MAX + 1.0);
    const double z0 = sqrt(-2.0 * log(u1)) * cos(2 * M_PI * u2);
    return z0 * stddev + mean;
}


//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_ParticleStellarStream
// Description :  Initialize all particle attributes for the particle test problem
//                --> Modified from "Par_Init_ByFile.cpp"
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_FindHomePatch_UniformGrid()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//                4. File format: plain C binary in the format [Number of particles][Particle attributes]
//                   --> [Particle 0][Attribute 0], [Particle 0][Attribute 1], ...
//                   --> Note that it's different from the internal data format in the particle repository,
//                       which is [Particle attributes][Number of particles]
//                   --> Currently it only loads particle mass, position x/y/z, and velocity x/y/z
//                       (and exactly in this order)
//
// Parameter   :  NPar_ThisRank   : Number of particles to be set by this MPI rank
//                NPar_AllRank    : Total Number of particles in all MPI ranks
//                ParMass         : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z     : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z     : Particle velocity array with the size of NPar_ThisRank
//                ParTime         : Particle time     array with the size of NPar_ThisRank
//                ParType         : Particle type     array with the size of NPar_ThisRank
//                AllAttributeFlt : Pointer array for all particle floating-point attributes
//                                  --> Dimension = [PAR_NATT_FLT_TOTAL][NPar_ThisRank]
//                                  --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                      to access the data
//                AllAttributeInt : Pointer array for all particle integer attributes
//                                  --> Dimension = [PAR_NATT_INT_TOTAL][NPar_ThisRank]
//                                  --> Use the attribute indices defined in Field.h to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttributeFlt, AllAttributeInt
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_ParticleStellarStream( const long NPar_ThisRank, const long NPar_AllRank,
                                       real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                       real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                       long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                       long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   srand(314); // or use time(NULL) for non-deterministic


// define the particle attribute arrays
   real_par *ParFltData_AllRank[PAR_NATT_FLT_TOTAL];
   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   ParFltData_AllRank[v] = NULL;
   long_par *ParIntData_AllRank[PAR_NATT_INT_TOTAL];
   for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   ParIntData_AllRank[v] = NULL;

// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 ) {

//    allocate memory for particle attribute arrays
      ParFltData_AllRank[PAR_MASS] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSX] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSY] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSZ] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELX] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELY] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELZ] = new real_par [NPar_AllRank];

      ParIntData_AllRank[PAR_TYPE] = new long_par [NPar_AllRank];

      long p = 0;

      if ( ParStream_Use_Massive ) {

         const double Stream_Length    = amr->BoxSize[0];  // kpc along X

         // stream is assumed in x diretion
         const double x0 = 0.5 * amr->BoxSize[0] - 0.5 * Stream_Length;
         const double y0 = 0.5 * amr->BoxSize[1];
         const double z0 = 0.5 * amr->BoxSize[2];
         
         #ifdef SUPPORT_GSL
         const gsl_rng_type *T;
         gsl_rng *rng;
         gsl_rng_env_setup();
         T   = gsl_rng_default;
         rng = gsl_rng_alloc(T);
         gsl_rng_set(rng, 314); // fixed seed
         #endif

      for (long p = 0; p < NPar_AllRank; p++) {

            double x, y, z;

            do {
               x = x0 + Stream_Length * (double)p / NPar_AllRank;  // uniform distribution along X
               x += rand_normal_parpos(0.0, 0.5);  // small noise around the streamline

            } while ( x < 0.0 || x >= amr->BoxSize[0] );

            do {
               y = rand_normal_parpos(y0, 0.5 * ParStream_Width);
            } while ( y < 0.0 || y >= amr->BoxSize[1] );

            do {
               z = rand_normal_parpos(z0, 0.5 * ParStream_Width);
            } while ( z < 0.0 || z >= amr->BoxSize[2] );


         // Velocities
         const double vx = rand_normal_parpos(ParStream_BulkSigmaX, ParStream_SigmaX);
         const double vy = rand_normal_parpos(0.0, ParStream_SigmaY);
         const double vz = rand_normal_parpos(0.0, ParStream_SigmaZ);

         // stream particles are assumed massless
         ParFltData_AllRank[PAR_MASS][p] = ParStream_Mass;
         ParFltData_AllRank[PAR_POSX][p] = real_par(x);
         ParFltData_AllRank[PAR_POSY][p] = real_par(y);
         ParFltData_AllRank[PAR_POSZ][p] = real_par(z);

         ParFltData_AllRank[PAR_VELX][p] = real_par(vx);
         ParFltData_AllRank[PAR_VELY][p] = real_par(vy);
         ParFltData_AllRank[PAR_VELZ][p] = real_par(vz);

         ParIntData_AllRank[PAR_TYPE][p] = PTYPE_GENERIC_MASSIVE;
      }

         #ifdef SUPPORT_GSL
         gsl_rng_free(rng);
         #endif
      }
   } // if ( MPI_Rank == 0 )

// send particle attributes from the master rank to all ranks
   Par_ScatterParticleData( NPar_ThisRank, NPar_AllRank, _PAR_MASS|_PAR_POS|_PAR_VEL, _PAR_TYPE,
                            ParFltData_AllRank, ParIntData_AllRank, AllAttributeFlt, AllAttributeInt );

// synchronize all particles to the physical time on the base level
   for (long p=0; p<NPar_ThisRank; p++)
      ParTime[p] = (real_par)Time[0];

// free resource
   if ( MPI_Rank == 0 )
   {
      for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   delete [] ParFltData_AllRank[v];
      for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   delete [] ParIntData_AllRank[v];
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_Particle



#endif // #ifdef PARTICLE
