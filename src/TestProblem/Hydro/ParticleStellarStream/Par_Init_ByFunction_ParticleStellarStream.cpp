#include "GAMER.h"

#ifdef PARTICLE
#ifdef SUPPORT_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif

extern int    ParStream_NPar[3];
extern double ParStream_Point_Mass;
extern double ParStream_Par_Sep;
extern bool   ParStream_Use_Tracers;
extern bool   ParStream_Use_Massive;




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


   long NPar_All = 0;
   if ( ParStream_Use_Massive ) NPar_All += 2;
   if ( ParStream_Use_Tracers ) NPar_All += ParStream_NPar[0]*ParStream_NPar[1]*ParStream_NPar[2];

   if ( NPar_All != NPar_AllRank )
      Aux_Error( ERROR_INFO, "total number of particles found [%ld] != expect [%ld] !!\n",
                 NPar_All, NPar_AllRank );

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

         const int    Stream_Nx        = ParStream_NPar[0];
         const int    Stream_Nz        = ParStream_NPar[2];
         const int    NStream          = Stream_Nx * Stream_Nz;

         const double Stream_Length    = 100.0;  // kpc along X
         const double Stream_Thickness = 20.0;   // kpc in Z
         const double Stream_Height    = 20.0;   // kpc in Y

         const double dx = Stream_Length / Stream_Nx;

         const double x0 = 0.5 * amr->BoxSize[0] - 0.5 * Stream_Length;
         const double y0 = 0.5 * amr->BoxSize[1];
         const double z0 = 0.5 * amr->BoxSize[2];

         const double BulkVelX = 150.0;     // km/s
         const double VelDisp  = 1.0;       // km/s

         #ifdef SUPPORT_GSL
         const gsl_rng_type *T;
         gsl_rng *rng;
         gsl_rng_env_setup();
         T   = gsl_rng_default;
         rng = gsl_rng_alloc(T);
         gsl_rng_set(rng, 314); // fixed seed
         #endif

         for (int ix = 0; ix < Stream_Nx; ix++)
         for (int iz = 0; iz < Stream_Nz; iz++) {
            const double x_kpc = x0 + (ix + 0.5) * dx;

            #ifdef SUPPORT_GSL
            const double dy = gsl_ran_gaussian(rng, Stream_Height / 2.0);
            const double dz = gsl_ran_gaussian(rng, Stream_Thickness / 2.0);
            const double dVx = gsl_ran_gaussian(rng, VelDisp);
            const double dVy = gsl_ran_gaussian(rng, VelDisp);
            const double dVz = gsl_ran_gaussian(rng, VelDisp);
            #else
            const double dy = 0.0;
            const double dz = 0.0;
            const double dVx = 0.0;
            const double dVy = 0.0;
            const double dVz = 0.0;
            #endif

            ParFltData_AllRank[PAR_MASS][p]  = 0.0;  // massless to avoid gravity
            ParFltData_AllRank[PAR_POSX][p]  = real_par(x_kpc);
            ParFltData_AllRank[PAR_POSY][p]  = real_par(y0 + dy);
            ParFltData_AllRank[PAR_POSZ][p]  = real_par(z0 + dz);

            ParFltData_AllRank[PAR_VELX][p]  = real_par(BulkVelX + dVx);
            ParFltData_AllRank[PAR_VELY][p]  = real_par(dVy);
            ParFltData_AllRank[PAR_VELZ][p]  = real_par(dVz);

            ParIntData_AllRank[PAR_TYPE][p] = PTYPE_GENERIC_MASSIVE;

            p++;
         }

         #ifdef SUPPORT_GSL
         gsl_rng_free(rng);
         #endif
      }
 // if ( ParStream_Use_Massive )

      if ( ParStream_Use_Tracers ) {

   const int    Stream_Nx        = ParStream_NPar[0];
   const int    Stream_Nz        = ParStream_NPar[2];
   const int    NStream          = Stream_Nx * Stream_Nz;

   const double Stream_Length    = 100.0;  // kpc along X
   const double Stream_Thickness = 20.0;    // kpc in Z
   const double Stream_Height    = 20.0;    // kpc in Y

   const double dx = Stream_Length / Stream_Nx;

   const double x0 = 0.5 * amr->BoxSize[0] - 0.5 * Stream_Length;
   const double y0 = 0.5 * amr->BoxSize[1];  // fixed central Y
   const double z0 = 0.5 * amr->BoxSize[2];  // centered in Z

   const double BulkVelX = 150.0;     // km/s bulk motion
   const double VelDisp  = 1.0;       // km/s velocity dispersion (cold)

   // Initialize RNG locally
   #ifdef SUPPORT_GSL
   const gsl_rng_type *T;
   gsl_rng *rng;
   gsl_rng_env_setup();
   T   = gsl_rng_default;
   rng = gsl_rng_alloc(T);
   gsl_rng_set(rng, 314); // fixed seed for reproducibility
   #endif

   for (int ix = 0; ix < Stream_Nx; ix++)
   for (int iz = 0; iz < Stream_Nz; iz++) {
      const double x_kpc = x0 + (ix + 0.5) * dx;

      #ifdef SUPPORT_GSL
      const double dy = gsl_ran_gaussian(rng, Stream_Height / 2.0);
      const double dz = gsl_ran_gaussian(rng, Stream_Thickness / 2.0);
      const double dVx = gsl_ran_gaussian(rng, VelDisp);
      const double dVy = gsl_ran_gaussian(rng, VelDisp);
      const double dVz = gsl_ran_gaussian(rng, VelDisp);
      #else
      const double dy = 0.0;
      const double dz = 0.0;
      const double dVx = 0.0;
      const double dVy = 0.0;
      const double dVz = 0.0;
      #endif

      ParFltData_AllRank[PAR_MASS][p]  = 0.0;
      ParFltData_AllRank[PAR_POSX][p]  = real_par(x_kpc);
      ParFltData_AllRank[PAR_POSY][p]  = real_par(y0 + dy);
      ParFltData_AllRank[PAR_POSZ][p]  = real_par(z0 + dz);

      ParFltData_AllRank[PAR_VELX][p]  = real_par(BulkVelX + dVx);
      ParFltData_AllRank[PAR_VELY][p]  = real_par(dVy);
      ParFltData_AllRank[PAR_VELZ][p]  = real_par(dVz);

      ParIntData_AllRank[PAR_TYPE][p] = PTYPE_TRACER;

      p++;
   }

   // Free RNG
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
