#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double ParStream_Dens_Bg;        // background mass density
static double ParStream_Pres_Bg;        // background pressure
static double ParStream_Ang_Freq;       // gas angular frequency
       int    ParStream_NStar;           // number of star particle
       double ParStream_Point_Mass;     // the mass of the active particles
       bool   ParStream_Use_Tracers;    // whether or not to include tracers
       bool   ParStream_Use_Massive;    // whether or not to include massive particles
       double ParStream_SigmaX;          // 1D velocity dispersion of X component in km/s
       double ParStream_SigmaY;          // 1D velocity dispersion of Y component in km/s
       double ParStream_SigmaZ;          // 1D velocity dispersion of Z component in km/s
       double ParStream_Width;          // width of stellar stream in kpc
       double ParStream_BulkSigmaX;        // bulk velocity dispersion of X component in km/s
       double ParStream_BulkSigmaY;        // bulk velocity dispersion of Y component in km/s
       double ParStream_BulkSigmaZ;        // bulk velocity dispersion of Z component in km/s
       double ParStream_Mass;             // stream particle mass (default 0)
// =======================================================================================

// problem-specific function prototypes
#ifdef PARTICLE
void Par_Init_ByFunction_ParticleStellarStream( const long NPar_ThisRank, const long NPar_AllRank,
                                       real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                       real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                       long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                       long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] );
#endif

bool Flag_ParticleStellarStream( const int i, const int j, const int k, const int lv,
                        const int PID, const double *Threshold );




//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

   if ( !OPT__FREEZE_FLUID )
      Aux_Error( ERROR_INFO, "OPT__FREEZE_FLUID must be enabled !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  LoadInputTestProb
// Description :  Read problem-specific runtime parameters from Input__TestProb and store them in HDF5 snapshots (Data_*)
//
// Note        :  1. Invoked by SetParameter() to read parameters
//                2. Invoked by Output_DumpData_Total_HDF5() using the function pointer Output_HDF5_InputTest_Ptr to store parameters
//                3. If there is no problem-specific runtime parameter to load, add at least one parameter
//                   to prevent an empty structure in HDF5_Output_t
//                   --> Example:
//                       LOAD_PARA( load_mode, "TestProb_ID", &TESTPROB_ID, TESTPROB_ID, TESTPROB_ID, TESTPROB_ID );
//
// Parameter   :  load_mode      : Mode for loading parameters
//                                 --> LOAD_READPARA    : Read parameters from Input__TestProb
//                                     LOAD_HDF5_OUTPUT : Store parameters in HDF5 snapshots
//                ReadPara       : Data structure for reading parameters (used with LOAD_READPARA)
//                HDF5_InputTest : Data structure for storing parameters in HDF5 snapshots (used with LOAD_HDF5_OUTPUT)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LoadInputTestProb( const LoadParaMode_t load_mode, ReadPara_t *ReadPara, HDF5_Output_t *HDF5_InputTest )
{

#  ifndef SUPPORT_HDF5
   if ( load_mode == LOAD_HDF5_OUTPUT )   Aux_Error( ERROR_INFO, "please turn on SUPPORT_HDF5 in the Makefile for load_mode == LOAD_HDF5_OUTPUT !!\n" );
#  endif

   if ( load_mode == LOAD_READPARA     &&  ReadPara       == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_READPARA and ReadPara == NULL !!\n" );
   if ( load_mode == LOAD_HDF5_OUTPUT  &&  HDF5_InputTest == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_HDF5_OUTPUT and HDF5_InputTest == NULL !!\n" );

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// --> LOAD_PARA() is defined in "include/TestProb.h"
// **************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",     &VARIABLE,               DEFAULT,      MIN,              MAX               );
// **************************************************************************************************************************
   LOAD_PARA( load_mode, "ParStream_Dens_Bg",     &ParStream_Dens_Bg,        1.0e-2,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "ParStream_Pres_Bg",     &ParStream_Pres_Bg,        1.0e-2,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "ParStream_Ang_Freq",    &ParStream_Ang_Freq,       0.00051668,   Eps_double,       1.0e-3            );
   LOAD_PARA( load_mode, "ParStream_NStar",       &ParStream_NStar,          1000,         1,                100000           );   
   LOAD_PARA( load_mode, "ParStream_Point_Mass",  &ParStream_Point_Mass,     1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "ParStream_Use_Tracers", &ParStream_Use_Tracers,    true,         Useless_bool,     Useless_bool      );
   LOAD_PARA( load_mode, "ParStream_Use_Massive", &ParStream_Use_Massive,    true,         Useless_bool,     Useless_bool      );
   LOAD_PARA( load_mode, "ParStream_SigmaX",      &ParStream_SigmaX,         1.0,          0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "ParStream_SigmaY",      &ParStream_SigmaY,         1.0,          0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "ParStream_SigmaZ",      &ParStream_SigmaZ,         1.0,          0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "ParStream_Width",       &ParStream_Width,           1.0e-2,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "ParStream_BulkSigmaX",  &ParStream_BulkSigmaX,    120.0,        0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "ParStream_BulkSigmaY",  &ParStream_BulkSigmaY,    120.0,        0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "ParStream_BulkSigmaZ",  &ParStream_BulkSigmaZ,    120.0,        0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "ParStream_Mass",        &ParStream_Mass,           0.0,          0.0,              NoMax_double      );


} // FUNCITON : LoadInputTestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
// (1-1) read parameters from Input__TestProb
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

   LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );

   ReadPara->Read( FileName );

   delete ReadPara;


// (2) check the runtime parameters
   if ( !ParStream_Use_Tracers  &&  !ParStream_Use_Massive )
      Aux_Error( ERROR_INFO,
                 "either ParStream_Use_Tracer, ParStream_Use_Massive, or both must be true !!\n" );

#  ifndef TRACER
   if ( ParStream_Use_Tracers )    Aux_Error( ERROR_INFO, "must enable TRACER for ParStream_Use_Tracers !!\n" );
#  endif

#  ifndef GRAVITY
   if ( ParStream_Use_Massive )    Aux_Error( ERROR_INFO, "must enable GRAVITY for ParStream_Use_Massive !!\n" );
#  endif


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = 2*M_PI/ParStream_Ang_Freq;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

// overwrite the total number of particles
#  ifdef PARTICLE
   amr->Par->NPar_Active_AllRank = 0;
   if ( ParStream_Use_Massive )    amr->Par->NPar_Active_AllRank += ParStream_NStar;
   PRINT_RESET_PARA( amr->Par->NPar_Active_AllRank, FORMAT_LONG, "(PAR_NPAR in Input__Parameter)" );
#  endif


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID            = %d\n",     TESTPROB_ID         );
      Aux_Message( stdout, "  total star particle number = %d\n",     amr->Par->NPar_Active_AllRank    );
      Aux_Message( stdout, "  background mass density    = %13.7e\n", ParStream_Dens_Bg     );
      Aux_Message( stdout, "  background pressure        = %13.7e\n", ParStream_Pres_Bg     );
      Aux_Message( stdout, "  angular frequency          = %13.7e\n", ParStream_Ang_Freq    );
      Aux_Message( stdout, "  active particle mass       = %13.7e\n", ParStream_Point_Mass  );
      Aux_Message( stdout, "  include tracer particles   = %d\n",     ParStream_Use_Tracers );
      Aux_Message( stdout, "  x velocity dispersion      = %13.7e\n", ParStream_SigmaX );
      Aux_Message( stdout, "  y velocity dispersion      = %13.7e\n", ParStream_SigmaY );
      Aux_Message( stdout, "  z velocity dispersion      = %13.7e\n", ParStream_SigmaZ );
      Aux_Message( stdout, "  x bulk velocity dispersion = %13.7e\n", ParStream_BulkSigmaX );
      Aux_Message( stdout, "  y bulk velocity dispersion = %13.7e\n", ParStream_BulkSigmaY );
      Aux_Message( stdout, "  z bulk velocity dispersion = %13.7e\n", ParStream_BulkSigmaZ );
      Aux_Message( stdout, "  stellar stream width       = %d\n",     ParStream_Width );
      Aux_Message( stdout, "  stellar particle mass      = %d\n",     ParStream_Mass );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{


   const double dr_Bg [2]  = { x - 0.25*amr->BoxSize[0], y - 0.25*amr->BoxSize[1] };
   const double dr_Mom[2]  = { x - 0.50*amr->BoxSize[0], y - 0.50*amr->BoxSize[1] };
   const double Radius_Bg  = sqrt( dr_Bg [0]*dr_Bg [0] + dr_Bg [1]*dr_Bg [1] );
   const double Radius_Mom = sqrt( dr_Mom[0]*dr_Mom[0] + dr_Mom[1]*dr_Mom[1] );

   const double Velocity = ParStream_Ang_Freq*Radius_Mom;

   const double Cos_theta = dr_Mom[0]/Radius_Mom;
   const double Sin_theta = dr_Mom[1]/Radius_Mom;

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

   Dens = ParStream_Dens_Bg * ( 1.0 + 5.0 * Radius_Bg / amr->BoxSize[0] );
   Pres = ParStream_Pres_Bg * ( 1.0 + 5.0 * Radius_Bg / amr->BoxSize[0] );
   MomX = -Dens*Velocity*Sin_theta;
   MomY = Dens*Velocity*Cos_theta;
   MomZ = 0.0;

   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );    // assuming EoS requires no passive scalars

   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );      // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ParticleStellarStream
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ParticleStellarStream()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
   Output_User_Ptr               = NULL;
   Flag_User_Ptr                 = Flag_ParticleStellarStream;
   Mis_GetTimeStep_User_Ptr      = NULL;
   Aux_Record_User_Ptr           = NULL;
   BC_User_Ptr                   = NULL;
   Flu_ResetByUser_Func_Ptr      = NULL;
   End_User_Ptr                  = NULL;
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr       = Par_Init_ByFunction_ParticleStellarStream;
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ParticleStellarStream



bool Flag_ParticleStellarStream( const int i, const int j, const int k, const int lv,
                                  const int PID, const double *Threshold )
{
    const double dh     = amr->dh[lv];
    const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,
                            amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                            amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh };

    const double MidY = 0.5 * amr->BoxSize[1];
    const double MidZ = 0.5 * amr->BoxSize[2];

    const double dy = Pos[1] - MidY;
    const double dz = Pos[2] - MidZ;

    // Allow full extent in X, thin in Y and Z
    bool Flag = (FABS(dy) < 0.25) && (FABS(dz) < 0.25);

    return Flag;
}