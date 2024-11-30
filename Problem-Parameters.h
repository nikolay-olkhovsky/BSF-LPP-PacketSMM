/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP AlfaMove (surface movement method)
Module: Problem-Parameters.h (Problem Parameters)
Prefix: PP
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/

//#include "_Problems05-1.h"
#include "_Problems-Miscellaneous.h"

#define PP_PROBLEM_NAME	"packet_settings"
#define PP_M 11		// Number of equations (number of rows in *.mtx)
#define PP_N 5		// Number of variables (number of cols in *.mtx)
#define PP_MAX_OBJ_VALUE 		9999.0
/*---------------------------------------------------------------------*/

#define PP_METHOD_NAME "AlFaMove-MPI"

//-------------------------- Compilation Modes -----------------------
//#define PP_DEBUG
//#define PP_MATRIX_OUTPUT	// To output Matrix
#define PP_DATABASE_OUTPUT
#define PP_SAVE_RESULT
#define PP_MPI
//#define PP_PATH "Problems/"
//#define PP_PATH "C:/Random-LP-Problems/"

//=========================== Problem Parameters =========================
#define PP_MM (2*(PP_M+PP_N))	// Maximal number of inequalities
//#define PP_KK 131072				// Maximal number of hyperedges that include surface point: 2^17 = 131072
//#define PP_KK 524288				// Maximal number of hyperedges that include surface point: 2^19 = 524288
//#define PP_KK 1048575				// Maximal number of hyperedges that include surface point: 2^20 = 1 048 576
//#define PP_KK 16777215				// Maximal number of hyperedges that include surface point: 2^25 = 16 777 216
#define PP_MAX_ITER_COUNT		10000000000		// Maximal count of iterations
#define PP_DBL_MAX				1E+308			// Highest value
#define PP_RND_MAX				32767			// This is necessary for compatibility with different compilers
//-------------------------- Input/Outpoot Parameters ---------------------------
#define PP_OUTPUT_LIMIT	30	// Number of Elements to output
#define PP_SETW 16
//------------------------- Matrix format ----------------
#define PP_INFINITY			1E+308		// Highest bound in *_hi.mtx
#define PP_MTX_PREFIX		"lp_"
#define PP_MTX_POSTFIX_A	".mtx"
#define PP_MTX_POSTFIX_B	"_b.mtx"
#define PP_MTX_POSTFIX_C	"_c.mtx"
#define PP_MTX_POSTFIX_HI	"_hi.mtx"
#define PP_MTX_POSTFIX_LO	"_lo.mtx"
#define PP_MTX_POSTFIX_SO	"_so.mtx"	// Solution point
#define PP_MTX_POSTFIX_SP		"_sp.mtx"	// Surface point
#define PP_MTX_POSTFIX_U0	"_u0.mtx"	// Starting point
//-------------------------- Jobs  -----------------------
#define PP_JOB_GET_BEST_DIRECTION	0 
//------------- Vector Projection Onto Halfspace exit codes -------------
#define PP_EXITCODE_DEGENERATE_INEQUALITY		0
#define PP_EXITCODE_ON_HYPERPLANE				1
#define PP_EXITCODE_INSIDE_HALFSPACE			2
#define PP_EXITCODE_PARALLEL					4
#define PP_EXITCODE_RECESSIVE					5
#define PP_EXITCODE_NONDEGENERATE_PROJECTING	9

//=========================== Problem Parameters =========================
#define PP_MAX_RND_INEQUALITIES	40
#define PP_MAX_N				20
#define PP_MAX_MTX_N			81
#define PP_MAX_M				81
#define PP_MAX_MTX_M			61

#define PP_FILE_INI "config.ini"
static std::string PP_PATH;
//static std::string PP_PROBLEM_NAME;
//static int PP_N;												// Space dimension
