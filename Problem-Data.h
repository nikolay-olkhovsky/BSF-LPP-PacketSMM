/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP AlfaMove (surface movement method)
Module: Problem-Data.h (Problem Data)
Prefix: PD
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-Types.h"		// Problem Parameters 
using namespace std;
//========================== Problem variables ====================================
static int PD_m;					// Current number of inequalities
static int PD_n;					// Space dimension
static int PD_mh;					// Number of hyperplanes that include surface point
static int PD_ma;					// Number of hyperplanes used for pseudoprojection
static int PD_K;					// Number of faces of all possible dimensions
static int PD_iterNo;				// Number of iterations
static double PD_objF_initialValue; // Initial value of objective function
static double PD_objF_u;			// Current value of objective function
static double PD_shiftLength;		// Shift length
//========================== Problem structures ====================================
static PT_matrix_T PD_A;			// Matrix of coefficients of inequalities
static PT_column_T PD_b;			// Column of the constant terms of the system Ax <= PD_b
static PT_vector_T PD_c;			// Gradient of Objective Function
static PT_vector_T PD_e_c;			// PD_e_c = PD_c / ||PD_c||
static PT_vector_T PD_u;			// Current surface point
static PT_vector_T PD_hi;			// Higher bound
static PT_vector_T PD_lo;			// Lower bound
static PT_vector_T PD_previous_u;	// Previous surface point
static PT_column_T PD_norm_a;		// Column of norms of matrix rows
static PT_vector_T PD_objVector;	// Used for pseudoprojecting
//static int PD_index_includingHyperplanes[PP_MM];	// Index of hyperplanes that include surface point u
//static int PD_faceCodeList[PP_KK];					// Face codes
//static int PD_index_activeHyperplanes[PP_MM];		// Index of hyperplanes used for pseudoprojection
static int PD_pointHyperplanes[PP_MM];	// Index of hyperplanes that include surface point u
static int PD_faceCodes[PP_KK];			// Face codes
static int PD_faceHyperplanes[PP_MM];	// Index of hyperplanes used for pseudoprojection

static PT_unsigned_T PD_index = 0;			// Index of current LPP in dataset
static PT_unsigned_T PD_packetSize;         // Number of LPPs in dataset
static double PD_time;                      // Overall time of execution

//========================== Files ==============================
#ifdef PP_DATABASE_OUTPUT
static auto storage = sqlite_orm::make_storage("C:/HS/dataset100000_3.sqlite3",
    sqlite_orm::make_table("problems",
        sqlite_orm::make_column("id", &Problem::id, sqlite_orm::primary_key()),
        sqlite_orm::make_column("N", &Problem::N),
        sqlite_orm::make_column("seed", &Problem::seed),
        sqlite_orm::make_column("high", &Problem::high),
        sqlite_orm::make_column("low", &Problem::low),
        sqlite_orm::make_column("c", &Problem::c)
    ),
    sqlite_orm::make_table("inequalities",
        sqlite_orm::make_column("id", &Inequality::id, sqlite_orm::primary_key()),
        sqlite_orm::make_column("coefficients", &Inequality::coefficients),
        sqlite_orm::make_column("b", &Inequality::b),
        sqlite_orm::make_column("problem_id", &Inequality::problem_id),
        sqlite_orm::foreign_key(&Inequality::problem_id).references(&Problem::id)
    ),
    sqlite_orm::make_table("surface_points",
        sqlite_orm::make_column("id", &SurfacePoint::id, sqlite_orm::primary_key()),
        sqlite_orm::make_column("coefficients", &SurfacePoint::coefficients),
        sqlite_orm::make_column("problem_id", &SurfacePoint::problem_id),
        sqlite_orm::foreign_key(&SurfacePoint::problem_id).references(&Problem::id)
    ),
    sqlite_orm::make_table("precedents",
        sqlite_orm::make_column("id", &Precedent::id, sqlite_orm::primary_key().autoincrement()),
        sqlite_orm::make_column("coefficients", &Precedent::coefficients),
        sqlite_orm::make_column("face", &Precedent::face),
        sqlite_orm::make_column("d", &Precedent::d),
        sqlite_orm::make_column("shift", &Precedent::shift),
        sqlite_orm::make_column("face_count", &Precedent::face_count),
        sqlite_orm::make_column("face_numbers", &Precedent::face_numbers),
        sqlite_orm::make_column("problem_id", &Precedent::problem_id),
        sqlite_orm::foreign_key(&Precedent::problem_id).references(&Problem::id)
    ),
    sqlite_orm::make_table("images",
        sqlite_orm::make_column("id", &Image::id, sqlite_orm::primary_key().autoincrement()),
        sqlite_orm::make_column("density", &Image::density),
        sqlite_orm::make_column("rank", &Image::rank),
        sqlite_orm::make_column("answer_vector", &Image::answer_vector),
        sqlite_orm::make_column("cosine_vector", &Image::cosine_vector),
        sqlite_orm::make_column("num_of_points", &Image::num_of_points),
        sqlite_orm::make_column("data", &Image::data),
        sqlite_orm::make_column("precedent_id", &Image::precedent_id),
        sqlite_orm::make_column("field_points", &Image::field_points),
        sqlite_orm::foreign_key(&Image::precedent_id).references(&Precedent::id)
    )
);
static std::vector<unsigned> PD_ids;
static std::vector<Problem> PD_DB_problems;
static Problem PD_problem;
static std::vector<Inequality> PD_inequalities;
static std::vector<Inequality> PD_DB_inequalities;
static std::vector<SurfacePoint> PD_surfacePoints;
static std::vector<SurfacePoint> PD_DB_surfacePoints;
static Precedent PD_precedent;
static std::vector<Precedent> PD_newPrecedents;
//static PT_unsigned_T PD_currentPointId;
static SurfacePoint PD_currentPoint;
//static SurfacePoint PD_newPoint;
//static std::vector<SurfacePoint> PD_newSurfacePoints;

#else
CProblem* PD_currentProblem;
CArray* PD_currentX0;

CMTXX0PacketWriter*		PD_packetWriter;
CMTXReader*				PD_packetReader;
CMTXReaderX0*			PD_readerX0;
FILE*					PD_stream_sp;
FILE*					PD_stream_pr;
static string			PD_filename;
//========================== Input/Output ====================================
static string PD_problemName;
//===================== Matrix format (with equations only) ===============================
// nor - number of rows; noc - number of columns; non - number of non-zero values
static string PD_MTX_File_A; /* File of matrix A (only non-zero elements):
------------ begin of file -------------
nor noc non // nor=m (number of inequalities); noc=n (dimension); non - number of non-zero values
i_1 j_1 A_{{i_1}{j_1}} // i_1 - index of row; j_1 - index of column
...
i_k j_k A_{{i_k}{j_k}}
------------ end of file----------------*/
static string PD_MTX_File_b; /* File of column b:
------------ begin of file -------------
nor noc // nor=m (number of inequalities); noc=1
b_1
...
b_{nor}
------------ end of file----------------*/
static string PD_MTX_File_lo; /* File of variable lower bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
lo_1
...
lo_{nor}
------------ end of file----------------*/
static string PD_MTX_File_hi; /* File of variable higher bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
lo_1
...
lo_{nor}
------------ end of file----------------*/
static string PD_MTX_File_c; /* File of variable higher bounds:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
c_1
...
c_{nor}
------------ end of file----------------*/

static string PD_MTX_File_sp;/* File of surface point in the following format:
------------ begin of file -------------
nor noc // nor=n (dimension); noc=1
x_1
...
x_{nor}
------------ end of file----------------*/
#endif // PP_DATABASE_OUTPUT
