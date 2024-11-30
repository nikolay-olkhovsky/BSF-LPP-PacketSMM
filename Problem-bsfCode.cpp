/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP AlfaMove (surface movement method)
Module: Problem-bsfCode.cpp (Implementation of the Problem)
Prefix: PC
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-Data.h"			// Problem Types 
#include "Problem-Forwards.h"		// Problem Function Forwards
#include "Problem-bsfParameters.h"	// BSF-skeleton parameters
#include "BSF-SkeletonVariables.h"	// Skeleton Variables
using namespace std;
using namespace SF;
using namespace PF;

void PC_bsf_CopyParameter(PT_bsf_parameter_T parameterIn, PT_bsf_parameter_T* parameterOutP) {
	Vector_Copy(parameterIn.x, parameterOutP->x);
}

void PC_bsf_Start(bool* success) {
	ini::IniFile config;
	config.load(PP_FILE_INI);
	PP_PATH = config["general"]["PP_PATH"].as<string>();
//	PP_PROBLEM_NAME = config["general"]["PP_PROBLEM_NAME"].as<string>();
//	PP_N = config["general"]["PP_N"].as<int>();
	PD_index = 0;
	PD_time = clock();

#ifdef PP_DATABASE_OUTPUT
	try {
		//if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
		storage.sync_schema(true);
		//PD_ids = storage.select(&Problem::id, sqlite_orm::where(sqlite_orm::between(&Problem::id, 130, 140)));
		PD_ids = storage.select(&Problem::id);
		PD_packetSize = PD_ids.size();
		//PD_DB_problems = storage.get_all<Problem>(sqlite_orm::where(sqlite_orm::between(&Problem::id, 130, 140)));
		PD_DB_problems = storage.get_all<Problem>();
		PD_DB_inequalities = storage.get_all<Inequality>();
		PD_DB_surfacePoints = storage.get_all<SurfacePoint>();
		cout << "[" << BSF_sv_mpiRank << "] : " << "DB read success." << endl;
	}
	catch (...) {
		cout << "[" << BSF_sv_mpiRank << "] : " << "DB read failure." << endl;
		try {
			PD_ids = storage.select(&Problem::id);
			//PD_ids = storage.select(&Problem::id, sqlite_orm::where(sqlite_orm::between(&Problem::id, 130, 140)));
			PD_packetSize = PD_ids.size();
			PD_DB_problems = storage.get_all<Problem>();
			//PD_DB_problems = storage.get_all<Problem>(sqlite_orm::where(sqlite_orm::between(&Problem::id, 130, 140)));
			PD_DB_inequalities = storage.get_all<Inequality>();
			PD_DB_surfacePoints = storage.get_all<SurfacePoint>();
			cout << "[" << BSF_sv_mpiRank << "] : " << "DB second read success." << endl;
		}
		catch (...) {
			cout << "[" << BSF_sv_mpiRank << "] : " << "DB second read failure." << endl;
		}
	}
	//storage.begin_transaction();
#else
	PD_packetReader = new CMTXReader(PP_PATH.c_str());
	PD_readerX0 = new CMTXReaderX0(PP_PATH.c_str());
	PD_packetSize = PD_packetReader->packetSize;

	if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
	{
		const char* mtxFile;
		PD_MTX_File_sp = PP_PATH;
		PD_MTX_File_sp += PP_MTX_PREFIX;
		PD_MTX_File_sp += PD_problemName;
		PD_MTX_File_sp += PP_MTX_POSTFIX_SP;
		mtxFile = PD_MTX_File_sp.c_str();
		PD_stream_sp = fopen(mtxFile, "w");
		if (PD_stream_sp == NULL) {
			cout << "Failure of opening file '" << mtxFile << "'.\n";
			*success = false;
			return;
		}

		PD_filename = PP_PATH + "precedents.txt";
		PD_stream_pr = fopen(PD_filename.c_str(), "w");
		if (PD_stream_pr == NULL) {
			cout << "Failure of opening file '" << PD_filename.c_str() << "'.\n";
			*success = false;
			return;
		}
	}
#endif // PP_DATABASE_OUTPUT

	//if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
	//{
	//	PD_packetWriter = new CMTXX0PacketWriter(PP_PATH.c_str(), PD_packetReader->packetSize);
	//	//PD_packetWriter->clearFolder();
	//	PD_packetWriter->open();
	//}
}

void PC_bsf_Init(bool* success) {
	for (int i = 0; i < PP_MM; i++) {
		PD_pointHyperplanes[i] = 0;
		PD_faceHyperplanes[i] = 0;
	}
	//cout << BSF_sv_mpiRank << " : PD_loadNewProblem = " << PD_loadNewProblem << endl;
#ifdef PP_DATABASE_OUTPUT
	vector<double> _buff;
	unsigned problem_id;
	//PD_problem = storage.get<Problem>(PD_ids[PD_index]);
	PD_problem = PD_DB_problems[PD_index];
	problem_id = PD_problem.id;
	//PD_inequalities = storage.get_all<Inequality>(sqlite_orm::where(sqlite_orm::c(&Inequality::problem_id) == PD_problem.id));
	PD_inequalities.clear();
	copy_if(PD_DB_inequalities.begin(), PD_DB_inequalities.end(), back_inserter(PD_inequalities), [&problem_id](const Inequality& item)->bool {return item.problem_id == problem_id; });
	//PD_surfacePoints = storage.get_all<SurfacePoint>(sqlite_orm::where(sqlite_orm::c(&SurfacePoint::problem_id) == PD_problem.id));
	PD_surfacePoints.clear();
	copy_if(PD_DB_surfacePoints.begin(), PD_DB_surfacePoints.end(), back_inserter(PD_surfacePoints), [&problem_id](const SurfacePoint& item)->bool {return item.problem_id == problem_id; });

	PD_currentPoint = PD_surfacePoints.back();

	PD_n = PD_problem.N;
	PD_m = 0;

	for (unsigned i = 0; i < PD_n; i++) {
		for (unsigned j = 0; j < PD_n; j++)
			if (i == j) PD_A[PD_m][j] = 1.;
			else		PD_A[PD_m][j] = 0.;
		PD_b[PD_m] = PD_problem.high;
		PD_m++;
	}
	for (auto& inequality : PD_inequalities) {
		_buff = charToDouble(inequality.coefficients);
		for (unsigned j = 0; j < PD_n; j++)
			PD_A[PD_m][j] = _buff[j];
		PD_b[PD_m] = inequality.b;
		PD_m++;
	}
	for (unsigned i = 0; i < PD_n; i++) {
		for (unsigned j = 0; j < PD_n; j++)
			if (i == j) PD_A[PD_m][j] = -1.;
			else		PD_A[PD_m][j] = 0.;
		PD_b[PD_m] = PD_problem.low;
		PD_m++;
	}

	_buff = charToDouble(PD_problem.c);
	for (unsigned i = 0; i < PD_n; i++) {
		PD_c[i] = _buff[i];
	}

//	PD_currentPointId = PD_currentPoint.id;
	_buff = charToDouble(PD_currentPoint.coefficients);
	for (unsigned i = 0; i < PD_n; i++) {
		PD_u[i] = _buff[i];
	}
#else
	PD_currentProblem = PD_packetReader->readProblem();
	PD_currentX0 = PD_readerX0->readX0();

	PD_m = PD_currentProblem->A->getRows();
	PD_n = PD_currentProblem->A->getCols();
	for (int j = 0; j < PD_n; j++) {
		PD_u[j] = PD_currentX0->getValue(j);
		//PD_x0[j] = 0;
	}
#endif // PP_DATABASE_OUTPUT
	PD_index++;
	MakeColumnOfNorms(PD_A, PD_norm_a);

	if (!PointBelongsPolytope(PD_u, PP_EPS_ZERO)) {
		if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
			cout << "Starting point does not belong to the feasible polytope with precision PP_EPS_ZERO = "
			<< PP_EPS_ZERO << "!!!\n";
		*success = false;
		return;
	}
	PD_iterNo = 0;
	Vector_MakeLike(PD_c, 1, PD_e_c);
//	if(BSF_sv_mpiRank == BSF_sv_mpiMaster)
//		PreparationForIteration(PD_u);
}

//void PC_bsf_MapInit(PT_bsf_parameter_T parameter) {
//	PreparationForIteration(parameter.x);
///**
//#ifdef PP_DEBUG
//	for(int i = 0; i < PD_mh; i++)
//		cout << "[" << BSF_sv_mpiRank << "] : face " << i << " => " << PD_faceCodeList[i] << endl;
//#endif //PP_DEBUG
///**/
//}

void PC_bsf_IterInit(PT_bsf_parameter_T parameter) {
	PreparationForIteration(parameter.x);
}

void PC_bsf_IterOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob) {

	cout << "# " << BSF_sv_iterCounter << "\tTime " << round(elapsedTime);
	cout << "\tx =";
	Print_Vector(parameter.x);
	cout << "\tF(x) = " << setw(PP_SETW) << ObjF(parameter.x) << endl;
}

void PC_bsf_IterOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

void PC_bsf_IterOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

void PC_bsf_IterOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter,
	double elapsedTime, int nextJob)
{
	// Not used
}

void PC_bsf_JobDispatcher(PT_bsf_parameter_T* parameter, int* job, bool* toExit, double t) {
	//cout << "Number of hyperplane subsets = " << PD_K << endl;
}

void PC_bsf_MapF(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T* reduceElem, int* success) {
	int faceCode = *mapElem->faceCode;
	PT_vector_T u;		// current surface point
	PT_vector_T v;		// v = u + PD_objVector (objVector = PP_OBJECTIVE_VECTOR_LENGTH*e_c)
	PT_vector_T w;		// pseudiprojection of v
	double objF_w = -PP_DBL_MAX; // F(w)

	for (int i = 0; i < PP_N; i++)
		reduceElem->face_count[i] = 0;

	if (faceCode == 0) {
		Vector_Zeroing((*reduceElem).d);
		reduceElem->objF_p = -PP_DBL_MAX;
		reduceElem->faceCode = 0;
		return;
	}

	/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
	cout << "------------------------------------ Map(" << PF_MAP_LIST_INDEX << ") ------------------------------------" << endl;
#endif // PP_DEBUG /**/

// Condition for breakpoint: PD_iterNo == 2 && (BSF_sv_addressOffset + BSF_sv_numberInSublist == 2)
	Vector_Copy(BSF_sv_parameter.x, u);
	double objF_u = ObjF(u);
	reduceElem->faceCode = faceCode;

	CodeToSubset(faceCode, PD_faceHyperplanes, &PD_ma);

	/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "Face hyperplanes: {";
	for (int i = 0; i < PD_ma - 1; i++) {
		cout << PD_faceHyperplanes[i] << ", ";
	}
	cout << PD_faceHyperplanes[PD_ma - 1] << "}.\n";
#endif // PP_DEBUG /**/

	int old_PD_ma = PD_ma;
	Vector_MultiplyByNumber(PD_e_c, PP_OBJECTIVE_VECTOR_LENGTH, PD_objVector);

	Vector_Addition(u, PD_objVector, v);

	PseudoprojectionOnFace(v, w, PP_EPS_PROJECTION_ROUND, success);
	if (!*success) {
		cout << "\n\nProcess " << BSF_sv_mpiRank
			<< ". Error in PC_bsf_MapF: Exceeded the maximum number of iterations when calculating pseudoprojection (PP_MAX_PROJECTING_ITER = "
			<< PP_MAX_PROJECTING_ITER << "). It is impossible to calculate Map function for element "
			<< BSF_sv_addressOffset + BSF_sv_numberInSublist << "!\n Perhaps you should decrease parameter PP_EPS_PROJECTION_ROUND.";
		return;
	}

	Vector_Round(w, PP_EPS_PROJECTION_ROUND * 10);

	objF_w = ObjF(w);

	Vector_Subtraction(w, u, (*reduceElem).d);

	double norm_d = Vector_Norm((*reduceElem).d);
	if (norm_d < PP_EPS_ZERO) {
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "\t\t\t\t\t\t\t\t\t\t\t\t\t||d|| = ||w-u|| < PP_EPS_ZERO ===>>> movement is impossible.\n";
#endif // PP_DEBUG /**/
		Vector_Zeroing((*reduceElem).d);
		reduceElem->objF_p = objF_u;
		return;
	}

	Vector_MultiplyEquals((*reduceElem).d, PP_PROBE_LENGTH / norm_d);
	Vector_Round((*reduceElem).d, PP_EPS_PROJECTION_ROUND);

	PT_vector_T p;
	Vector_Addition(u, (*reduceElem).d, p);

	if (!PointBelongsPolytope(p, PP_EPS_POINT_IN_HALFSPACE)) {
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "Shifted point p = ";
		Print_Vector(p);
		cout << "\tnot in feasible polytope ===>>> movement is impossible." << endl;
#endif // PP_DEBUG /**/
		Vector_Zeroing((*reduceElem).d);
		reduceElem->objF_p = objF_u;
		return;
	}

	reduceElem->objF_p = ObjF(p);

	if (RelativeError(objF_u, reduceElem->objF_p) < PP_EPS_ZERO) {
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "u =\t    ";
		Print_Vector(u);
		cout << "\tF(u) = " << setw(PP_SETW) << objF_u << endl;
		cout << "|F(u1)-F(u2)|/|F(u1)| = " << RelativeError(objF_u, reduceElem->objF_p) << " < PP_EPS_ZERO = " << PP_EPS_ZERO << " ===>>> movement is impossible.\n";
#endif // PP_DEBUG /**/

		Vector_Zeroing((*reduceElem).d);
		reduceElem->objF_p = objF_u;
		return;
}

	if (reduceElem->objF_p < objF_u) {
		/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
		cout << "u =\t    ";
		Print_Vector(u);
		cout << "\tF(w) = " << setw(PP_SETW) << reduceElem->objF_p << " < F(u) = " << objF_u << " ===>>> movement is impossible.\n";
#endif // PP_DEBUG /**/
		Vector_Zeroing((*reduceElem).d);
		reduceElem->objF_p = objF_u;
		return;
	}
	reduceElem->face_count[old_PD_ma - 1]++;
	/*DEBUG PC_bsf_MapF**
#ifdef PP_DEBUG
	cout << "d = ";
	Print_Vector((*reduceElem).d);
	cout << endl;
	cout << "Shifted point p = ";
	Print_Vector(p);
	cout << "\tF(p) ="
		<< setw(PP_SETW) << reduceElem->objF_p << "\t\t---> Movement is possible." << endl;
#endif // PP_DEBUG /**/
} // end PC_bsf_MapF

void PC_bsf_MapF_1(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_1* reduceElem, int* success) {
	// Not used
}

void PC_bsf_MapF_2(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_2* reduceElem, int* success) {
	// Not used
}

void PC_bsf_MapF_3(PT_bsf_mapElem_T* mapElem, PT_bsf_reduceElem_T_3* reduceElem, int* success) {
	// Not used
}

void PC_bsf_ParametersOutput(PT_bsf_parameter_T parameter) {
/**
	cout << "=================================================== " << PP_METHOD_NAME << " ====================================================" << endl;
	cout << "Problem name: " << PP_PROBLEM_NAME << endl;

#ifdef PP_MPI
	cout << "Number of Workers: " << BSF_sv_numOfWorkers << endl;
#else
	cout << "No MPI" << endl;
#endif // PP_MPI

#ifdef PP_BSF_OMP
#ifdef PP_BSF_NUM_THREADS
	cout << "Number of Threads: " << PP_BSF_NUM_THREADS << endl;
#else
	cout << "Number of Threads: " << omp_get_num_procs() << endl;
#endif // PP_BSF_NUM_THREADS
#else
	cout << "OpenMP is turned off!" << endl;
#endif // PP_BSF_OMP

#ifdef PP_BSF_FRAGMENTED_MAP_LIST
	cout << "Map List is Fragmented" << endl;
#else
	cout << "Map List is not Fragmented" << endl;
#endif
	cout << "In parameters: m =\t" << PP_M << "\tn = " << PP_N << endl;
	cout << "Actually read: m =\t" << PD_m << "\tn = " << PD_n << endl;
	cout << "PP_EPS_ZERO\t" << PP_EPS_ZERO << endl;
	cout << "PP_EPS_POINT_IN_HALFSPACE\t" << PP_EPS_POINT_IN_HALFSPACE << endl;
	cout << "PP_EPS_MOVING_ON_POLYTOPE\t" << PP_EPS_MOVING_ON_POLYTOPE << endl;
	cout << "PP_EPS_PROJECTION_ROUND\t\t" << PP_EPS_PROJECTION_ROUND << endl;
	cout << "PP_OBJECTIVE_VECTOR_LENGTH\t" << PP_OBJECTIVE_VECTOR_LENGTH << endl;
	cout << "PP_PROBE_LENGTH\t\t\t" << PP_PROBE_LENGTH << endl;
	cout << "PP_START_SHIFT_LENGTH\t\t" << PP_START_SHIFT_LENGTH << endl;
	cout << "--------------- Data ---------------\n";
#ifdef PP_MATRIX_OUTPUT
	cout << "------- Matrix PD_A & Column PD_b -------" << endl;
	Print_Inequalities();
#endif // PP_MATRIX_OUTPUT
	cout << "Vector of norms:\t";
	for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_m); j++) cout << PD_norm_a[j] << ' ';
	if (PP_OUTPUT_LIMIT < PD_m) cout << "	...";
	cout << endl;
	cout << "Obj Function:\t\t";
	Print_Vector(PD_c);
	cout << endl << "\t\t\t";
	Print_Vector(PD_e_c);
	cout << endl;
	cout << "u0 =\t\t";
	Print_Vector(PD_u); cout << "\tF(x) = " << setw(PP_SETW) << ObjF(PD_u) << endl;

#ifdef PP_DEBUG
	if (!PointBelongsPolytope(PD_u, PP_EPS_POINT_IN_HALFSPACE))
		cout << "u0 is outside feasible polytope!!!\n";
	else
		cout << "u0 is belongs feasible polytope.\n";
	cout << "Including hyperplanes:\t"; Print_HyperplanesIncludingPoint(PD_u, PP_EPS_POINT_IN_HALFSPACE);
	cout << "\tPD_mh = " << PD_mh << "\tPD_K = " << PD_K << endl;
	cout << "Including faces:\t"; Print_Number_of_faces(PD_u);
#endif // PP_DEBUG
/**/
}

void PC_bsf_ProblemOutput(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	bool lastPoint = fabs(PD_objF_u - PD_objF_initialValue) < PP_EPS_ZERO;
	if (PD_index % 500 == 0) {
		//storage.insert_range(PD_newSurfacePoints.begin(), PD_newSurfacePoints.end());
		//PD_newSurfacePoints.clear();
		storage.insert_range(PD_newPrecedents.begin(), PD_newPrecedents.end());
		PD_newPrecedents.clear();
		cout << PD_index << " problems processed. Time: " << (clock() - PD_time) / CLOCKS_PER_SEC << endl;
	}
	//cout << setprecision(PP_SETW / 2);

	//PD_objF_u = ObjF(PD_u);

	//cout << "================================================" << endl;
	//cout << "// Elapsed time: " << t << endl;
	//cout << "// Number of iterations: " << PD_iterNo << endl;
	//cout << "// Computed objective value: " << setprecision(16) << PD_objF_u << endl;
	//cout << "// Maximal objective value:  " << PP_MAX_OBJ_VALUE << endl;
	//cout << "// Relative error = " << setprecision(3) << RelativeError(PP_MAX_OBJ_VALUE, PD_objF_u) << setprecision(PP_SETW / 2) << endl;
	//cout << "================================================" << endl;

	//cout << "=============================================" << endl;
	//cout << "Elapsed time: " << t << endl;
	//cout << "Current objective value: " << setprecision(16) << PD_objF_u << endl;
	//cout << "Optimal objective value: " << PP_OPTIMAL_OBJ_VALUE << endl;
	//cout << "Relative error = " << setprecision(PP_SETW / 2) << relativeError(PP_OPTIMAL_OBJ_VALUE, PD_objF_u) << endl;
	//cout << "=============================================" << endl;

	//if (lastPoint) {
	//	cout << "Value of the objective function cannot be refined!\n";
	////	return;
	//}

	//CodeToSubset(reduceResult->subsetCode, PD_index_activeHalfspaces, &PD_ma);
	//cout << "Face dimension: " << PD_n - PD_ma << ".\tGenerating hyperplanes: {";
	//for (int i = 0; i < PD_ma - 1; i++)
	//	cout << PD_index_activeHalfspaces[i] << ", ";
	//cout << PD_index_activeHalfspaces[PD_ma - 1]
	//	<< "}.\tShift = " << PD_shiftLength << "\tF(x) = " << PD_objF_u << endl;

	//cout << "Surface point:\t";
	//for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << PD_u[j];
	//if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	//cout << endl;
	//cout << "Polytope residual: " << PolytopeResidual(PD_u) << endl;

	//cout << "Solution point:\t";
	//Print_Vector(PD_u);	cout << endl;
#ifdef PP_DEBUG
	cout << "Distance to polytope: " << Distance_PointToPolytope(PD_u) << endl;
#endif // PP_DEBUG

}

void PC_bsf_ProblemOutput_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	//ProblemOutput(t);
}

void PC_bsf_ProblemOutput_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Not used
}

void PC_bsf_ProblemOutput_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T parameter, double t) {
	// Not used
}

void PC_bsf_ProcessResults(PT_bsf_reduceElem_T* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* toExit) {
	PT_vector_T u_moved;
	Vector_Copy(PD_u, PD_previous_u);

	if (Vector_Norm(reduceResult->d) < PP_EPS_ZERO) {
		reduceResult->faceCode = 0;
		PD_shiftLength = 0;
		SavePrecedent(reduceResult);
		//if (PD_index < 20) {
		if (PD_index < PD_packetSize) {
			*nextJob = BD_JOB_RESET;
			return;
		}
		else {
			//storage.commit();
			storage.insert_range(PD_newPrecedents.begin(), PD_newPrecedents.end());
			PD_newPrecedents.clear();
			*toExit = true;
			return;
		}
	}

	MovingOnPolytope(PD_u, reduceResult->d, u_moved, PP_EPS_MOVING_ON_POLYTOPE);
	PD_shiftLength = Distance_PointToPoint(PD_u, u_moved);

#ifdef PP_DEBUG
	if (PD_shiftLength < PP_EPS_ZERO) {
		cout << "\nError in MovingOnPolytope(): Shift length < PP_EPS_ZERO! Possible you should increase PP_EPS_MOVING_ON_POLYTOPE or decrease PP_EPS_ZERO.\n";
		exit(1);
	}
#endif // PP_DEBUG

	Vector_Copy(u_moved, PD_u);

	PD_iterNo++;

	PD_objF_u = ObjF(PD_u);

#ifdef PP_DEBUG
	CodeToSubset(reduceResult->faceCode, PD_faceHyperplanes, &PD_ma);
	//
	//
	//
	//
	cout << "_________________________________________________ " << PD_iterNo << " _____________________________________________________" << endl;
	//PD_objF_u = ObjF(PD_u);
	cout << "Face code:" << reduceResult->faceCode << ".\tDimension: " << PD_n - PD_ma << endl;
	cout << "Point:\t";
	Print_Vector(PD_previous_u); cout << "\tF(x) = " << ObjF(PD_previous_u) << endl;
	cout << endl;
	cout << "Face Hyperplanes:\t{";
	for (int i = 0; i < PD_ma - 1; i++)
		cout << PD_faceHyperplanes[i] << ", ";
	cout << PD_faceHyperplanes[PD_ma - 1]
		<< "}.\tShift = " << PD_shiftLength << "\tF(x) = " << PD_objF_u << endl;

	cout << "New surface point:\t";
	Print_Vector(PD_u); cout << "\tF(x) = " << ObjF(PD_u) << endl;
	cout << "Including hyperplanes:\t"; Print_HyperplanesIncludingPoint(PD_u, PP_EPS_POINT_IN_HALFSPACE); cout << endl;
#endif
	//SavePoint();
	SavePrecedent(reduceResult);
	Vector_Copy(PD_u, parameter->x);
}

void PC_bsf_ProcessResults_1(PT_bsf_reduceElem_T_1* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* toExit) {
	// Not used
}

void PC_bsf_ProcessResults_2(PT_bsf_reduceElem_T_2* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* toExit) {
	// Not used
}

void PC_bsf_ProcessResults_3(PT_bsf_reduceElem_T_3* reduceResult, int reduceCounter, PT_bsf_parameter_T* parameter, int* nextJob, bool* toExit) {
	// Not used
}

void PC_bsf_ReduceF(PT_bsf_reduceElem_T* x, PT_bsf_reduceElem_T* y, PT_bsf_reduceElem_T* z) { // z = x + y
	if (x->objF_p > y->objF_p) {
		z->objF_p = x->objF_p;
		for (int j = 0; j < PD_n; j++)
			(*z).d[j] = (*x).d[j];
		z->faceCode = x->faceCode;
	}
	else {
		z->objF_p = y->objF_p;
		for (int j = 0; j < PD_n; j++)
			(*z).d[j] = (*y).d[j];
		z->faceCode = y->faceCode;
	}
	
	for (int i = 0; i < PD_n - 1; i++)
		z->face_count[i] = x->face_count[i] + y->face_count[i];
}

void PC_bsf_ReduceF_1(PT_bsf_reduceElem_T_1* x, PT_bsf_reduceElem_T_1* y, PT_bsf_reduceElem_T_1* z) {
	// Not used
}

void PC_bsf_ReduceF_2(PT_bsf_reduceElem_T_2* x, PT_bsf_reduceElem_T_2* y, PT_bsf_reduceElem_T_2* z) {
	// Not used
}

void PC_bsf_ReduceF_3(PT_bsf_reduceElem_T_3* x, PT_bsf_reduceElem_T_3* y, PT_bsf_reduceElem_T_3* z) {
	// Not used
}

void PC_bsf_SetInitParameter(PT_bsf_parameter_T* parameter) {
	Vector_Copy(PD_u, parameter->x);
}

void PC_bsf_SetListSize(int* listSize) {
	*listSize = PP_KK;
}

void PC_bsf_SetMapListElem(PT_bsf_mapElem_T* elem, int i) {
	elem->faceCode = &(PD_faceCodes[i]);
}

//----------------------- Assigning Values to BSF-skeleton Variables (Do not modify!) -----------------------
void PC_bsfAssignAddressOffset(int value) { BSF_sv_addressOffset = value; }
void PC_bsfAssignIterCounter(int value) { BSF_sv_iterCounter = value; }
void PC_bsfAssignJobCase(int value) { BSF_sv_jobCase = value; }
void PC_bsfAssignMpiMaster(int value) { BSF_sv_mpiMaster = value; }
void PC_bsfAssignMpiRank(int value) { BSF_sv_mpiRank = value; }
void PC_bsfAssignNumberInSublist(int value) { BSF_sv_numberInSublist = value; }
void PC_bsfAssignNumOfWorkers(int value) { BSF_sv_numOfWorkers = value; }
void PC_bsfAssignParameter(PT_bsf_parameter_T parameter) { PC_bsf_CopyParameter(parameter, &BSF_sv_parameter); }
void PC_bsfAssignSublistLength(int value) { BSF_sv_sublistLength = value; }

//---------------------------------- Shared Functions -------------------------
namespace SF {

	//static void SF_Conversion() { // Transformation to inequalities & dimensionality reduction
	//	int m_equation = PD_m;
	//	int m_inequality;
	//	int m_lowerBound;
	//	int m_higherBound;
	//
	//	for (int i = 0; i < m_equation; i++) { // Conversion to inequalities
	//		for (int j = 0; j < PD_n; j++)
	//			PD_A[PD_m][j] = -PD_A[i][j];
	//		PD_b[PD_m] = -PD_b[i];
	//		PD_m++; assert(PD_m <= PP_MM);
	//	}
	//
	//	for (int i = 0; i < PD_m; i++) // Remove negative sign for zero value
	//		for (int j = 0; j < PD_n; j++)
	//			if (PD_A[i][j] == 0)
	//				PD_A[i][j] = 0;
	//	m_inequality = PD_m;
	//
	//	for (int i = 0; i < PD_n; i++) { // Adding lower bound conditions
	//		for (int j = 0; j < PD_n; j++)
	//			PD_A[i + PD_m][j] = 0;
	//		PD_A[i + PD_m][i] = -1;
	//		if (PD_lo[i] == 0)
	//			PD_b[i + PD_m] = 0;
	//		else
	//			PD_b[i + PD_m] = -PD_lo[i];
	//	}
	//	PD_m += PD_n; assert(PD_m <= PP_MM);
	//	m_lowerBound = PD_m;
	//
	//	for (int i = 0; i < PD_n; i++) { // Adding higher bound conditions
	//		if (PD_hi[i] != PP_INFINITY) {
	//			for (int j = 0; j < PD_n; j++)
	//				PD_A[PD_m][j] = 0;
	//			PD_A[PD_m][i] = 1;
	//			PD_b[PD_m] = PD_hi[i];
	//			PD_m++; assert(PD_m <= PP_MM);
	//		}
	//	}
	//	m_higherBound = PD_m;
	//
	//	/**
	//	cout << "-----------------------------------------------------\n";
	//	Print_Inequalities();
	//	cout << "-----------------------------------------------------\n";
	//	cout << "PD_c: "; Print_Vector(PD_c); cout << endl;/**/
	//
	//	RemoveFreeVariables(m_equation, m_inequality, m_lowerBound, m_higherBound);
	//
	//	/**
	//	cout << "-----------------------------------------------------\n";
	//	Print_Inequalities();
	//	cout << "-----------------------------------------------------\n";/**/
	//}

	//inline void SF_ConversionSimple() { // Transformation to inequalities
	//	int m_equation = PD_m;
	//
	//	for (int i = 0; i < m_equation; i++) { // Conversion to inequalities
	//		for (int j = 0; j < PD_n; j++)
	//			PD_A[PD_m][j] = -PD_A[i][j];
	//		PD_b[PD_m] = -PD_b[i];
	//		PD_m++; assert(PD_m <= PP_MM);
	//	}
	//
	//	for (int i = 0; i < PD_m; i++) // Remove negative sign for zero value
	//		for (int j = 0; j < PD_n; j++)
	//			if (PD_A[i][j] == 0)
	//				PD_A[i][j] = 0;
	//
	//	for (int i = 0; i < PD_n; i++) { // Adding lower bound conditions
	//		for (int j = 0; j < PD_n; j++)
	//			PD_A[i + PD_m][j] = 0;
	//		PD_A[i + PD_m][i] = -1;
	//		if (PD_lo[i] == 0)
	//			PD_b[i + PD_m] = 0;
	//		else
	//			PD_b[i + PD_m] = -PD_lo[i];
	//	}
	//	PD_m += PD_n; assert(PD_m <= PP_MM);
	//
	//	for (int i = 0; i < PD_n; i++) { // Adding higher bound conditions
	//		if (PD_hi[i] != PP_INFINITY) {
	//			for (int j = 0; j < PD_n; j++)
	//				PD_A[PD_m][j] = 0;
	//			PD_A[PD_m][i] = 1;
	//			PD_b[PD_m] = PD_hi[i];
	//			PD_m++; assert(PD_m <= PP_MM);
	//		}
	//	}
	//}

	inline double Distance_PointToHalfspace_i(PT_vector_T x, int i) {
		double a_dot_z_minus_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];

		if (PD_norm_a[i] < PP_EPS_ZERO) //Degenerate equation
			return 0;

		if (a_dot_z_minus_b < 0) // Point belongs to halfspace
			return 0;

		return a_dot_z_minus_b / PD_norm_a[i];
	}

	inline double Distance_PointToHyperplane_i(PT_vector_T x, int i) {
		if (PD_norm_a[i] < PP_EPS_ZERO) //Degenerate equation
			return 0;
		else
			return fabs(Vector_DotProduct(PD_A[i], x) - PD_b[i]) / PD_norm_a[i];
	}

	inline double Distance_PointToPoint(PT_vector_T x, PT_vector_T y) {
		PT_vector_T z;
		Vector_Subtraction(x, y, z);
		return Vector_Norm(z);
	}

	inline double Distance_PointToPolytope(PT_vector_T x) { // Measure of distance from point to polytope
		double maxDistance = 0;
		double distance;

		for (int i = 0; i < PD_m; i++) {
			distance = Distance_PointToHalfspace_i(x, i);
			if (distance > 0)
				maxDistance = PF_MAX(maxDistance, distance);
		}
		return maxDistance;
	}

	inline void MakeColumnOfNorms(PT_matrix_T A, PT_column_T norm_a) {
		for (int i = 0; i < PD_m; i++)
			norm_a[i] = Vector_Norm(A[i]);
	}

	inline void MakeListOfNotIncludingHalfspaces(PT_vector_T x, int* notIncludingHalfspacesList, double eps) {
		int mo = 0;
		for (int i = 0; i < PD_m; i++)
			if (!PointBelongsHalfspace_i(x, i, eps)) {
				notIncludingHalfspacesList[mo] = i;
				mo++;
			}
		if (mo < PD_m)
			notIncludingHalfspacesList[mo] = -1;
	}

	inline void MakePointHyperplaneList(PT_vector_T u, int* pointHyperplaneList, int* mh, double eps) {
		*mh = 0;
		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsHyperplane_i(u, i, eps)) {
				pointHyperplaneList[*mh] = i;
				(*mh)++;
			}
		}
	}

	inline void MovingOnPolytope(PT_vector_T startPoint, PT_vector_T directionVector, PT_vector_T finishPoint, double epsMoving) {
		double leftBound = 0;
		double rightBound = PP_DBL_MAX;
		double factor = 1;
		double delta;

		assert(Vector_Norm(directionVector) >= PP_EPS_ZERO);

		delta = factor / 2;

		while (rightBound - leftBound >= PP_EPS_ZERO && delta > 0) {
			Shift(startPoint, directionVector, factor, finishPoint);
			if (PointBelongsPolytope(finishPoint, PP_EPS_POINT_IN_HALFSPACE)) {
				leftBound = factor;
				delta *= 2;
				factor += delta;
			}
			else {
				rightBound = factor;
				delta /= 2;
				factor -= delta;
			}
		}

		Shift(startPoint, directionVector, factor, finishPoint);
		delta = epsMoving;
		while (!PointBelongsPolytope(finishPoint, epsMoving) && delta > 0) {
			factor -= delta;
			delta *= 2;
			Shift(startPoint, directionVector, factor, finishPoint);
		}
	}

	inline void MovingToPolytope(PT_vector_T startPoint, PT_vector_T directionVector, PT_vector_T finishPoint, double epsMoving) {
		double leftBound = 0;
		double rightBound = PP_DBL_MAX;
		double factor = 1;
		double delta;
		static int outerHalspace_i[PP_MM];	// Index of out half-spaces
		int mo;								// Number of out half-spaces
		bool pointInsideCone;

		assert(Vector_Norm(directionVector) >= PP_EPS_ZERO);

		mo = 0;
		for (int i = 0; i < PD_m; i++)
			if (!PointBelongsHalfspace_i(startPoint, i, PP_EPS_POINT_IN_HALFSPACE)) {
				outerHalspace_i[mo] = i;
				mo++;
			}

		delta = factor / 2;

		while (rightBound - leftBound >= PP_EPS_ZERO && delta > 0) {
			Shift(startPoint, directionVector, factor, finishPoint);

			pointInsideCone = true;
			for (int i = 0; i < mo; i++)
				if (PointBelongsHalfspace_i(finishPoint, outerHalspace_i[i], PP_EPS_POINT_IN_HALFSPACE)) {
					pointInsideCone = false;
					break;
				}
			if (pointInsideCone) {
				leftBound = factor;
				delta *= 2;
				factor += delta;
			}
			else {
				rightBound = factor;
				delta /= 2;
				factor -= delta;
				assert(factor > 0);
			}
		}

		Shift(startPoint, directionVector, factor, finishPoint);
		delta = epsMoving;
		do {
			pointInsideCone = false;
			for (int i = 0; i < mo; i++)
				if (!PointBelongsHalfspace_i(finishPoint, outerHalspace_i[i], epsMoving)) {
					pointInsideCone = true;
					factor -= delta;
					delta *= 2;
					assert(factor > 0);
					Shift(startPoint, directionVector, factor, finishPoint);
					break;
				}
		} while (pointInsideCone && delta > 0);
	}
	inline double ObjF(PT_vector_T x) {
		double s = 0;
		for (int j = 0; j < PD_n; j++)
			s += PD_c[j] * x[j];
		return s;
	}

	inline void	ObliqueProjectingVectorOntoHalfspace_i(PT_vector_T z, int i, PT_vector_T g, PT_vector_T o, double eps, int* exitCode) {
		// Oblique projecting vector o of point z onto Half-space H_i with respect to vector g
		double a_dot_g;	// <a,g>
		double a_dot_z_minus_b;	// <a,z> - b
		double factor;	// (b - <a,z>) / <a,g>

		if (PD_norm_a[i] < PP_EPS_ZERO) {
			Vector_Zeroing(o);
			*exitCode = 0; // #define PP_EXITCODE_DEGENERATE_INEQUALITY		0
			return;
		}

		a_dot_z_minus_b = Vector_DotProduct(PD_A[i], z) - PD_b[i]; // <a,z> - b

		if (fabs(a_dot_z_minus_b) < PP_EPS_ZERO) { // <a,z> - b < 0
			*exitCode = 1; // #define PP_EXITCODE_ON_HYPERPLANE				1
			Vector_Zeroing(o);
			return;
		}

		if (a_dot_z_minus_b <= -PP_EPS_ZERO) { // <a,z> - b < 0
			*exitCode = 2; // #define PP_EXITCODE_INSIDE_HALFSPACE			2
			Vector_Zeroing(o);
			return;
		}

		a_dot_g = Vector_DotProduct(PD_A[i], g); // <a,g>


		if (fabs(a_dot_g) < PP_EPS_ZERO) {
			*exitCode = 4; // #define PP_EXITCODE_PARALLEL				4
			Vector_Zeroing(o);
			return;
		}

		if (a_dot_g >= PP_EPS_ZERO) {
			*exitCode = 5; // #define PP_EXITCODE_RECESSIVE				5
			Vector_Zeroing(o);
			return;
		}

		factor = a_dot_z_minus_b / a_dot_g; // (<a,z> - b) / <a,g>

		// Oblique projection vector: o = -(<a,z> - b)g/<a, g> = -factor * g
		Vector_MultiplyByNumber(g, -factor, o);

		*exitCode = 9; // #define PP_EXITCODE_NONDEGENERATE_PROJECTING	9
		return;
	}

	inline void OrthogonalProjectingVectorOntoHalfspace_i(PT_vector_T z, int i, PT_vector_T r, double eps, int* exitCode) {
		double factor;
		double a_dot_z_minus_b = Vector_DotProduct(PD_A[i], z) - PD_b[i]; // <a,z> - b
		double distance = fabs(a_dot_z_minus_b) / PD_norm_a[i];

		if (PD_norm_a[i] < PP_EPS_ZERO) {
			Vector_Zeroing(r);
			*exitCode = 0; // #define PP_EXITCODE_DEGENERATE_INEQUALITY		0
			return;
		}

		if (distance < eps) {
			Vector_Zeroing(r);
			*exitCode = 1; // #define PP_EXITCODE_ON_HYPERPLANE				1
			return;
		}

		if (a_dot_z_minus_b < 0) { // <a,z> - b < 0
			Vector_Zeroing(r);
			*exitCode = 2; // #define PP_EXITCODE_INSIDE_HALFSPACE			2
			return;
		}

		factor = -a_dot_z_minus_b / (PD_norm_a[i] * PD_norm_a[i]); // (b - <z,a>) / ||a||^2
		Vector_MultiplyByNumber(PD_A[i], factor, r); // r = a(b - <z,a>) / ||a||^2
		*exitCode = 9; // #define PP_EXITCODE_NONDEGENERATE_PROJECTING	9
	}

	inline void OrthogonalProjectingVectorOntoHyperplane_i(PT_vector_T x, int i, PT_vector_T p) {
		assert(Vector_NormSquare(PD_A[i]));
		Vector_MultiplyByNumber(PD_A[i], -(Vector_DotProduct(PD_A[i], x) - PD_b[i]) / Vector_NormSquare(PD_A[i]), p);
	}

	inline bool PointBelongsHalfspace_i(PT_vector_T x, int i, double eps) {
		if (PD_norm_a[i] < eps) //Degenerate equation
			return true;
		double a_dot_x_minus_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];
		double distanceToHyperplane = fabs(a_dot_x_minus_b) / PD_norm_a[i];
		if (distanceToHyperplane < eps)
			return true;
		if (a_dot_x_minus_b < 0)
			return true;
		return false;
	}

	inline bool PointBelongsHyperplane_i(PT_vector_T x, int i, double eps) {
		if (Distance_PointToHyperplane_i(x, i) < eps)
			return true;
		else
			return false;
	}

	inline bool PointBelongsPolytope(PT_vector_T x, double eps) { // If the point belongs to the polytope with prescigion of eps
		for (int i = 0; i < PD_m; i++)
			if (!PointBelongsHalfspace_i(x, i, eps))
				return false;
		return true;
	}

	inline bool PointBelongsOuterCone(PT_vector_T x, int* notIncludingHalfspacesList, double eps) { // If the point belongs to the outer cone with prescigion of eps
		for (int i = 0; i < PD_m && notIncludingHalfspacesList[i] >= 0; i++)
			if (PointBelongsHalfspace_i(x, i, eps))
				return false;
		return true;
	}

	inline void PointHomothety(PT_vector_T x, PT_vector_T center, double ratio) { // https://en.wikipedia.org/wiki/Homothety
		if (ratio == 1)
			return;
		assert(ratio > 0);
		for (int j = 0; j < PD_n; j++)
			x[j] = ratio * x[j] - (ratio - 1) * center[j];
	}

	inline bool PointInsideHalfspace_i(PT_vector_T x, int i, double eps) {
		if (PD_norm_a[i] < PP_EPS_ZERO) //Degenerate equation
			return true;
		double a_dot_x_minus_b = Vector_DotProduct(PD_A[i], x) - PD_b[i];
		double distanceToHyperplane = fabs(a_dot_x_minus_b) / PD_norm_a[i];
		if (distanceToHyperplane < eps)
			return false;
		if (a_dot_x_minus_b < 0)
			return true;
		return false;
	}

	inline void PolytopeHomothety(PT_vector_T center, double ratio) { // https://en.wikipedia.org/wiki/Homothety
		if (ratio == 1)
			return;
		assert(ratio > 0);

		for (int i = 0; i < PD_m; i++) {
			PD_b[i] = ratio * PD_b[i] - (ratio - 1) * Vector_DotProduct(PD_A[i], center);
		}
	}

	inline void Print_Inequalities() {
		for (int i = 0; i < PD_m; i++) {
			cout << i << ")";
			for (int j = 0; j < PD_n; j++)
				cout << setw(PP_SETW) << PD_A[i][j];
			cout << "\t<=" << setw(PP_SETW) << PD_b[i] << endl;
		}
	}

	inline void Print_HalfspacesIncludingPoint(PT_vector_T x, double eps) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsHalfspace_i(x, i, eps)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	inline void Print_HalfspacesOutOfPoint(PT_vector_T x, double eps) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (!PointBelongsHalfspace_i(x, i, eps)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	inline void Print_HyperplanesIncludingPoint(PT_vector_T x, double eps) {
		bool comma = false;

		cout << "{";

		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsHyperplane_i(x, i, eps)) {
				if (comma)
					cout << ", ";
				else
					comma = true;
				cout << i;
			}
		}

		cout << "}";
	}

	inline void Print_Vector(PT_vector_T x) {
		for (int j = 0; j < PF_MIN(PP_OUTPUT_LIMIT, PD_n); j++) cout << setw(PP_SETW) << x[j];
		if (PP_OUTPUT_LIMIT < PD_n) cout << "	...";
	}

	inline double RelativeError(double trueValue, double calculatedValue) {
		if (fabs(trueValue) >= PP_EPS_ZERO)
			return fabs(calculatedValue - trueValue) / fabs(trueValue);
		else
			return fabs(calculatedValue - trueValue);
	}

	inline void Shift(PT_vector_T point, PT_vector_T shiftVector, double factor, PT_vector_T shiftedPoint) {
		for (int j = 0; j < PD_n; j++)
			shiftedPoint[j] = point[j] + shiftVector[j] * factor;
	}

	inline void Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x + y
		for (int j = 0; j < PD_n; j++)
			z[j] = x[j] + y[j];
	}

	inline void Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint) { // toPoint = fromPoint
		for (int j = 0; j < PD_n; j++)
			toPoint[j] = fromPoint[j];
	}

	inline void Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = x/r
		for (int j = 0; j < PD_n; j++)
			y[j] = x[j] / r;
	}

	inline void Vector_DivideEquals(PT_vector_T x, double r) {  // x = x/r
		for (int j = 0; j < PD_n; j++)
			x[j] /= r;
	}

	inline double Vector_DotProduct(PT_vector_T x, PT_vector_T y) {
		double sum = 0;
		for (int j = 0; j < PD_n; j++)
			sum += x[j] * y[j];
		return sum;
	}

	inline bool Vector_Is_Tiny(PT_vector_T x, double eps) {
		return Vector_Norm(x) < eps;
	}

	inline void Vector_MakeLike(PT_vector_T x, double lengthOfLikeVector, PT_vector_T likeVector) {
		double norm_x = Vector_Norm(x);
		if (norm_x == 0)
			Vector_Zeroing(likeVector);
		else
			Vector_MultiplyByNumber(x, lengthOfLikeVector / norm_x, likeVector);
	}

	void Vector_MakeMinus_e(PT_vector_T minus_e) {
		for (int j = 0; j < PD_n; j++)
			minus_e[j] = -1;
	}

	inline void Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector) { // equalPoint += minusVector
		for (int j = 0; j < PD_n; j++)
			equalPoint[j] -= minusVector[j];
	}

	inline void Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y) {  // y = r*x
		for (int j = 0; j < PD_n; j++)
			y[j] = x[j] * r;
	}

	inline void Vector_MultiplyEquals(PT_vector_T x, double r) {  // x = r*x
		for (int j = 0; j < PD_n; j++)
			x[j] *= r;
	}

	inline double Vector_Norm(PT_vector_T x) {
		double norm_x = sqrt(Vector_NormSquare(x));
		return norm_x;
	}

	inline double Vector_NormSquare(PT_vector_T x) {
		double sum = 0;

		for (int j = 0; j < PD_n; j++) {
			sum += x[j] * x[j];
		}
		return sum;
	}

	inline void Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector) { // equalVector += plusVector
		for (int j = 0; j < PD_n; j++)
			equalVector[j] += plusVector[j];
	}

	inline void Vector_Round(PT_vector_T x, double eps) {
		double floorValue;
		double fractionalPart;
		double sign;
		double absValue;

		for (int j = 0; j < PD_n; j++) {
			if (fabs(x[j]) < eps) {
				x[j] = 0;
				continue;
			}
			absValue = fabs(x[j]);
			sign = x[j] > 0 ? 1 : -1;
			floorValue = floor(absValue);
			fractionalPart = absValue - floorValue;
			if (1 - fractionalPart < eps) {
				x[j] = sign * (floorValue + 1);
				continue;
			}
			if (fractionalPart < eps)
				x[j] = sign * floorValue;
		}
	}

	inline void Vector_SetValue(PT_vector_T x, double v) {  // x = (v,...,v)
		for (int j = 0; j < PD_n; j++) x[j] = v;
	}

	inline void Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z) {  // z = x - y
		for (int j = 0; j < PD_n; j++)
			z[j] = x[j] - y[j];
	}

	inline void Vector_Zeroing(PT_vector_T x) {  // x = 0
		for (int j = 0; j < PD_n; j++) x[j] = 0;
	}
}

//---------------------------------- Private Functions -------------------------
namespace PF {
	using namespace SF;

	inline void CodeToSubset(int code, int subset[PP_MM], int* ma) {
		*ma = 0;
		for (int i = 0; i < PD_mh; i++) {
			if (code % 2 == 1) {
				subset[*ma] = PD_pointHyperplanes[i];
				(*ma)++; assert(*ma <= PP_MM);
			}
			code /= 2;
			if (code == 0)
				break;
		}
	}

	inline void MakeFaceList(int* faceCodeList, int K) {
		int index;

		for (int k = 0; k < PP_KK; k++) {
			faceCodeList[k] = 0;
		}

		if (PP_KK <= (PP_RND_MAX + 1)) {
			for (int k = 0; k < K; k++) {
				index = rand() % PP_KK;
				if (faceCodeList[index] == 0)
					faceCodeList[index] = k;
				else {
					for (int ki = 1; ki < K; ki++) {
						if (faceCodeList[(index + ki) % PP_KK] == 0) {
							faceCodeList[(index + ki) % PP_KK] = k;
							break;
						}
					}
				}
			}
			return;
		}

		assert(K <= PP_KK);

		int segmentCount = PP_KK / (PP_RND_MAX + 1);
		assert(PP_KK % (PP_RND_MAX + 1) == 0);

		if (K < segmentCount) {
			for (int k = 1; k < K; k++)
				faceCodeList[(k - 1) * (PP_RND_MAX + 1)] = k;
			return;
		}

		int segmentK = K / segmentCount;
		assert(K % segmentK == 0);

		for (int l = 0; l < segmentCount; l++) {
			for (int k = l * segmentK; k < (l + 1) * segmentK; k++) {
				index = rand() % (PP_RND_MAX + 1) + l * (PP_RND_MAX + 1);
				if (faceCodeList[index] == 0)
					faceCodeList[index] = k;
				else {
					for (int ki = 1; ki < (PP_RND_MAX + 1); ki++) {
						if (faceCodeList[l * (PP_RND_MAX + 1) + (index + ki) % (PP_RND_MAX + 1)] == 0) {
							faceCodeList[l * (PP_RND_MAX + 1) + (index + ki) % (PP_RND_MAX + 1)] = k;
							break;
						}
					}
				}
			}
		}
	}

	inline void PreparationForIteration(PT_vector_T u) {
#ifdef PP_DEBUG
		cout << "[" << BSF_sv_mpiRank << "] : prepare for iteration point = ";
		Print_Vector(u);
		cout << endl;
		if (!PointBelongsPolytope(u, PP_EPS_POINT_IN_HALFSPACE)) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster) {
				cout << "\nPoint does not belong to the feasible polytope with precision PP_EPS_ZERO = "
					<< PP_EPS_POINT_IN_HALFSPACE << "!!! You should decrease this parameter.\n";
				abort();
			}
		}
#endif PP_DEBUG

		MakePointHyperplaneList(u, PD_pointHyperplanes, &PD_mh, PP_EPS_POINT_IN_HALFSPACE);
		assert(PD_mh <= PP_MM);

		if (!PointBelongsPolytope(u, PP_EPS_POINT_IN_HALFSPACE)) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster) {
				cout << "\nPoint does not belong to the feasible polytope with precision PP_EPS_ZERO = "
					<< PP_EPS_POINT_IN_HALFSPACE << "!!! You should decrease this parameter.\n";
				exit(1);
			}
		}

		PD_K = (int)pow(2, (double)PD_mh);
		if (PD_K > PP_KK) {
			if (BSF_sv_mpiRank == BSF_sv_mpiMaster)
				cout << "Parameter PP_KK = " << PP_KK << " must be greater than or equal to " << PD_K << "\n";
			exit(1);
		}

		//
		MakeFaceList(PD_faceCodes, PD_K);
		PD_objF_u = ObjF(PD_u);
		PD_objF_initialValue = PD_objF_u;
		PD_shiftLength = 0;
	}

	inline void Print_Number_of_faces(PT_vector_T x) {
		int mh;
		int K;

		mh = 0;
		for (int i = 0; i < PD_m; i++) {
			if (PointBelongsHyperplane_i(x, i, PP_EPS_POINT_IN_HALFSPACE))
				mh++;
		}

		//
		//
		//
		//
		//
		//
		//
		//
		//
		// 
		// 

		K = (int)pow(2, (double)mh) - 1;
		cout << K << endl;
		//
	}

	static inline void PseudoprojectionOnFace(PT_vector_T v, PT_vector_T w, double eps, int* success) {
		PT_vector_T r;
		PT_vector_T w_previous;
		double dist;
		int iterCount = 0;

		Vector_Copy(v, w);

		do {
			Vector_Zeroing(r);
			Vector_Copy(w, w_previous);

			for (int i = 0; i < PD_ma; i++) {
				int ia = PD_faceHyperplanes[i];
				PT_vector_T p;
				OrthogonalProjectingVectorOntoHyperplane_i(w, ia, p);
				Vector_PlusEquals(r, p);
			}

			Vector_DivideEquals(r, PD_ma);
			Vector_Round(r, eps);
			Vector_PlusEquals(w, r);

			dist = Distance_PointToPoint(w, w_previous);
			iterCount++;
			if (iterCount > PP_MAX_PROJECTING_ITER) {
				*success = false;
				break;
			}
		} while (dist >= eps / 10);
	}
	//---------------------------------- Old Problem functions -------------------------
#ifdef PP_DATABASE_OUTPUT
//void SavePoint() {
//	vector<double> _buff(PD_n);
//	for (unsigned i = 0; i < PD_n; i++)
//		_buff[i] = PD_u[i];
//
//	PD_newPoint.id = 0;
//	PD_newPoint.problem_id = PD_problem.id;
//	PD_newPoint.coefficients = doubleToChar(_buff);
//	PD_newSurfacePoints.push_back(PD_newPoint);
//	//if (newPoint.id)
//	//	cout << "Calculated surface point is saved into the database!" << endl;
//}

	void SavePrecedent(PT_bsf_reduceElem_T* reduceResult) {
		vector<double> _buff(PD_n);

		CodeToSubset(reduceResult->faceCode, PD_faceHyperplanes, &PD_ma);

		for (unsigned i = 0; i < PD_n; i++)
			_buff[i] = reduceResult->face_count[i];
		PD_precedent.face_count = doubleToChar(_buff);

		PD_precedent.id = 0;

		PD_precedent.problem_id = PD_problem.id;

		for (unsigned i = 0; i < PD_n; i++)
			_buff[i] = PD_previous_u[i];
		PD_precedent.coefficients = doubleToChar(_buff);

		for (unsigned i = 0; i < PD_n; i++)
			_buff[i] = reduceResult->d[i];

		PD_precedent.d = doubleToChar(_buff);

		PD_precedent.shift = PD_shiftLength;

		//		PD_precedent.face = PD_n - PD_ma;
		PD_precedent.face = PD_ma;
		_buff.clear();
		for (unsigned i = 0; i < PD_ma; i++)
			_buff.push_back(PD_faceHyperplanes[i]);
		PD_precedent.face_numbers = doubleToChar(_buff);

		PD_newPrecedents.push_back(PD_precedent);

		//if(PD_precedent.id)
		//	cout << "New precedent is saved into the database!" << endl;

	}

	std::vector<double> charToDouble(std::vector<char> _In) {
		std::vector<double> _Out;
		size_t size = _In.size();
		_Out.resize(size / sizeof(double));
		std::memcpy(_Out.data(), _In.data(), size);
		return _Out;
	}

	std::vector<char> doubleToChar(std::vector<double> _In) {
		std::vector<char> _Out;
		size_t size = _In.size() * sizeof(double);
		_Out.resize(size);
		std::memcpy(_Out.data(), _In.data(), size);
		return _Out;
	}

	void printLppForm(Problem problem, std::vector<Inequality> inequalities) {
		int M = 2 * problem.N + inequalities.size();
		int width = 10;
		std::vector<double> _vec;
		std::cout << problem.id << '\t' << problem.N << '\t' << M << std::endl;
		for (int i = 0; i < problem.N; i++) {
			for (int j = 0; j < problem.N; j++)
				if (i == j) std::cout << std::setw(width) << double(1);
				else std::cout << std::setw(width) << double(0);
			std::cout << std::setw(width) << problem.high << std::endl;
		}
		for (int i = 0; i < inequalities.size(); i++) {
			_vec = charToDouble(inequalities[i].coefficients);
			for (int j = 0; j < problem.N; j++)
				std::cout << std::setw(width) << _vec[j];
			std::cout << std::setw(width) << inequalities[i].b << std::endl;
		}
		for (int i = 0; i < problem.N; i++) {
			for (int j = 0; j < problem.N; j++)
				if (i == j) std::cout << std::setw(width) << double(-1);
				else std::cout << std::setw(width) << double(0);
			std::cout << std::setw(width) << problem.low << std::endl;
		}
		_vec = charToDouble(problem.c);
		std::copy(_vec.begin(), _vec.end(), std::ostream_iterator<double>(std::cout, "\t"));
		std::cout << std::endl;
	}
#endif // PP_DATABASE_OUTPUT
}
