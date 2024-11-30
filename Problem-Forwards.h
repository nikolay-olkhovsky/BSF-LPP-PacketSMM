/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP AlfaMove (surface movement method)
Module: Problem-Forwards.h (Problem Function Forwards)
Author: Nikolay A. Olkhovsky
This source code has been produced with using BSF-skeleton (https://github.com/leonid-sokolinsky/BSF-skeleton)
==============================================================================*/
#include "Problem-bsfTypes.h"
#include "Problem-Types.h"
//====================== Private Functions ===========================
namespace PF {
	void	CodeToSubset(int code, int subset[PP_MM], int* ma);
	//
	void	MakeFaceList(int* faceCodeList, int K);
//	void	MakeHyperplaneList(PT_vector_T u, int* index_includingHyperplanes, int* mh, double eps);
	void	PreparationForIteration(PT_vector_T u);
	void	Print_Number_of_faces(PT_vector_T x);
	void	PseudoprojectionOnFace(PT_vector_T v, PT_vector_T w, double eps, int* success);
//	void	SavePoint();
	void	SavePrecedent(PT_bsf_reduceElem_T* reduceResult);

#ifdef PP_DATABASE_OUTPUT
	std::vector<double> charToDouble(std::vector<char> _In);
	std::vector<char> doubleToChar(std::vector<double> _In);
	void printLppForm(Problem problem, std::vector<Inequality> inequalities);
#endif //PP_DATABASE_OUTPUT
}
//====================== Shared Functions ===========================
namespace SF {
//	void	Conversion();
//	void	ConversionSimple();
	double	Distance_PointToHalfspace_i(PT_vector_T x, int i);
	double	Distance_PointToHyperplane_i(PT_vector_T x, int i);
	double	Distance_PointToPoint(PT_vector_T x, PT_vector_T y);
	double	Distance_PointToPolytope(PT_vector_T x);
	void	MakeColumnOfNorms(PT_matrix_T A, PT_column_T norm_a);
	void	MakeListOfNotIncludingHalfspaces(PT_vector_T x, int* notIncludingHalfspacesList, double eps);
	void	MakePointHyperplaneList(PT_vector_T u, int* pointHyperplaneList, int* mh, double eps);
	void	MovingOnPolytope(PT_vector_T startPoint, PT_vector_T directionVector, PT_vector_T finishPoint, double epsMoving);
	void	MovingToPolytope(PT_vector_T startPoint, PT_vector_T directionVector, PT_vector_T finishPoint, double epsMoving);
	double	ObjF(PT_vector_T x);
	void	ObliqueProjectingVectorOntoHalfspace_i(PT_vector_T z, int i, PT_vector_T g, PT_vector_T o, double eps, int* exitCode);
	void	OrthogonalProjectingVectorOntoHalfspace_i(PT_vector_T z, int i, PT_vector_T r, double eps, int* exitcode);
	void	OrthogonalProjectingVectorOntoHyperplane_i(PT_vector_T x, int i, PT_vector_T p);
	bool	PointBelongsHalfspace_i(PT_vector_T point, int i, double eps);
	bool	PointBelongsHyperplane_i(PT_vector_T z, int i, double eps);
	bool	PointBelongsOuterCone(PT_vector_T x, int* notIncludingHalfspacesList, double eps);
	bool	PointBelongsPolytope(PT_vector_T x, double eps);
	void	PointHomothety(PT_vector_T x, PT_vector_T center, double ratio);
	bool	PointInsideHalfspace_i(PT_vector_T x, int i, double eps);
	void	PolytopeHomothety(PT_vector_T center, double ratio);
	void	Print_Inequalities();
	void	Print_HalfspacesIncludingPoint(PT_vector_T x, double eps);
	void	Print_HalfspacesOutOfPoint(PT_vector_T x, double eps);
	void	Print_HyperplanesIncludingPoint(PT_vector_T x, double eps);
	void	Print_Vector(PT_vector_T x);
	double	RelativeError(double trueValue, double calculatedValue);
	//void	RemoveFreeVariables(int m_equation, int m_inequality, int m_lowerBound, int m_higherBound);
	void	Shift(PT_vector_T point, PT_vector_T shiftVector, double factor, PT_vector_T shiftedPoint);
	void	Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z);
	void	Vector_Copy(PT_vector_T x, PT_vector_T y);
	void	Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y);
	void	Vector_DivideEquals(PT_vector_T x, double r);
	double	Vector_DotProduct(PT_vector_T x, PT_vector_T y);
	bool	Vector_Is_Tiny(PT_vector_T x, double eps);
	void	Vector_MakeLike(PT_vector_T x, double lengthOfLikeVector, PT_vector_T likeVector);
	void	Vector_MakeMinus_e(PT_vector_T minus_e);
	void	Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector);
	void	Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y);
	void	Vector_MultiplyEquals(PT_vector_T x, double r);
	double	Vector_Norm(PT_vector_T x);
	double	Vector_NormSquare(PT_vector_T x);
	void	Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector);
	void	Vector_Round(PT_vector_T x, double eps);
	void	Vector_SetValue(PT_vector_T x, double v);
	void	Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z);
	void	Vector_Zeroing(PT_vector_T x);
////====================== Old Problem Functions ===========================
//void	AddOppositeInequality(int hyperplaneIndex, int m);
//void	CodeToSubset(int code, int subset[PP_MM], int* ma);
//bool	Conversion();
//void	DirVectorCleanup(PT_vector_T x, double eps);
//double	Distance(PT_vector_T x, PT_vector_T y);
//void	MakeHyperplaneList(int* mh);
//bool	MakeHyperplaneSubsetCodeList(int mh, int* K);
//void	MakeObjVector(PT_vector_T c, double length, PT_vector_T objVector);
//bool	MovingOnSurface(PT_vector_T directionVector, PT_vector_T point);
//bool	MTX_Load__Problem();
//bool	MTX_Load_A(int* nor, int* noc, int* non, int* noe);
//bool	MTX_Load_b(int* nor, int* noc, int* noe);
//bool	MTX_Load_c(int* nor, int* noc);
//bool	MTX_Load_hi(int* nor, int* noc);
//bool	MTX_Load_lo(int* nor, int* noc);
//bool	MTX_Load_sp(int* nor, int* noc);
//bool	MTX_Save_sp(PT_vector_T x, double elapsedTime);
//double	ObjF(PT_vector_T x);
//bool	PointInHalfspace(PT_vector_T point, PT_vector_T a, double b, double eps);
//bool	PointInPolytope(PT_vector_T x);
//double	PolytopeResidual(PT_vector_T x);
//void	Print_VectorOnActiveHyperplanes(PT_vector_T x);
//void	Print_VectorOnHyperplanes(PT_vector_T x);
//double	ProblemScale();
//void	PseudoprojectionOnPolytope(PT_vector_T v, PT_vector_T w);
//void	PseudorojectionOnFace(PT_vector_T v, PT_vector_T w, double eps);
//double	relativeError(double trueValue, double calcValue);
//void	Shift(PT_vector_T basePoint, PT_vector_T direction, double PD_shiftLength, PT_vector_T endPoint);
//void	SkipComments(FILE* stream);
//void	Vector_Addition(PT_vector_T x, PT_vector_T y, PT_vector_T z);
//void	Vector_Copy(PT_vector_T fromPoint, PT_vector_T toPoint);
//double	Vector_DistanceToHalfspace(PT_vector_T z, PT_vector_T a, double b);
//void	Vector_DivideByNumber(PT_vector_T x, double r, PT_vector_T y);
//void	Vector_DivideEquals(PT_vector_T x, double r);
//double	Vector_DotProduct(PT_vector_T x, PT_vector_T y);
//bool	Vector_Equal(PT_vector_T x, PT_vector_T y, double eps);
//double	Vector_Norm(PT_vector_T x);
//double	Vector_NormSquare(PT_vector_T x);
//void	Vector_MinusEquals(PT_vector_T equalPoint, PT_vector_T minusVector);
//void	Vector_MultiplyByNumber(PT_vector_T x, double r, PT_vector_T y);
//void	Vector_MultiplyEquals(PT_vector_T x, double r);
//void	Vector_PlusEquals(PT_vector_T equalVector, PT_vector_T plusVector);
//void	Vector_ObliqueProjectionOntoHalfspace(PT_vector_T z, PT_vector_T a, double b, PT_vector_T g, PT_vector_T o, int* exitCode);
//bool	Vector_OnHyperplane(PT_vector_T point, PT_vector_T a, double b, double eps, double* residual);
//double	Vector_OrthogonalProjectionOntoHalfspace(PT_vector_T z, PT_vector_T a, double b, PT_vector_T r, double eps, int* exitCode);
//void	Vector_Round(PT_vector_T x, double eps);
//void	Vector_Subtraction(PT_vector_T x, PT_vector_T y, PT_vector_T z);
//void	Vector_Unit(PT_vector_T vector);
//void	Vector_Zero(PT_vector_T x);
}
//====================== Macros ================================
#define PF_MIN(x,y) (x<y?x:y)
#define PF_MAX(x,y) (x>y?x:y)
#define PF_MAP_LIST_INDEX (BSF_sv_addressOffset + BSF_sv_numberInSublist)