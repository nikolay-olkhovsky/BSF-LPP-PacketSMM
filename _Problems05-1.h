/*==============================================================================
Project: LiFe - New Linear Programming Solvers
Theme: Packet LPP AlfaMove (surface movement method)
Module: Problems05-1.h (LP problems of dimension 5 with 1 randome inequality)
Prefix: PP
Author: Nikolay A. Olkhovsky
This include file is part of Problem-Parameters.h
LP problems were obtained using BSF-LPP-Generator.
Initial surface points for these problems were calculated using BSF-Apex-Quest.
==============================================================================*/
#pragma once

//=========================== Method Parameters =========================
#define PP_EPS_ZERO						1E-8	// Accuracy for comparison with zero
#define PP_EPS_P_PROJ_ON_FASE			1E-12	// Precision for calculating pseudoprojection onto edge
#define PP_EPS_P_PROJ_ON_POLYTOPE		1E-9	// Precision for calculating pseudoprojection onto polytope
#define PP_OBJECTIVE_VECTOR_LENGTH		100000	// Starting length of Objective Vector
#define PP_EPS_MAKE_H_PLANE_LIST		1E-5	// Precision for MakeHyperplaneList()
#define PP_MAX_B_NO_CORRECT				200		// Maximum b that does not require correction
#define PP_RND_EPS_POINT_IN_POLYTOPE	1E-6	// Precision for random inequality in PointInPolytope()
#define PP_PROBE_LENGTH					0.003	// length of probe shift

#define PP_M 4		// Number of equations (number of rows in *.mtx)
#define PP_N 7		// Number of variables (number of cols in *.mtx)
#define PP_OPTIMAL_OBJ_VALUE 			0.0

/*============================== rnd5-0 LP problem ==============================*
// Start point:	200               0             0            0             0
// Exact solution:   100   200   200   200   200
// Face dimension: 4.      Generating hyperplanes: {0}.			Shift = 97979.59        F(x) = 2360
// Face dimension: 3.      Generating hyperplanes: {0, 4}.		Shift = 17950.549       F(x) = 2650
// Face dimension: 2.      Generating hyperplanes: {0, 3, 4}.	Shift = 12018.504       F(x) = 2780
// Face dimension: 2.      Generating hyperplanes: {3, 4, 5}.	Shift = 9428.0902       F(x) = 2820
// Face dimension: 1.      Generating hyperplanes: {2, 3, 4, 5}.Shift = 37712.362       F(x) = 2900
#define PP_PROBLEM_NAME	"rnd5-0"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 11		// Number of variables (number of cols in *.mtx)
#define PP_OPTIMAL_OBJ_VALUE 		2900
//------------------------------------------------------------------/**/

/*============================== rnd5-1-1 LP problem ==============================*
// Start point:	200               0               0             200             200
// Face dimension: 2.      Generating hyperplanes: {0, 3, 4}.      Shift = 26721.7		F(x) = 2579.2546
// Face dimension: 1.      Generating hyperplanes: {0, 3, 4, 5}.   Shift = 54967.646	F(x) = 2584.3495
#define PP_PROBLEM_NAME	"rnd5-1-1"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 11		// Number of variables (number of cols in *.mtx)
#define PP_OPTIMAL_OBJ_VALUE	2584.3495
//------------------------------------------------------------------/**/

/*============================== rnd5-1-2 LP problem ==============================*
// Start point:	0       199.99694       183.63843       188.16601       191.28208
// Face dimension: 3.      Generating hyperplanes: {5, 6}.		Shift = 1.0520607       F(x) = 2647.8187
// Face dimension: 2.      Generating hyperplanes: {1, 5, 6}.	Shift = 3497.6578       F(x) = 2652.9365
// Face dimension: 1.      Generating hyperplanes: {1, 4, 5, 6}.Shift = 12053.422       F(x) = 2657.5256
//#define PP_PROBLEM_NAME	"rnd5-1-2"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 11		// Number of variables (number of cols in *.mtx)
#define PP_OPTIMAL_OBJ_VALUE 			2657.5256116
//------------------------------------------------------------------/**/

/*============================== rnd5-1-3 LP problem ==============================*
// Start point:	200           200           200             0             0
// Face dimension: 2.      Generating hyperplanes: {0, 1, 2}.      Shift = 21948.028       F(x) = 2408.2173
// Face dimension: 1.      Generating hyperplanes: {0, 1, 2, 5}.   Shift = 8096.5775       F(x) = 2424.9192
#define PP_PROBLEM_NAME	"rnd5-1-3"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 11		// Number of variables (number of cols in *.mtx)
#define PP_OPTIMAL_OBJ_VALUE 2424.9191538
//------------------------------------------------------------------/**/

/*============================== rnd5-1-4 LP problem ==============================*
// Start point:	200             0             0           200           200
// Face dimension: 2.      Generating hyperplanes: {0, 3, 4}.      Shift = 22190.147       F(x) = 2274.477
// Face dimension: 1.      Generating hyperplanes: {0, 3, 4, 5}.   Shift = 6409.9716       F(x) = 2300.1128
#define PP_PROBLEM_NAME	"rnd5-1-4"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 11		// Number of variables (number of cols in *.mtx)
#define PP_OPTIMAL_OBJ_VALUE 2300.1127587
//------------------------------------------------------------------/**/

/*============================== rnd5-1-5 LP problem ==============================*
// Start point:	0           200           200             0             0
// Face dimension: 3.      Generating hyperplanes: {1, 2}.			Shift = 65129.188       F(x) = 2531.0733
// Face dimension: 2.      Generating hyperplanes: {1, 2, 5}.		Shift = 21689.588       F(x) = 2618.1681
// Face dimension: 1.      Generating hyperplanes: {1, 2, 5, 10}.	Shift = 27943.645       F(x) = 2626.4732
#define PP_PROBLEM_NAME	"rnd5-1-5"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 11		// Number of variables (number of cols in *.mtx)
#define PP_OPTIMAL_OBJ_VALUE 2626.4731647
//------------------------------------------------------------------/**/

/*============================== rnd5-1-6 LP problem ==============================*
// Start point:	200           200           200             0             0
// Face dimension: 2.      Generating hyperplanes: {0, 1, 2}.      Shift = 31072.332       F(x) = 2608.4395
// Face dimension: 1.      Generating hyperplanes: {0, 1, 2, 5}.   Shift = 54925.396       F(x) = 2675.352
#define PP_PROBLEM_NAME	"rnd5-1-6"
#define PP_M 6		// Number of equations (number of rows in *.mtx)
#define PP_N 11		// Number of variables (number of cols in *.mtx)
#define PP_OPTIMAL_OBJ_VALUE 2675.351984
//------------------------------------------------------------------/**/