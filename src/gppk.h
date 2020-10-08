/* gppk.h */

/******************************************************************************
 *  This code is part of GPPK (The Generalized Pooling Problem Kit).
 *
 *  Copyright (C) 2009, 2010 Mohammed Alfaki, Department of Informatics,
 *  University of Bergen, Bergen, Norway. All rights reserved. E-mail:
 *  <mohammeda@ii.uib.no>.
 *
 *  GPPK is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 ******************************************************************************/

#ifndef GPPK_H
#define GPPK_H

#include <limits.h>
#include <float.h>
#include <list>
#include <vector>
#include <iostream>
#include <signal.h>

#include <ilcplex/ilocplex.h>
ILOSTLBEGIN

#define inf IloInfinity
#define RGAP 1.0e-6
/* relative optimality gap */

#define AGAP 1.0e-6
/* absolute optimality gap */

#define TIME 3600
/* Maximum time allowable */

using namespace std;

/* library version numbers: */
#define GPP_MAJOR_VERSION 0
#define GPP_MINOR_VERSION 1

#define BRNCH 1
/* select branch rule: 1 min-max-infeas, 2 max-infeas */

static jmp_buf buf;
/* Jump buffer */

#define STRING_LENGTH 2048
/* maximum string length */

/******************************************************************************/
/* gpplib.cpp methods */

string get_str(const char *fmt, ...);
/* Write formatted output to char string */

void xerror(const char *fmt, ...);
/* print error message and terminate processing */

string file_name (const char *file);
/* return file name for a full path of file */

int str2int(const char *str, int *_val);
/* convert character string to value of int type */

int str2num(const char *str, double *_val);
/* convert character string to value of double type */

char *cdate();
/* return the current date and time char string */

string clocktime(double ttime);
/* convert time in seconds to hh:mm:ss:ms */

int n_choose_k(int n, int k);
/* compute n choose k */

/******************************************************************************/
/* gppbnd.cpp, Branch and bound methods */

/* node in the branch and bound tree */
class node {
public:
	double obj;
	/* ub on the obj */

	vector <double> *fl;
	/* variable lb */

	vector <double> *fu;
	/* solution values */

	vector <double> *wyl;
	/* variable ub */

	vector <double> *wyu;
	/* variable ub */

	node() { obj = -IloInfinity; };
	/* constructor */

	friend bool operator<(const node &n1, const node &n2);
	/* sort the nodes by its objective value */

	friend bool prune_by_bound(const node &n);
	/* prune a node by bound */
};

/******************************************************************************/
class vertex;
/* vertex class */

class arc;
/* arc class */

/* the generalized pooling problem class */
class graph {

	IloEnv env;
	/* cplex environment object */

	IloModel model;
	/* cplex model object */

	IloCplex cplex;
	/* cplex object */

public:
	FILE *fo;
	/* solution output file  */

	string name;
	/* name assigned to the graph */

	int na;
	/* number of arcs in the graph, na >= 0 */

	int nv;
	/* number of vertices in the graph, nv >= 0 */

	int ns;
	/* number of source vertices in the graph */

	int nt;
	/* number of terminal vertices in the graph */

	int nk;
	/* number of qualities in the problem */

	bool isSpp;
	/* check if the network is standard pooling problem */

	int mdl;
	/* formulation indicator: 1 MCF-model, 2 MCF^t-model,
	 * 3 STP-model and 4 P-model */

	int nlp_consts;
	/* number of original constraints */

	int lconsts;
	/* retrieve number of linear constraints */

	int nlconsts;
	/* retrieve number of nonlinear constraints */

	int vars;
	/* retrieve number of variables */

	int intvars;
	/* retrieve number of binary variables */

	int nlterms;
	/* retrieve number of distinct nonlinear terms */

	double ttime;
	/* total time elapsed */

	int solstat;
	/* optimizer status */

	int nnode;
	/* number of B&B nodes */

	double gpp_ub;
	/* upper for the optimal solution */

	double gpp_lb;
	/* lower for the optimal solution */

	list <node> tree;
	/* branch and bound tree */

	IloObjective obj;
	/* objective function object */

	IloRangeArray cons;
	/* constraints object */

	vector <int> *p_indeg;
	/* number of sources associated with each pool */

	vector <int> num_p;
	/* number of duplications to each pools */

	vector <double> *Y;
	/* discretized proportional variables */

	vector <vertex> N;
	/* set of all nodes in the graph */

	vector <arc> *A; /* A[nv];
    /* A[i], list of all arcs (i,j)\in A */

	int trngl[3];
	/* optimal rectangle for splitting */

	bool source_McCormick;
	/* the type of McCormick envelopes */

	int read_graph(const char *fname, int mod);
	/* read the pooling problem data in DIMACS format */

	int random_graph(int *val, double p, const char *fname);
	/* generate random pooling problem data */

	int write_graph(const char *fname);
	/* write the pooling problem data in DIMACS format */

	int write_gams(const char *fname);
	/* write the pooling problem data in GAMS file */

	int write_matlab(const char *fname, int flag);
	/* write the pooling problem data in Matlab format */

	int rnd_node(int i, int net);
	/* generate node j randomly to the arc (i,j) */

	int is_arc(int i, int j);
	/* check whether the arc (i,j) in the network or not */

	int i_p(int i, int j, int flg);
	/* return the position of variable in the extended network */

	int s_p(int s, int i);
	/* work the same as is_path(s, i) */

	int i_j(int j);
	/* return the parent of j */

	bool is_path(int s, int l);
	/* check whether there is a path between s and i or not */

	IloEnv get_env() const { return env; }
	/* get the cplex environment object */

	void build_relaxation();
	/* graph constructor */

	void def_objectivefunc();
	/* define the objective function */

	bool nlp_feasibility();
	/* check the solution of the current node if it feasible */

	int check_feasibility();
	/* check if the pooling problem is feasible or not */

	void branch_and_bound(clock_t stat);
	/* branch and bound */

	void set_bounds(node *curr_n, int itr);
	/* set the variable bounds */

	node create_child(node *curr_n, double obj);
	/* create a child node */

	void add_child_nodes(node &curr_n, double obj);
	/* create child nodes for the current node */

	int discretize_proportions(int n, int num_s, int j);
	/* discretize the proportional variables */

	double discrete_model(int n, int flg, double tm);
	/* solve the discrete model */

	vector<double> *one_terminal_mdl(int tau, list<int> Tp, double tm);
	/* solving the pooling problem with only one terminal */

	void printFp();

	void heuristic_alg();
	/* heuristic algorithm to find a feasible solution to the pooling problem */

	void problem_info();
	/* print problem information */

	void print_summary(int sgnl);
	/* print solution summary */

	graph(): env() {
		gpp_lb = -inf, gpp_ub = inf;
		solstat = 0, nnode = 0, vars = 0, nlterms = 0;
		lconsts = 0, nlconsts = 0, intvars = 0;
	};
	/* graph constructor */

	~graph() {
		delete[] A;
		if (mdl == 6) delete[] Y;
		delete[] p_indeg;
		env.end();
	};
	/* graph destructor */
};

/* vertex class */
class vertex {
public:
	graph *G;
	/* the owner graph of this vertex*/

	int v;
	/* vertex index */

	double bl;
	/* lower bound capacity associated with the vertex */

	double bu;
	/* upper bound capacity associated with the vertex */

	vector <double> q;
	/* the upper (and lower if any) bound quality vector for a terminal */

	IloNumVarArray wy;
	/* quality (or proportion) variables associated with the arc */

	void def_flowcapconsts();
	/* define the flow capacity constraints */

	void def_qualrequrmntconsts();
	/* define the quality requirement constraints */

	void def_rltflowcap();
	/* add the the RLT for the flow capacity */

	vertex(graph *g, int i, double l, double u);
	/* vertex constructor */

	vertex() { };
	/* vertex constructor */

	~vertex() { };
	/* vertex destructor */
};

/* arc class */
class arc {
public:
	graph *G;
	/* the owner graph of this arc*/

	int i;
	/* the head endpoint */

	int j;
	/* tail endpoint */

	double c;
	/* unit cost associated with the arc */

	IloNumVarArray fx;
	/* flow (and path flow) variable associated with the arc */

	string name;
	/* flow variable associated with the arc */

	IloNum F;

	void def_rltpropor();
	/* add the the RLT for the proportion balances */

	void def_McCormickenvlp();
	/* define the McCormick envelopes */

	int set_McCormickenvlp(int c);
	/* set the McCormick envelopes */

	arc(graph *g, int u, int v, double cst);
	/* arc constructor */

	arc() { };
	/* arc constructor */

	~arc() { };
	/* arc destructor */

};

/******************************************************************************/
/* gppmain.cpp */

void print_usage(const char* pkg);
/* print the package usage instructions */

/******************************************************************************/
#endif

/* eof */
