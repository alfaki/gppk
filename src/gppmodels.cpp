/* gppmodels.cpp */

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

#include "gppk.h"
#include "gdxco.hpp"
#include <sys/time.h>

using namespace GAMS;

static string Indx[GMS_MAX_INDEX_DIM];
static gdxValues_t Values;

/* vertex constructor */
vertex::vertex(graph *g, int i, double l, double u): G(g), v(i), bl(l), bu(u) {
	if (!(1 <= G->mdl && G->mdl <= 5)) return;
	if (G->ns+1 <= v && v <= G->nv-G->nt) {
		IloEnv env = G->get_env();
		wy = IloNumVarArray(env);
		if (G->mdl == 5) {
			for (int k = 1; k <= G->nk; k++)
				wy.add(IloNumVar(env, -inf, inf, get_str("w(%d,%d)", k, v).c_str()));
		}
		else {
			for (int st = 1; st <= G->nv; st++)
				wy.add(IloNumVar(env, 0, 1, get_str("y(%d,%d)", st, v).c_str()));
		}
	}
}

/* arc constructor */
arc::arc(graph *g, int u, int v, double cst): G(g), i(u), j(v), c(cst) {
	if (G->mdl == 7) F = 0.0;
	if (!(1 <= G->mdl && G->mdl <= 5)) return;
	IloEnv env = G->get_env();
	fx = IloNumVarArray(env);
	name = get_str("(%d,%d)", i, j);
	fx.add(IloNumVar(env, 0, fmin(G->N[i].bu, G->N[j].bu), ("f"+name).c_str()));
	if (G->ns+1 <= i && i <= G->nv-G->nt) {
		if (G->mdl == 5) {
			for (int k = 1; k <= G->nk; k++)
				fx.add(IloNumVar(env, -inf, inf, get_str("x(%d,%d,%d)", k, i, j).c_str()));
		}
		else {
			for (int s = 1; s <= G->ns; s++) {
				double ub = fmin(G->N[s].bu, fmin(G->N[i].bu, G->N[j].bu));
				fx.add(IloNumVar(env, 0, ub, get_str("x(%d,%d,%d)", s, i, j).c_str()));
			}
		}
	}
}

/* graph constructor */
void graph::build_relaxation() {
	gpp_lb = -inf;
	gpp_ub = inf;
	obj = IloMinimize(env);
	def_objectivefunc();
	for (int i = 1; i <= nv; i++)
		N[i].def_flowcapconsts();
	for (int i = ns+1; i <= nv; i++)
		N[i].def_qualrequrmntconsts();
	for (int i = ns+1; i <= nv-nt; i++) {
		for (int j = 0; j < A[i].size(); j++)
			A[i][j].def_rltpropor();
	}
	for (int i = ns+1; i <= nv-nt; i++)
		N[i].def_rltflowcap();
	nlp_consts = cons.getSize();
	for (int i = ns+1; i <= nv-nt; i++) {
		for (int j = 0; j < A[i].size(); j++)
			A[i][j].def_McCormickenvlp();
	}
	model.add(cons);
	cplex.extract(model);
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	cplex.setError(env.getNullStream());
	vars = cplex.getNcols()-1;
	lconsts = cplex.getNrows()-4*nlterms;
	intvars = cplex.getNintVars();
	//cplex.exportModel(("lps/"+name+".lp").c_str());
}

/* define the objective function */
void graph::def_objectivefunc() {
	IloExpr expr(env);
	int nlt = 0;
	if (mdl == 4 || mdl == 5) { /* P- and MCF-formulation objective function */
		for (int i = 1; i <= nv-nt; i++) {
			for (int j = 0; j < A[i].size(); j++)
				expr += A[i][j].c*A[i][j].fx[0];
		}
		if (mdl == 4) {
			for (int s = 1; s <= ns; s++) {
				for (int i = ns+1; i <= nv-nt; i++) {
					if (s_p(s,i)>=0)
						nlt += A[i].size();
				}
			}
		}
	}
	else { /* PQ-, TP-, STP-formulation objective function */
		for (int s = 1; s <= ns; s++) {
			for (int i = ns+1; i <= nv-nt; i++) {
				for (int j = 0; j < A[i].size(); j++) {
					int ip = is_arc(s,i);
					if (ip>=0) {
						expr += (A[s][ip].c+A[i][j].c)*A[i][j].fx[s];
						nlt++;
					}
				}
			}
			for (int j = 0; j < A[s].size(); j++) {
				int t = A[s][j].j;
				if (nv-nt < t && t <= nv)
					expr += A[s][j].c*A[s][j].fx[0];
			}
		}
	}
	if (mdl != 5) {
		if (mdl == 3)
			nlt = 2*nlt;
		nlterms = nlt, nlconsts = nlt;
	}
	obj.setName("objective");
	obj.setExpr(expr);
	model.add(obj);
	expr.end();
}

/* define the flow capacity constraints */
void vertex::def_flowcapconsts() {
	IloEnv env = G->get_env();
	int c = G->cons.getSize();
	IloExpr expr(env);
	if (1 <= G->mdl && G->mdl <= 3) {
		if (v <= G->ns) {
			for (int l = 0; l < G->A[v].size(); l++) {
				int i = G->A[v][l].j;
				if (i <= G->nv-G->nt) {
					for (int j = 0; j < G->A[i].size(); j++)
						expr -= G->A[i][j].fx[v];
				}
				else
					expr -= G->A[v][l].fx[0];
			}
		}
		else if (v > G->nv-G->nt) {
			for (int s = 1; s <= G->ns; s++) {
				for (int j = G->ns+1; j <= G->nv-G->nt; j++) {
					int vp = G->is_arc(j,v);
					if (G->s_p(s,j)>=0 && vp>=0)
						expr -= G->A[j][vp].fx[s];
				}
				int vp = G->is_arc(s,v);
				if (vp>=0)
					expr -= G->A[s][vp].fx[0];
			}
		}
		else {
			for (int s = 1; s <= G->ns; s++) {
				if (G->s_p(s,v)>=0) {
					for (int j = 0; j < G->A[v].size(); j++)
						expr -= G->A[v][j].fx[s];
				}
			}
		}
	}
	if (G->mdl == 4 || G->mdl == 5) {
		if (v <= G->ns) {
			for (int j = 0; j < G->A[v].size(); j++)
				expr -= G->A[v][j].fx[0];
		}
		else {
			for (int j = 1; j <= G->nv; j++) {
				int vp = G->is_arc(j,v);
				if (vp>=0)
					expr -= G->A[j][vp].fx[0];
			}
		}
	}
	if (bl > 0 && !(G->ns < v && v <= G->nv-G->nt)) {
		G->cons.add(IloRange(env, -inf, -bl, get_str("flowcaplb(%d)", v).c_str()));
		G->cons[c++].setExpr(expr);
	}

	if (bu < inf) {
		G->cons.add(IloRange(env, -inf, bu, get_str("flowcapub(%d)", v).c_str()));
		G->cons[c++].setExpr(-1*expr);
		expr.end();
	}
	if (G->ns < v && v <= G->nv-G->nt) {
		if (G->mdl == 4) {
			for (int s = 1; s <= G->ns; s++) {
				if (G->s_p(s,v)>=0) {
					IloExpr expr(env);
					for (int j = 0; j < G->A[v].size(); j++)
						expr += G->A[v][j].fx[s];
					int vp = G->is_arc(s,v);
					if (vp>=0)
						expr -= G->A[s][vp].fx[0];
					for (int j = G->ns+1; j <= G->nv-G->nt; j++) {
						int vp = G->is_arc(j,v);
						if (G->is_path(s,j) && vp>=0)
							expr -= G->A[j][vp].fx[s];
					}
					G->cons.add(IloRange(env, 0.0, 0.0, get_str("flowpathblnc(%d,%d)", s, v).c_str()));
					G->cons[c++].setExpr(expr);
					expr.end();
				}
			}
		}
		if (G->mdl == 5) {
			IloExpr expr(env);
			for (int j = 1; j <= G->nv; j++) {
				int vp = G->is_arc(j,v);
				if (vp>=0)
					expr += G->A[j][vp].fx[0];
			}
			for (int j = 0; j < G->A[v].size(); j++)
				expr -= G->A[v][j].fx[0];
			G->cons.add(IloRange(env, 0.0, 0.0, get_str("flowmassblnc(%d)", v).c_str()));
			G->cons[c++].setExpr(expr);
			expr.end();
		}
	}
}

/* define the quality requirement constraints */
void vertex::def_qualrequrmntconsts() {
	if (v <= G->ns) return;
	IloEnv env = G->get_env();
	int c = G->cons.getSize();
	if (1 <= G->mdl && G->mdl <= 4) {
		for (int k = 1; k <= G->nk; k++) {
			if (v > G->nv-G->nt) {
				IloExpr expr(env);
				for (int s = 1; s <= G->ns; s++) {
					for (int j = G->ns+1; j <= G->nv-G->nt; j++) {
						int vp = G->is_arc(j,v);
						if (G->s_p(s,j)>=0 && vp>=0)
							expr += (G->N[s].q[k]-q[k])*G->A[j][vp].fx[s];
					}
					int vp = G->is_arc(s,v);
					if (vp>=0)
						expr += (G->N[s].q[k]-q[k])*G->A[s][vp].fx[0];
				}
				if (0 < fabs(q[k]) && fabs(q[k]) < inf) {
					G->cons.add(IloRange(env, -inf, 0.0, get_str("qualub(%d,%d)", v, k).c_str()));
					G->cons[c++].setExpr(expr);
					expr.end();
				}
			}
		}
		if (G->ns < v && v <= G->nv-G->nt) {
			if (G->mdl == 1 || G->mdl == 3 || G->mdl == 4) {
				IloExpr expr(env);
				for (int s = 1; s <= G->ns; s++) {
					if (G->s_p(s,v)>=0)
						expr += wy[s-1];
				}
				G->cons.add(IloRange(env, 1.0, 1.0, get_str("spropblnc(%d)", v).c_str()));
				G->cons[c++].setExpr(expr);
				expr.end();
			}
			if (G->mdl == 2 || G->mdl == 3) {
				IloExpr expr(env);
				for (int t = G->nv-G->nt+1; t <= G->nv; t++) {
					if (G->is_path(v,t))
						expr += wy[t-1];
				}
				G->cons.add(IloRange(env, 1.0, 1.0, get_str("tpropblnc(%d)", v).c_str()));
				G->cons[c++].setExpr(expr);
				expr.end();
			}
		}
	}
	if (G->mdl == 5) {
		if (G->ns < v && v <= G->nv-G->nt) {
			for (int k = 1; k <= G->nk; k++) {
				IloExpr expr(env);
				for (int s = 1; s <= G->ns; s++) {
					int vp = G->is_arc(s,v);
					if (vp>=0)
						expr += G->N[s].q[k]*G->A[s][vp].fx[0];
				}
				for (int l = G->ns+1; l <= G->nv-G->nt; l++) {
					int vp = G->is_arc(l,v);
					if (vp>=0)
						expr += G->A[l][vp].fx[k];
				}
				for (int j = 0; j < G->A[v].size(); j++) {
					expr -= G->A[v][j].fx[k];
				}
				G->cons.add(IloRange(env, 0.0, 0.0, get_str("qualblnc(%d,%d)", v, k).c_str()));
				G->cons[c++].setExpr(expr);
				expr.end();
			}
		}
		if (v > G->nv-G->nt) {
			for (int k = 1; k <= G->nk; k++) {
				IloExpr expr(env);
				for (int s = 1; s <= G->ns; s++) {
					int vp = G->is_arc(s,v);
					if (vp>=0)
						expr += (G->N[s].q[k]-q[k])*G->A[s][vp].fx[0];
				}
				for (int l = G->ns+1; l <= G->nv-G->nt; l++) {
					int vp = G->is_arc(l,v);
					if (vp>=0)
						expr += G->A[l][vp].fx[k];
				}
				for (int j = G->ns+1; j <= G->nv-G->nt; j++) {
					int vp = G->is_arc(j,v);
					if (vp>=0)
						expr -= q[k]*G->A[j][vp].fx[0];
				}
				if (fabs(q[k]) < inf) {
					G->cons.add(IloRange(env, -inf, 0.0, get_str("qualub(%d,%d)", v, k).c_str()));
					G->cons[c++].setExpr(expr);
					expr.end();
				}
			}
		}
	}
}

/* add the the RLT for the proportion balances */
void arc::def_rltpropor() {
	if (G->mdl == 5 || !(G->ns < i && i <= G->nv-G->nt)) return;
	IloEnv env = G->get_env();
	int c = G->cons.getSize();
	if (G->mdl == 1 || G->mdl == 3 || G->mdl == 4) {
		IloExpr expr(env);
		for (int s = 1; s <= G->ns; s++) {
			if (G->s_p(s,i)>=0)
				expr += fx[s];
		}
		expr -= fx[0];
		G->cons.add(IloRange(env, 0.0, 0.0, ("srlt1"+name).c_str()));
		G->cons[c++].setExpr(expr);
		expr.end();

	}
	if (G->N[i].bu < inf && (G->mdl == 2 || G->mdl == 3)) {
		IloExpr expr(env);
		if (G->mdl == 2) {
			for (int s = 1; s <= G->ns; s++) {
				if (G->is_arc(s,i)>=0)
					expr += fx[s];
			}
		}
		else
			expr += fx[0];
		expr -= G->N[i].bu*G->N[i].wy[j-1];
		G->cons.add(IloRange(env, -inf, 0.0, ("trlt2"+name).c_str()));
		G->cons[c++].setExpr(expr);
		expr.end();
	}
}

/* add the the RLT for the flow capacity */
void vertex::def_rltflowcap() {
	if (G->mdl == 5 || !(G->ns < v && v <= G->nv-G->nt)) return;
	IloEnv env = G->get_env();
	int c = G->cons.getSize();
	if (G->mdl == 2 || G->mdl == 3) {
		for (int s = 1; s <= G->ns; s++) {
			int vp = G->is_arc(s,v);
			if (vp>=0) {
				IloExpr expr(env);
				for (int j = 0; j < G->A[v].size(); j++)
					expr += G->A[v][j].fx[s];
				expr -=G->A[s][vp].fx[0];
				G->cons.add(IloRange(env, 0.0, 0.0, get_str("trlt1(%d,%d)", s, v).c_str()));
				G->cons[c++].setExpr(expr);
				expr.end();
			}
		}
	}
	if (bu < inf && (G->mdl == 1 || G->mdl == 3 || G->mdl == 4)) {
		for (int s = 1; s <= G->ns; s++) {
			int vp = G->is_arc(s,v);
			if (vp>=0) {
				IloExpr expr(env);
				if (G->mdl == 1 || G->mdl == 4) {
					for (int j = 0; j < G->A[v].size(); j++)
						expr += G->A[v][j].fx[s];
				}
				else
					expr += G->A[s][vp].fx[0];
				expr -= bu*wy[s-1];
				G->cons.add(IloRange(env, -inf, 0.0, get_str("srlt2(%d,%d)", s, v).c_str()));
				G->cons[c++].setExpr(expr);
				expr.end();
			}
		}
	}
}

/* define the McCormick envelopes */
void arc::def_McCormickenvlp() {
	if (!(G->ns < i && i <= G->nv-G->nt)) return;
	int ske = G->ns;
	if (G->mdl == 5) ske = G->nk;
	IloEnv env = G->get_env();
	int c = G->cons.getSize();
	if (G->mdl != 2) {
		for (int s = 1; s <= ske; s++) {
			if (G->mdl == 5 || G->s_p(s,i)>=0) {
				double wyl, wyu, fxl, fxu;
				wyl = G->N[i].wy[s-1].getLB();
				wyu = G->N[i].wy[s-1].getUB();
				fxl = fx[0].getLB();
				fxu = fx[0].getUB();
				G->cons.add(IloRange(env, -inf, wyl*fxl, get_str("svexlb(%d,%d,%d)", s, i, j).c_str()));
				G->cons[c].setLinearCoef(fx[s], -1.0);
				G->cons[c].setLinearCoef(fx[0], wyl);
				G->cons[c++].setLinearCoef(G->N[i].wy[s-1], fxl);
				//==========================================
				G->cons.add(IloRange(env, -inf, wyu*fxu, get_str("svexub(%d,%d,%d)", s, i, j).c_str()));
				G->cons[c].setLinearCoef(fx[s], -1.0);
				G->cons[c].setLinearCoef(fx[0], wyu);
				G->cons[c++].setLinearCoef(G->N[i].wy[s-1], fxu);
				//==========================================
				G->cons.add(IloRange(env, -inf, -wyl*fxu, get_str("scavlb(%d,%d,%d)", s, i, j).c_str()));
				G->cons[c].setLinearCoef(fx[s], 1.0);
				G->cons[c].setLinearCoef(fx[0], -wyl);
				G->cons[c++].setLinearCoef(G->N[i].wy[s-1], -fxu);
				//==========================================
				G->cons.add(IloRange(env, -inf, -wyu*fxl, get_str("scavub(%d,%d,%d)", s, i, j).c_str()));
				G->cons[c].setLinearCoef(fx[s], 1.0);
				G->cons[c].setLinearCoef(fx[0], -wyu);
				G->cons[c++].setLinearCoef(G->N[i].wy[s-1], -fxl);
			}
		}
	}
	if (G->mdl == 2 || G->mdl == 3 ) {
		for (int s = 1; s <= G->ns; s++) {
			int v = G->is_arc(s,i);
			if (v>=0) {
				double wyl, wyu, fxl, fxu;
				wyl = G->N[i].wy[j-1].getLB();
				wyu = G->N[i].wy[j-1].getUB();
				fxl = G->A[s][v].fx[0].getLB();
				fxu = G->A[s][v].fx[0].getUB();
				G->cons.add(IloRange(env, -inf, fxl*wyl, get_str("tvexlb(%d,%d,%d)", s, i, j).c_str()));
				G->cons[c].setLinearCoef(fx[s], -1.0);
				G->cons[c].setLinearCoef(G->N[i].wy[j-1], fxl);
				G->cons[c++].setLinearCoef(G->A[s][v].fx[0], wyl);
				//==========================================
				G->cons.add(IloRange(env, -inf, fxu*wyu, get_str("tvexub(%d,%d,%d)", s, i, j).c_str()));
				G->cons[c].setLinearCoef(fx[s], -1.0);
				G->cons[c].setLinearCoef(G->N[i].wy[j-1], fxu);
				G->cons[c++].setLinearCoef(G->A[s][v].fx[0], wyu);
				//==========================================
				G->cons.add(IloRange(env, -inf, -fxl*wyu, get_str("tcavlb(%d,%d,%d)", s, i, j).c_str()));
				G->cons[c].setLinearCoef(fx[s], 1.0);
				G->cons[c].setLinearCoef(G->N[i].wy[j-1], -fxl);
				G->cons[c++].setLinearCoef(G->A[s][v].fx[0], -wyu);
				//==========================================
				G->cons.add(IloRange(env, -inf, -fxu*wyl, get_str("tcavub(%d,%d,%d)", s, i, j).c_str()));
				G->cons[c].setLinearCoef(fx[s], 1.0);
				G->cons[c].setLinearCoef(G->N[i].wy[j-1], -fxu);
				G->cons[c++].setLinearCoef(G->A[s][v].fx[0], -wyl);
			}
		}
	}
}

/* define the McCormick envelopes */
int arc::set_McCormickenvlp(int c) {
	if (i <= G->ns) return c;
	int ske = G->ns;
	if (G->mdl == 5) ske = G->nk;
	if (G->mdl != 2) {
		for (int s = 1; s <= ske; s++) {
			if (G->mdl == 5 || G->s_p(s,i)>=0) {
				double wyl, wyu, fxl, fxu;
				wyl = G->N[i].wy[s-1].getLB();
				wyu = G->N[i].wy[s-1].getUB();
				fxl = fx[0].getLB();
				fxu = fx[0].getUB();
				G->cons[c].setBounds(-inf, wyl*fxl);
				G->cons[c].setLinearCoef(fx[0], wyl);
				G->cons[c++].setLinearCoef(G->N[i].wy[s-1], fxl);
				//==========================================
				G->cons[c].setBounds(-inf, wyu*fxu);
				G->cons[c].setLinearCoef(fx[0], wyu);
				G->cons[c++].setLinearCoef(G->N[i].wy[s-1], fxu);
				//==========================================
				G->cons[c].setBounds(-inf, -wyl*fxu);
				G->cons[c].setLinearCoef(fx[0], -wyl);
				G->cons[c++].setLinearCoef(G->N[i].wy[s-1], -fxu);
				//==========================================
				G->cons[c].setBounds(-inf, -wyu*fxl);
				G->cons[c].setLinearCoef(fx[0], -wyu);
				G->cons[c++].setLinearCoef(G->N[i].wy[s-1], -fxl);
			}
		}
	}
	if (G->mdl == 2 || G->mdl == 3 ) {
		for (int s = 1; s <= G->ns; s++) {
			int v = G->is_arc(s,i);
			if (v>=0) {
				double wyl, wyu, fxl, fxu;
				wyl = G->N[i].wy[j-1].getLB();
				wyu = G->N[i].wy[j-1].getUB();
				fxl = G->A[s][v].fx[0].getLB();
				fxu = G->A[s][v].fx[0].getUB();
				G->cons[c].setBounds(-inf, fxl*wyl);
				G->cons[c].setLinearCoef(G->N[i].wy[j-1], fxl);
				G->cons[c++].setLinearCoef(G->A[s][v].fx[0], wyl);
				//==========================================
				G->cons[c].setBounds(-inf, fxu*wyu);
				G->cons[c].setLinearCoef(G->N[i].wy[j-1], fxu);
				G->cons[c++].setLinearCoef(G->A[s][v].fx[0], wyu);
				//==========================================
				G->cons[c].setBounds(-inf, -fxl*wyu);
				G->cons[c].setLinearCoef(G->N[i].wy[j-1], -fxl);
				G->cons[c++].setLinearCoef(G->A[s][v].fx[0], -wyu);
				//==========================================
				G->cons[c].setBounds(-inf, -fxu*wyl);
				G->cons[c].setLinearCoef(G->N[i].wy[j-1], -fxu);
				G->cons[c++].setLinearCoef(G->A[s][v].fx[0], -wyl);
			}
		}
	}
	return c;
}

/* compute for-loop s in the discretization of proportional variables */
int graph::discretize_proportions(int n, int k, int j) {
	int c[k+3], q[k+1], x;
	memset(c, 0, (k+3)*sizeof(int));
	memset(q, 0, (k+1)*sizeof(int));
	for (int e = 1; e <= k; e++)
		c[e] = e-1;
	c[k+1] = n;
	c[k+2] = 0;
	int i = k;
	visit:
	for (int s = 0; s <= k; s++) {
		if (s == k)
			q[s] = n-c[s]-1;
		else if (s == 0)
			q[s] = c[s+1];
		else
			q[s] = c[s+1]-c[s]-1;
		Y[j].push_back((double)q[s]/(n-k));
	}
	j++;
	if (i > 0) {
		x = i;
		goto increase_c_i;
	}
	if (c[1] + 1 < c[2]) {
		c[1] += 1;
		goto visit;
	}
	else
		i = 2;
	find_i:
	c[i-1] = i-2;
	x = c[i] + 1;
	if (x == c[i+1]) {
		i++;
		goto find_i;
	}
	if (i > k) return j;
	increase_c_i:
	c[i] = x;
	i--;
	goto visit;
}

/* return the position of variable in the extended network */
int graph::i_p(int ij, int jp, int j) {
	int vj = 0, vi = 0;
	if (j == 0) {
		for (int r = 0; r <= jp; r++) {
			int l = A[ij][r].j;
			if (!(ns < l && l <= nv-nt))
				vi += 1;
			else vi += num_p[l-ns]+1;
			if (jp < vi) {
				vj = r;
				break;
			}
		}
	}
	else {
		for (int r = 0; r < jp; r++) {
			int l = A[ij][r].j;
			if (!(ns < l && l <= nv-nt)) vj += 1;
			else vj += num_p[l-ns]+1;
		}
		int l = A[ij][jp].j;
		if (ns < l && l <= nv-nt) {
			int vi = ns;
			for (int u = 1; u <= l-(ns+1); u++)
				vi += num_p[u]+1;
			vj += j-vi-1;
		}
	}
	return vj;
}

/* return the parent of j */
int graph::i_j(int j) {
	if (j <= ns)
		return j;
	else {
		int vi = ns;
		for (int i = ns+1; i <= nv-nt; i++) {
			vi += num_p[i-ns]+1;
			if (j <= vi) return i;
		}
		return j-(vi+nt-nv);
	}
}

/* work the same as is_path(s, i) */
int graph::s_p(int s, int i) {
	int ret = -1;
	for (int j = 0; j < p_indeg[i-ns].size(); j++) {
		if (p_indeg[i-ns][j] == s) {
			ret = j;
			break;
		}
	}
	return ret;
}

/* solve the discrete model */
double graph::discrete_model(int n, int flg, double tm) {
	double pt;
	int ni = nv-ns-nt, i_new = ni;
	num_p.assign(ni+1, 0);
	for (int i = 1; i <= nv-ns-nt; i++) {
		int k = p_indeg[i].size();
		num_p[i] = n_choose_k(n+k-1, k-1)-1;
		i_new += num_p[i];
	}
	int nv_new = ns+i_new+nt;
	int var_c = i_new;
	for (int s = 1; s <= ns; s++) {
		for (int l = 0; l < A[s].size(); l++) {
			if (A[s][l].j > nv-nt)
				var_c++;
		}
	}
	for (int r = ns+1; r <= nv_new-nt; r++) {
		int i = i_j(r);
		for (int l = 0; l < A[i].size(); l++) {
			int j = A[i][l].j;
			if (!(ns < j && j <= nv-nt)) var_c++;
			else var_c += num_p[j-ns]+1;
		}
	}
	if (var_c > 1400000)
		exit(0);
	clock_t start, end;
	start = clock();
	Y = new vector<double>[i_new];
	int j = 0;
	for (int i = 1; i <= nv-ns-nt; i++) {
		int k = p_indeg[i].size();
		j = discretize_proportions(n+k-1, k-1, j);
	}
	if (flg == 1) {
		end = clock();
		pt = (end-start)/(double)CLOCKS_PER_SEC;
		fprintf(fo, "\nDiscretization time (hh:mm:ss:ms) = %s\n", clocktime(pt).c_str());
	}
	double bestfeasible = inf;
	IloEnv envr;
	try {
		vector <IloNumVar> *f;
		f = new vector<IloNumVar>[nv_new-nt+1];
		for (int r = 1; r <= nv_new-nt; r++) {
			int i = i_j(r);
			for (int l = 0; l < A[i].size(); l++) {
				int j = A[i][l].j;
				double ub = fmin(N[i].bu, N[j].bu);
				if (!(ns < j && j <= nv-nt)) {
					IloNumVar var(envr, 0.0, ub, get_str("f(%d,%d)", r, (j <= ns) ? j:j+i_new-ni).c_str());
					f[r].push_back(var);
				}
				else {
					int vj = 0;
					for (int v = 0; v <= j-(ns+1); v++)
						vj += num_p[v];
					for (int v = vj+j; v <= num_p[j-ns]+vj+j; v++) {
						IloNumVar var(envr, 0.0,ub, get_str("f(%d,%d)", r, v).c_str());
						f[r].push_back(var);
					}
				}
			}
		}
		IloIntVarArray p(envr, i_new, 0, 1);
		for (int j = ns+1; j <= nv_new-nt; j++)
			p[j-(ns+1)].setName(get_str("p(%d)", j).c_str());
		/* TODO: reduce the for-loops work, specially the call to other
		 * functions. And this appears to be in all constraints except
		 * pflowcapub and poolcuts.
		 */
		IloModel model(envr);
		IloRangeArray cons(envr);
		int c = cons.getSize();
		IloExpr obj(envr);
		for (int j = 1; j <= nv_new-nt; j++) {
			int i = i_j(j);
			for (int l = 0; l < f[j].size(); l++) {
				int jp = i_p(i,l,0);
				obj += A[i][jp].c*f[j][l];
			}
		}
		model.add(IloMinimize(envr, obj, "objective"));
		obj.end();
		for (int s = 1; s <= ns; s++) {
			IloExpr expr(envr);
			for (int l = 0; l < f[s].size(); l++)
				expr -= f[s][l];
			if (N[s].bl > 0) {
				cons.add(IloRange(envr, -inf, -N[s].bl, get_str("sflowcaplb(%d)", s).c_str()));
				cons[c++].setExpr(expr);
			}
			if (N[s].bu < inf) {
				cons.add(IloRange(envr, -inf, N[s].bu, get_str("sflowcapub(%d)", s).c_str()));
				cons[c++].setExpr(-1*expr);
				expr.end();
			}
		}
		for (int t = ns+1; t <= nv_new; t++) {
			IloExpr expr(envr);
			int l = i_j(t);
			for (int j = 1; j <= nv_new-nt; j++) {
				int r = i_j(j);
				int tp = is_arc(r,l);
				if (tp>=0)
					expr -= f[j][i_p(r,tp,t)];
			}
			double ub = N[l].bu;
			if (t <= nv_new-nt) {
				ub = 0.0;
				expr += N[l].bu*p[t-(ns+1)];
			}
			if (N[l].bl > 0 && t > nv_new-nt) {
				cons.add(IloRange(envr, -inf, -N[l].bl, get_str("tflowcaplb(%d)", t).c_str()));
				cons[c++].setExpr(expr);
			}
			if (N[l].bu < inf) {
				cons.add(IloRange(envr, -inf, ub, get_str("tpflowcapub(%d)", t).c_str()));
				cons[c++].setExpr(-1*expr);
				expr.end();
			}
		}
		for (int s = 1; s <= ns; s++) {
			for (int j = ns+1; j <= nv_new-nt; j++) {
				int i = i_j(j);
				int sp = s_p(s,i);
				if (sp>=0) {
					IloExpr expr(envr);
					for (int l = 0; l < f[j].size(); l++)
						expr += Y[j-(ns+1)][sp]*f[j][l];
					int jp = is_arc(s,i);
					if (jp>=0)
						expr -= f[s][i_p(s,jp,j)];
					for (int l = ns+1; l <= nv_new-nt; l++) {
						int r = i_j(l);
						int ip = s_p(s,r);
						int jp = is_arc(r,i);
						if (ip>=0 && jp>=0)
							expr -= Y[l-(ns+1)][ip]*f[l][i_p(r,jp,j)];
					}
					cons.add(IloRange(envr, 0.0, 0.0, get_str("flowpathblnc(%d,%d)", s, j).c_str()));
					cons[c++].setExpr(expr);
					expr.end();
				}
			}
		}
		for (int k = 1; k <= nk; k++) {
			for (int t = nv_new-nt+1; t <= nv_new; t++) {
				IloExpr expr(envr);
				int l = i_j(t);
				for (int j = 1; j <= nv_new-nt; j++) {
					int r = i_j(j);
					int tp = is_arc(r,l);
					double w_jk = 0.0;
					if (tp>=0) {
						if (r <= ns) w_jk = N[r].q[k];
						else {
							for (int s = 1; s <= ns; s++) {
								int sp = s_p(s,r);
								if (sp >= 0)
									w_jk += N[s].q[k]*Y[j-(ns+1)][sp];
							}
						}
						expr += (w_jk-N[l].q[k])*f[j][i_p(r,tp,t)];
					}
				}
				if (0 < fabs(N[l].q[k]) && fabs(N[l].q[k]) < inf) {
					cons.add(IloRange(envr, -inf, 0.0, get_str("qualub(%d,%d)", k, t).c_str()));
					cons[c++].setExpr(expr);
					expr.end();
				}
			}
		}
		delete[] f, Y;
		int vj = 0;
		for (int i = ns+1; i <= nv-nt; i++) {
			IloExpr expr(envr);
			vj += num_p[i-(ns+1)];
			for (int j = vj+i; j <= num_p[i-ns]+vj+i; j++)
				expr += p[j-(ns+1)];
			cons.add(IloRange(envr, 1.0, 1.0, get_str("poolcut(%d)", i).c_str()));
			cons[c++].setExpr(expr);
			expr.end();
		}
		p.end();
		model.add(cons);
		cons.end();
		IloCplex cplex(envr);
		if (fo != stdout || flg == 0) {
			cplex.setOut(envr.getNullStream());
			cplex.setWarning(envr.getNullStream());
			cplex.setError(envr.getNullStream());
		}
		cplex.extract(model);
		if (flg == 1) {
			vars = cplex.getNcols()-1;
			lconsts = cplex.getNrows();
			intvars = cplex.getNbinVars();
			name = get_str("%s%d", name.c_str(), n).c_str();
			fprintf(fo, "\n");
		}
		//cplex.exportModel(("lps/"+name+"_dsc.lp").c_str());
		cplex.setParam(IloCplex::TiLim, tm);
		cplex.setParam(IloCplex::EpGap, RGAP);
		cplex.setParam(IloCplex::WorkMem, 1500);

		if (cplex.solve()) {
			bestfeasible = cplex.getObjValue();
			if (flg == 1) {
				solstat = 1;
				gpp_ub = cplex.getObjValue();
				gpp_lb = cplex.getBestObjValue();
				ttime = cplex.getTime();
				nnode = cplex.getNnodes();
			}
		}
	}
	catch (IloException& ex) {
		cerr << "Error: " << ex << endl;
	}
	envr.end();
	return bestfeasible;
}

class term{
public:
	int t;
	double obj;
	term(int trm) {
		t = trm;
		obj = inf;
	};
	friend bool operator<(const term &t1, const term &t2) {
		return t1.obj < t2.obj;
	};
};

/* heuristic algorithm to find a feasible solution to the pooling problem */
void graph::heuristic_alg() {
	int p = 1, tau_p;
	list<int> Tp;
	list<term> T_diff_Tp;
	for (int t = nv-nt+1; t <= nv; t++) {
		term trm(t);
		T_diff_Tp.push_back(trm);
	}
	double z_p, totalTime = (double)TIME, time_p = totalTime/(double)T_diff_Tp.size();
	struct timeval end;
	struct timeval start;
	gettimeofday(&start, NULL);
	gpp_ub = 0.0;
	fprintf(fo, "\nIteration %2d:\n", p);
	list<term>::iterator it_p;
	for (list<term>::iterator it = T_diff_Tp.begin(); it != T_diff_Tp.end(); it++) {
		vector <double> *f = one_terminal_mdl(it->t, Tp, time_p);
		fprintf(fo, "  Term %2d, STAT %2.0f, Obj = %g\n", it->t, f[0][1], f[0][0]);
		if (f[0][1] <= 2) {
			it->obj = f[0][0];
			if (f[0][0] < gpp_ub) {
				gpp_ub = f[0][0];
				it_p = it;
				for (int i = 1; i <= nv-nt; i++) {
					for (int j = 0; j < A[i].size(); j++)
						A[i][j].F = f[i][j];
				}
			}
		}
		delete[] f;
	}
	fprintf(fo, "                                 SOL = %g\n", gpp_ub);
	Tp.push_back(it_p->t);
	T_diff_Tp.erase(it_p);
	T_diff_Tp.sort();
	do {
		fprintf(fo, "Iteration %2d:\n", p+1);
		tau_p = T_diff_Tp.front().t;
		gettimeofday(&end, NULL);
		totalTime = (double)TIME-((end.tv_sec-start.tv_sec)+
				(end.tv_usec-start.tv_usec)/(double)CLOCKS_PER_SEC);
		if (totalTime < 0.0) break;
		time_p = totalTime/(double)T_diff_Tp.size();
		vector <double> *f = one_terminal_mdl(tau_p, Tp, time_p);
		fprintf(fo, "  Term %2d, STAT %2.0f, Obj = %g\n", tau_p, f[0][1], f[0][0]);
		if (f[0][1] <= 2 && f[0][0] < 0.0) {
			//printFp();
			gpp_ub = 0.0;
			for (int i = 1; i <= nv-nt; i++) {
				for (int j = 0; j < A[i].size(); j++) {
					A[i][j].F += f[i][j];
					gpp_ub += A[i][j].c*A[i][j].F;
				}
			}
			fprintf(fo, "                                 SOL = %g\n", gpp_ub);
			Tp.push_back(tau_p);
		}
		T_diff_Tp.pop_front();
		p = p+1;
		delete[] f;
	} while (p != nt && ttime < TIME);
	gettimeofday(&end, NULL);
	ttime = (end.tv_sec-start.tv_sec)+(end.tv_usec-start.tv_usec)/(double)CLOCKS_PER_SEC;
	solstat = 1; nnode = p;
}

void graph::printFp() {
	string str;
	fprintf(fo, "parameter Fp(i,j)");
	fprintf(fo, "\n/ ");
	int l = 2;
	for (int i = 1; i <= nv-nt; i++) {
		for (int j = 0; j < A[i].size(); j++) {
			int jp = A[i][j].j;
			if (A[i][j].F==0.0) continue;
			if (j == A[i].size()-1 && i == nv-nt)
				str = get_str("%d.%d %g /;\n\n", i, jp, A[i][j].F);
			else
				str = get_str("%d.%d %g, ", i, jp, A[i][j].F);
			if (l+str.length() <= 80) {
				fprintf(fo, "%s", str.c_str());
				l += str.length();
			}
			else {
				fprintf(fo, "\n  %s", str.c_str());
				l = 2+str.length();
			}
		}
	}
	fprintf(fo, "\n");
}

/* solving the pooling problem with only one terminal */
vector<double> *graph::one_terminal_mdl(int tau, list<int> Tp, double tm) {
	string str; double val;
	vector<double> *fp = new vector<double>[nv-nt+1];
	fp[0].assign(2, 0.0);
	for (int i = 1; i <= nv-nt; i++)
		fp[i].assign(A[i].size(), 0.0);
	if (Tp.empty()) {
		IloEnv envr;
		try {
			vector <IloNumVar> *f;
			f = new vector<IloNumVar>[nv-nt+1];
			for (int i = 1; i <= nv-nt; i++) {
				for (int jp = 0; jp < A[i].size(); jp++) {
					int j = A[i][jp].j;
					double ub = fmin(N[i].bu, N[j].bu);
					str = get_str("f(%d,%d)", i, j);
					IloNumVar var(envr, 0.0, ub, str.c_str());
					f[i].push_back(var);
				}
			}
			IloModel model(envr);
			IloRangeArray cons(envr);
			IloExpr obj(envr);
			for (int i = 1; i <= nv-nt; i++) {
				for (int j = 0; j < A[i].size(); j++) {
					if (is_path(A[i][j].j,tau))
						obj += A[i][j].c*f[i][j];
				}
			}
			model.add(IloMinimize(envr, obj, "objective"));
			obj.end();
			int c = cons.getSize();
			for (int i = 1; i <= nv; i++) {
				if (is_path(i,tau) || i == tau) {
					IloExpr expr(envr); val = 0;
					if (i <= ns) {
						for (int j = 0; j < A[i].size(); j++) {
							if (is_path(A[i][j].j, tau))
								expr -= f[i][j];
						}
					}
					else {
						for (int j = 1; j <= nv; j++) {
							int ip = is_arc(j,i);
							if (ip>=0)
								expr -= f[j][ip];
						}
					}
					if (N[i].bl > 0 && !(ns < i && i <= nv-nt)) {
						str = get_str("flowcaplb(%d)", i);
						cons.add(IloRange(envr, -inf, -N[i].bl, str.c_str()));
						cons[c++].setExpr(expr);
					}
					if (N[i].bu < inf) {
						str = get_str("flowcapub(%d)", i);
						cons.add(IloRange(envr, -inf, N[i].bu, str.c_str()));
						cons[c++].setExpr(-1*expr);
					}
					expr.end();
				}
			}
			for (int i = ns+1; i <= nv-nt; i++) {
				if (is_path(i,tau)) {
					IloExpr expr(envr);
					for (int j = 1; j <= nv-nt; j++) {
						int ip = is_arc(j,i);
						if (ip>=0)
							expr += f[j][ip];
					}
					for (int j = 0; j < A[i].size(); j++) {
						if (is_path(A[i][j].j,tau))
							expr -= f[i][j];
					}
					str = get_str("flowmassblnc(%d)", i);
					cons.add(IloRange(envr, 0.0, 0.0, str.c_str()));
					cons[c++].setExpr(expr); expr.end();
				}
			}
			for (int k = 1; k <= nk; k++) {
				IloExpr expr(envr);
				for (int s = 1; s <= ns; s++) {
					int ip = is_arc(s,tau);
					if (ip>=0)
						expr += (N[s].q[k]-N[tau].q[k])*f[s][ip];
					for (int jp = 0; jp < A[s].size(); jp++) {
						int j = A[s][jp].j;
						if (j<= nv-nt && is_path(j,tau))
							expr += (N[s].q[k]-N[tau].q[k])*f[s][jp];
					}
				}
				if (fabs(N[tau].q[k])>0 && fabs(N[tau].q[k])<inf) {
					str = get_str("qualub(%d,%d)", tau, k);
					cons.add(IloRange(envr, -inf, 0.0, str.c_str()));
					cons[c++].setExpr(expr); expr.end();
				}
			}
			model.add(cons);
			cons.end();
			IloCplex cplex(envr);
			cplex.extract(model);
			cplex.setOut(envr.getNullStream());
			cplex.setWarning(envr.getNullStream());
			cplex.setError(envr.getNullStream());
			cplex.setParam(IloCplex::EpGap, 1.0e-3);
			//cplex.exportModel(get_str("lps/%s%d.lp", name.c_str(),tau).c_str());
			if (cplex.solve()) {
				fp[0][0] = cplex.getObjValue();
				fp[0][1] = 1;
				for (int i = 1; i <= nv-nt; i++) {
					for (int jp = 0; jp < A[i].size(); jp++) {
						if (is_path(A[i][jp].j,tau))
							fp[i][jp] = cplex.getValue(f[i][jp]);
					}
				}
			}
			delete[] f;
		}
		catch (IloException& ex) {
			cerr << "Error: " << ex << endl;
		}
		envr.end();
	}
	else {
		string str;
		FILE *ft = fopen("terminal.gms", "w");
		fprintf(ft, "scalar ttime; ttime = %g;\n", tm);
		fprintf(ft, "parameter Fp(i,j)");
		fprintf(ft, "\n/ ");
		int l = 2;
		for (int i = 1; i <= nv-nt; i++) {
			for (int j = 0; j < A[i].size(); j++) {
				int jp = A[i][j].j;

				if (j == A[i].size()-1 && i == nv-nt)
					str = get_str("%d.%d %g /;\n\n", i, jp, A[i][j].F);
				else
					str = get_str("%d.%d %g, ", i, jp, A[i][j].F);
				if (l+str.length() <= 80) {
					fprintf(ft, "%s", str.c_str());
					l += str.length();
				}
				else {
					fprintf(ft, "\n  %s", str.c_str());
					l = 2+str.length();
				}
			}
		}

		int ErrNr, VarNr, NrRecs, Dim, VarTyp;
		string sysdir = getenv("GAMS_HOME");
		string gams = sysdir+"/gams", msg, VarName;

		GDX gdx(sysdir, msg);
		/*gdx.OpenWrite("Fdata.gdx", "", ErrNr);

		gdx.DataWriteStrStart("ttime", "", 1, GMS_DT_PAR, 0);
		Indx[0] = get_str("%d", 1);
		Values[GMS_VAL_LEVEL] = tm;
		gdx.DataWriteStr(Indx, Values);
		gdx.DataWriteDone();

		gdx.DataWriteStrStart("Fp", "", 2, GMS_DT_PAR, 0);
		for (int i = 1; i <= nv-nt; i++) {
			for (int j = 0; j < A[i].size(); j++) {
				int jp = A[i][j].j;
				if (A[i][j].F == 0.0) continue;
				Indx[0] = get_str("%d", i);
				Indx[1] = get_str("%d", jp);
				Values[GMS_VAL_LEVEL] = A[i][j].F;
				gdx.DataWriteStr(Indx, Values);
			}
		}
		gdx.DataWriteDone();
		//gdx.Close();

		FILE *ft = fopen("terminal.gms", "w");*/
		fprintf(ft, "set Tp(t) /");
		for (list<int>::iterator p = Tp.begin(); p != Tp.end(); p++) {
			if (p == Tp.begin())
				fprintf(ft, "%d", *p);
			else
				fprintf(ft, ", %d", *p);
		}
		fprintf(ft, "/;\n");
		fprintf(ft, "set tau(t) /%d/;", tau);
		fclose(ft);

		ErrNr = system((gams+" gms/"+name+".gms o=gms/lst/"+name+".lst lo=0 suppress=1").c_str());

		//GDX gdx(sysdir, msg);
		gdx.OpenRead("results.gdx", ErrNr);

		VarName = "z_p";
		gdx.FindSymbol(VarName, VarNr);
		gdx.SymbolInfo(VarNr,VarName,Dim,VarTyp);
		gdx.DataReadStrStart(VarNr,NrRecs);
		gdx.DataReadStr(Indx,Values,NrRecs);
		fp[0][0] = Values[GMS_VAL_LEVEL];
		gdx.DataReadDone();

		VarName = "mdlstat";
		gdx.FindSymbol(VarName, VarNr);
		gdx.SymbolInfo(VarNr,VarName,Dim,VarTyp);
		gdx.DataReadStrStart(VarNr,NrRecs);
		gdx.DataReadStr(Indx,Values,NrRecs);
		fp[0][1] = Values[GMS_VAL_LEVEL];
		gdx.DataReadDone();

		VarName = "f_p";
		gdx.FindSymbol(VarName, VarNr);
		gdx.SymbolInfo(VarNr,VarName,Dim,VarTyp);
		gdx.DataReadStrStart(VarNr,NrRecs);
		while (gdx.DataReadStr(Indx,Values,NrRecs)) {
			if (Values[GMS_VAL_LEVEL] == 0.0) continue;
			int i = strtol(Indx[0].c_str(), NULL, 10);
			int j = strtol(Indx[1].c_str(), NULL, 10);
			int jp = is_arc(i,j);
			fp[i][jp] = Values[GMS_VAL_LEVEL];
		}
		gdx.DataReadDone();

    VarName = "psiz";
		gdx.FindSymbol(VarName, VarNr);
		gdx.SymbolInfo(VarNr,VarName,Dim,VarTyp);
		gdx.DataReadStrStart(VarNr,NrRecs);
		while (gdx.DataReadStr(Indx,Values,NrRecs)) {
			int i = strtol(Indx[0].c_str(), NULL, 10);
			int val = (int)Values[GMS_VAL_LEVEL];
			if (i==1) vars = max(vars, val);
			if (i==2) nlterms = max(nlterms, val);
			if (i==3) nlconsts = max(nlconsts, val);
			if (i==4) lconsts = max(lconsts, val);
		}
		gdx.DataReadDone();
		gdx.Close();
	}
	return fp;
}

/* eof */
