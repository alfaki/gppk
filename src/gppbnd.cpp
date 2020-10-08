/* gppbnd.cpp */

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

double upper_bnd;

/* sort the nodes by its objective value */
bool operator<(const node &n1, const node &n2) {
	return n1.obj < n2.obj;
}

/* prune a node by bound */
bool prune_by_bound(const node &n) {
	return (upper_bnd-n.obj < AGAP);
}

/* set the variable bounds */
void graph::set_bounds(node *curr_n, int itr) {
	int sks = 1, ske = ns;
	if (mdl == 2) {
		sks = nv-nt+1; ske = nv;
	}
	if (mdl == 3) ske = nv;
	if (mdl == 5) ske = nk;
	for (int s = sks; s <= ske; s++) {
		for (int i = 1; i <= nv-ns-nt; i++)
			N[ns+i].wy[s-1].setBounds(curr_n->wyl[s][i], curr_n->wyu[s][i]);
	}
	if (mdl != 2 ) {
		for (int i = ns+1; i <= nv; i++) {
			for (int j = 0; j < A[i].size(); j++)
				A[i][j].fx[0].setBounds(curr_n->fl[i][j], curr_n->fu[i][j]);
		}
	}
	if (mdl != 1) {
		for (int s = 1; s <= ns; s++) {
			for (int j = 0; j < A[s].size(); j++) {
				if (A[s][j].j > nv-nt) continue;
				A[s][j].fx[0].setBounds(curr_n->fl[s][j],curr_n->fu[s][j]);
			}
		}
	}
	int c = nlp_consts;
	for (int i = 1; i <= nv; i++) {
		for (int j = 0; j < A[i].size(); j++)
			c = A[i][j].set_McCormickenvlp(c);
	}
	model.add(cons);
	cplex.extract(model);
	/*if (itr < 5)
		cplex.exportModel(get_str("lps/%s_%d.lp", name, itr).c_str());*/
}

/* delete the data in node n */
void delete_node(node *n) {
    delete[] n->fl;
    delete[] n->fu;
    delete[] n->wyl;
    delete[] n->wyu;
}

/* check the the problem feasibility */
int graph::check_feasibility() {
	int sks = 1, ske = ns;
	if (mdl == 2) {
		sks = nv-nt+1; ske = nv;
	}
	if (mdl == 3) ske = nv;
	if (mdl == 5) ske = nk;
	node root;
	root.obj = -inf;
	root.wyl = new vector<double>[ske+1];
	root.wyu = new vector<double>[ske+1];
	for (int s = sks; s <= ske; s++) {
		root.wyl[s].assign(nv-ns-nt+1, 0.0);
		root.wyu[s].assign(nv-ns-nt+1, 0.0);
		for (int i = 1; i <= nv-ns-nt ; i++) {
			root.wyl[s][i] = N[ns+i].wy[s-1].getLB();
			root.wyu[s][i] = N[ns+i].wy[s-1].getUB();
		}
	}
	root.fl = new vector<double>[nv-nt+1];
	root.fu = new vector<double>[nv-nt+1];
	for (int i = 1; i <= nv-nt; i++) {
		root.fl[i].assign(A[i].size(), 0.0);
		root.fu[i].assign(A[i].size(), 0.0);
		for (int j = 0; j < A[i].size(); j++) {
			root.fu[i][j] = A[i][j].fx[0].getLB();
			root.fu[i][j] = A[i][j].fx[0].getUB();
		}
	}
	intvars = 0;
	problem_info();
	fprintf(fo, "\n Iteration   Open nodes   Best objective    Best feasible        Gap\n");
	if (cplex.solve()) {
		gpp_lb = cplex.getObjValue();
		if (nlp_feasibility()) {
			gpp_ub = cplex.getObjValue();
			fprintf(fo, "%10d %12d %+16.6f %+16.6f %10.6f%%\n", 1, 0, gpp_lb, gpp_ub, 0.0);
            delete_node(&root);
			return 1;
		}
		/*double val = discrete_model(1, 0, 60);
		if (val < gpp_ub)
			gpp_ub = val;*/
		/* TODO: implement the OBBT preprocessing */
		double rgap = fabs(gpp_lb-gpp_ub)/fabs(gpp_lb);
		if (rgap < RGAP) {
			fprintf(fo, "%10d %12d %+16.6f %+16.6f %10.6f%%\n", 1, 0, gpp_lb, gpp_ub, rgap);
            delete_node(&root);
			return 1;
		}
		fprintf(fo, "%10d %12d %+16.6f %+16.6f %10.6f%%\n", 1, 4, gpp_lb, gpp_ub, rgap);
		/* "splitting the rectangle" to create child nodes */
		add_child_nodes(root, gpp_lb);
	}
	else {
		fprintf(fo, "Problem %s is infeasible\n", name.c_str());
        delete_node(&root);
		return 0;
	}
	delete_node(&root);
	return 2;
}

/* check the solution of the current node if it feasible */
bool graph::nlp_feasibility() {
	int ske = ns;
	if (mdl == 5) ske = nk;
	bool feas;
	double mfeas_s = 0.0, mfeas_t = 0.0;
	int ss = 0, is = 0, js = 0, st = 0, it = 0, jt = 0;
	if (mdl != 2) {
		for (int s = 1; s <= ske; s++) {
			for (int i = ns+1; i <= nv-nt; i++) {
				for (int j= 0; j < A[i].size(); j++) {
					if (mdl == 5 || is_path(s,i)) {
						double wyl, wyu, wys, fxl, fxu, fxs;
						wys = cplex.getValue(N[i].wy[s-1]);
						wyl = N[i].wy[s-1].getLB();
						wyu = N[i].wy[s-1].getUB();
						fxs = cplex.getValue(A[i][j].fx[0]);
						fxl = A[i][j].fx[0].getLB();
						fxu = A[i][j].fx[0].getUB();
						double vex = wys*fxs-fmax(wyl*fxs+wys*fxl-wyl*fxl, wyu*fxs+wys*fxu-wyu*fxu);
						double cav = fmin(wyl*fxs+wys*fxu-wyl*fxu, wyu*fxs+wys*fxl-wyu*fxl)-wys*fxs;
						if (fmax(vex,cav) > mfeas_s) {
							mfeas_s = fmax(vex,cav);
							ss = s; is = i; js = A[i][j].j;
						}
					}
				}
			}
		}
		if (mdl != 3) {
			if (fabs(mfeas_s) < AGAP)
				feas = true;
			else {
				feas = false;
				source_McCormick = true;
				trngl[0] = ss; trngl[1] = is; trngl[2] = js;
			}
		}
	}
	if (mdl == 2 || mdl == 3) {
		for (int s = 1; s <= ns; s++) {
			for (int j = 0; j < A[s].size(); j++) {
				for (int t= nv-nt+1; t <= nv; t++) {
					int i = A[s][j].j;
					if (is_arc(i,t)>=0) {
						double wyl, wyu, wys, fxl, fxu, fxs;
						wys = cplex.getValue(N[i].wy[t-1]);
						wyl = N[i].wy[t-1].getLB();
						wyu = N[i].wy[t-1].getUB();
						fxs = cplex.getValue(A[s][j].fx[0]);
						fxl = A[s][j].fx[0].getLB();
						fxu = A[s][j].fx[0].getUB();
						double vex = wys*fxs-fmax(wyl*fxs+wys*fxl-wyl*fxl, wyu*fxs+wys*fxu-wyu*fxu);
						double cav = fmin(wyl*fxs+wys*fxu-wyl*fxu, wyu*fxs+wys*fxl-wyu*fxl)-wys*fxs;
						if (fmax(vex,cav) > mfeas_t) {
							mfeas_t = fmax(vex,cav);
							st = s; it = i; jt = t;
						}
					}
				}
			}
		}
		if (mdl == 2) {
			if (fabs(mfeas_t) < AGAP)
				feas = true;
			else {
				feas = false;
				source_McCormick = false;
				trngl[0] = st; trngl[1] = it; trngl[2] = jt;
			}
		}
		if (mdl == 3) {
			if (BRNCH == 1) { /* min-max-infeas */
				if (fabs(mfeas_s) < AGAP || fabs(mfeas_t) < AGAP)
					feas = true;
				else {
					feas = false;
					if (fabs(mfeas_s) < fabs(mfeas_t)) {
						source_McCormick = true;
						trngl[0] = ss; trngl[1] = is; trngl[2] = js;
					}
					else {
						source_McCormick = false;
						trngl[0] = st; trngl[1] = it; trngl[2] = jt;
					}
				}
			}
			if (BRNCH == 2) { /* max-infeas */
				if (fabs(mfeas_s) < AGAP && fabs(mfeas_t) < AGAP)
					feas = true;
				else {
					feas = false;
					if (fabs(mfeas_s) > fabs(mfeas_t)) {
						source_McCormick = true;
						trngl[0] = ss; trngl[1] = is; trngl[2] = js;
					}
					else {
						source_McCormick = false;
						trngl[0] = st; trngl[1] = it; trngl[2] = jt;
					}
				}
			}
		}
	}
	return feas;
}

/* branch and bound */
void graph::branch_and_bound(time_t start) {
	clock_t end;
	ttime = 0.0; solstat = 1;
	int i = 1, si = 0;
	double rgap = fabs(gpp_ub-gpp_lb)/fabs(gpp_lb);
	while (!(tree.empty() || rgap <= RGAP || ttime > TIME)) {
		i = i+1;
		/* Node selection */
		tree.sort();
		node curr_n = tree.front();
		gpp_lb = curr_n.obj;
		/* remove current Node from waiting list */
		tree.pop_front();
		/* create subproblem */
		set_bounds(&curr_n, i);
		/* check subproblem "optimality" */
		if (cplex.solve()) {
			/* check "feasibility" */
			double obj = cplex.getObjValue();
			if (nlp_feasibility()) {
				if (obj < gpp_ub) {
					/* improved the best found Solution */
					gpp_ub = obj;
					/* keep this solution */
					//finalsol = sol;
					si = i;
				}
				/* this Node pruned by bound */
			}
			else {
				/*double val = discrete_model(fo, 1, 0);
				if (val < gpp_ub)
					gpp_ub = val;*/
				/* "splitting the rectangle" to create child nodes */
				add_child_nodes(curr_n, obj);
			}
		}
		delete_node(&curr_n);
		/* remove all waiting Nodes with Node.obj > gpp_ub */
		for(list<node>::iterator n=tree.begin(); n!=tree.end(); ++n) {
			if(gpp_ub - n->obj < AGAP)
				delete_node(&(*n));
		}
		upper_bnd = gpp_ub;
		tree.remove_if(prune_by_bound);
		if (tree.empty())
			gpp_lb = gpp_ub;
		/* this Node is fathomed by infeasibility */
		rgap = fabs(gpp_ub-gpp_lb)/fabs(gpp_lb);
		if (si == i || i%100 == 0 || tree.empty() || rgap <= RGAP)
			fprintf(fo, "%10d %12d %+16.6f %+16.6f %10.6f%%\n", i, (int)tree.size(), gpp_lb, gpp_ub, rgap);
		nnode = i;
		end = clock();
		ttime = (end-start)/(double)CLOCKS_PER_SEC;
	}
	end = clock();
	ttime = (end-start)/(double)CLOCKS_PER_SEC;
}

/* create a child node */
node graph::create_child(node *curr_n, double obj) {
	int sks = 1, ske = ns;
	if (mdl == 2) {
		sks = nv-nt+1; ske = nv;
	}
	if (mdl == 3) ske = nv;
	if (mdl == 5) ske = nk;
	node chld;
	chld.obj = obj;
	chld.wyl = new vector<double>[ske+1];
	chld.wyu = new vector<double>[ske+1];
	for (int s = sks; s <= ske; s++) {
		chld.wyl[s].assign(nv-ns-nt+1, 0.0);
		chld.wyu[s].assign(nv-ns-nt+1, 0.0);
		for (int i = 1; i <= nv-ns-nt; i++) {
			chld.wyl[s][i] = curr_n->wyl[s][i];
			chld.wyu[s][i] = curr_n->wyu[s][i];
		}
	}
	chld.fl = new vector<double>[nv-nt+1];
	chld.fu = new vector<double>[nv-nt+1];
	for (int i = 1; i <= nv-nt; i++) {
		chld.fl[i].assign(A[i].size(), 0.0);
		chld.fu[i].assign(A[i].size(), 0.0);
		for (int j = 0; j < A[i].size(); j++) {
			chld.fl[i][j] = curr_n->fl[i][j];
			chld.fu[i][j] = curr_n->fu[i][j];
		}
	}
	return chld;
}

/* create child nodes for the current node */
void graph::add_child_nodes(node &curr_n, double obj) {
	/* "splitting the rectangle" to create child nodes */
	node child;
	if (source_McCormick) {
		int s = trngl[0], i = trngl[1], j = is_arc(i, trngl[2]);
		/* child 1 */
		child = create_child(&curr_n, obj);
		child.wyu[s][i-ns] = cplex.getValue(N[i].wy[s-1]);
		child.fu[i][j] = cplex.getValue(A[i][j].fx[0]);
		tree.push_back(child);
		/* child 2 */
		child = create_child(&curr_n, obj);
		child.wyl[s][i-ns] = cplex.getValue(N[i].wy[s-1]);
		child.fu[i][j] = cplex.getValue(A[i][j].fx[0]);
		tree.push_back(child);
		/* child 3 */
		child = create_child(&curr_n, obj);
		child.wyl[s][i-ns] = cplex.getValue(N[i].wy[s-1]);
		child.fl[i][j] = cplex.getValue(A[i][j].fx[0]);
		tree.push_back(child);
		/* child 4 */
		child = create_child(&curr_n, obj);
		child.wyu[s][i-ns] = cplex.getValue(N[i].wy[s-1]);
		child.fl[i][j] = cplex.getValue(A[i][j].fx[0]);
		tree.push_back(child);
	}
	else {
		int s = trngl[0], i = trngl[1], i_p = is_arc(s, trngl[1]), t = trngl[2];
		/* child 1 */
		child = create_child(&curr_n, obj);
		child.wyu[t][i-ns] = cplex.getValue(N[i].wy[t-1]);
		child.fu[s][i_p] = cplex.getValue(A[s][i_p].fx[0]);
		tree.push_back(child);
		/* child 2 */
		child = create_child(&curr_n, obj);
		child.wyl[t][i-ns] = cplex.getValue(N[i].wy[t-1]);
		child.fu[s][i_p] = cplex.getValue(A[s][i_p].fx[0]);
		tree.push_back(child);
		/* child 3 */
		child = create_child(&curr_n, obj);
		child.wyl[t][i-ns] = cplex.getValue(N[i].wy[t-1]);
		child.fl[s][i_p] = cplex.getValue(A[s][i_p].fx[0]);
		tree.push_back(child);
		/* child 4 */
		child = create_child(&curr_n, obj);
		child.wyu[t][i-ns] = cplex.getValue(N[i].wy[t-1]);
		child.fl[s][i_p] = cplex.getValue(A[s][i_p].fx[0]);
		tree.push_back(child);
	}
}

/* eof */
