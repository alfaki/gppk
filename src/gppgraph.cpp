/* gppgraph.cpp (reading/writing pooling problem data in DIMACS format ) */

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
 *****************************************************************************/
#include <set>
#include "gppk.h"

/* check whether the arc (i,j) in the network or not */
int graph::is_arc(int i, int j) {
	int isarc = -1;
	for (int it = 0; it < A[i].size(); it++) {
		if (A[i][it].j == j) {
			isarc = it;
			break;
		}
	}
	return isarc;
}

/* check whether there is a path between s and l or not */
bool graph::is_path(int s, int l) {
	/* Breadth First Search (BFS) algorithm */
	bool istherepath = false;
	list<int> openlst;
	list<int> closedlst;
	vector<bool> inclosed(nv);
	openlst.push_back(s);
	while (!openlst.empty()) {
		int u = openlst.front();
		openlst.pop_front();
		if (u == l) {
			istherepath = true;
			break;
		}
		else {
			closedlst.push_back(u);
			inclosed[u] = true;
			for (int j = 0; j < A[u].size(); j++) {
				if (!inclosed[A[u][j].j]) {
					openlst.push_back(A[u][j].j);
					openlst.unique();
				}
			}
		}
	}
	return istherepath;
}

/* common storage area */
struct csa {
	const char *fname;
	/* name of input text file */
	FILE *fp;
	/* stream assigned to input text file */
	int count;
	/* line count */
	int c;
	/* current character */
	char field[255+1];
	/* data field */
	int empty;
	/* warning 'empty line ignored' was printed */
	int nonint;
	/* warning 'non-integer data detected' was printed */
};

/* print error message and terminate processing */
static void error(struct csa *csa, const char *fmt, ...) {
	va_list arg;
	printf("%s:%d: error: ", csa->fname, csa->count);
	va_start(arg, fmt);
	vprintf(fmt, arg);
	va_end(arg);
	printf("\n");
	longjmp(buf, 1);
	/* no return */
}

/* print warning message and continue processing */
static void warning(struct csa *csa, const char *fmt, ...) {
	va_list arg;
	printf("%s:%d: warning: ", csa->fname, csa->count);
	va_start(arg, fmt);
	vprintf(fmt, arg);
	va_end(arg);
	printf("\n");
	return;
}

/* read character from input text file */
static void read_char(struct csa *csa) {
	int c;
	if (csa->c == '\n') csa->count++;
	c = fgetc(csa->fp);
	if (c < 0) {
		if (ferror(csa->fp))
			error(csa, "read error - ");
		else if (csa->c == '\n')
			error(csa, "unexpected end of file");
		else {
			warning(csa, "missing final end of line");
			c = '\n';
		}
	}
	else if (c == '\n');
	else if (isspace(c))
		c = ' ';
	else if (iscntrl(c))
		error(csa, "invalid control character 0x%02X", c);
	csa->c = c;
	return;
}

/* read one-character line designator */
static void read_designator(struct csa *csa) {
	assert(csa->c == '\n');
	read_char(csa);
	for (;;) {
		/* skip preceding white-space characters */
		while (csa->c == ' ')
			read_char(csa);
		if (csa->c == '\n') {
			/* ignore empty line */
			if (!csa->empty) {
				warning(csa, "empty line ignored");
				csa->empty = 1;
			}
			read_char(csa);
		}
		else if (csa->c == 'c') {
			/* skip comment line */
			while (csa->c != '\n')
				read_char(csa);
			read_char(csa);
		}
		else {
			/* hmm... looks like a line designator */
			csa->field[0] = (char)csa->c, csa->field[1] = '\0';
			/* check that it is followed by a white-space character */
			read_char(csa);
			if (!(csa->c == ' ' || csa->c == '\n'))
				error(csa, "line designator missing or invalid");
			break;
		}
	}
	return;
}

/* read data field */
static void read_field(struct csa *csa) {
	int len = 0;
	/* skip preceding white-space characters */
	while (csa->c == ' ')
		read_char(csa);
	/* scan data field */
	if (csa->c == '\n')
		error(csa, "unexpected end of line");
	while (!(csa->c == ' ' || csa->c == '\n')) {
		if (len == sizeof(csa->field)-1)
			error(csa, "data field '%.15s...' too long", csa->field);
		csa->field[len++] = (char)csa->c;
		read_char(csa);
	}
	csa->field[len] = '\0';
	return;
}

/* skip white-space characters until end of line */
static void end_of_line(struct csa *csa) {
	while (csa->c == ' ')
		read_char(csa);
	if (csa->c != '\n')
		error(csa, "too many data fields specified");
	return;
}

/* print a warning if non-integer data are detected */
static void check_int(struct csa *csa, double num) {
	if (!csa->nonint && num != floor(num)) {
		warning(csa, "non-integer data detected");
		csa->nonint = 1;
	}
	return;
}

/******************************************************************************
 *  NAME
 *  read_graph - read pooling problem data in DIMACS format
 *
 *  SYNOPSIS
 *  int read_graph(const char *fname);
 *
 *  DESCRIPTION
 *  The routine read_graph reads pooling problem data in DIMACS
 *  format from a text file.
 *
 *  RETURNS
 *  If the operation was successful, the routine returns zero. Otherwise
 *  it prints an error message and returns non-zero.
 *****************************************************************************/

int graph::read_graph(const char *fname, int mod) {
	mdl = mod; isSpp = true;
	if (1 <= mdl && mdl <= 5) {
		cons = IloRangeArray(env);
		model = IloModel(env);
		cplex = IloCplex(env);
	}
	struct csa _csa, *csa = &_csa;
	int ret = 0;
	double c, bl, bu, q;
	csa->fname = fname;
	csa->fp = NULL;
	csa->count = 0;
	csa->c = '\n';
	csa->field[0] = '\0';
	csa->empty = csa->nonint = 0;
	fprintf(fo, "\nReading data from '%s'...\n", fname);
	csa->fp = fopen(fname, "r");
	if (csa->fp == NULL) {
		printf("Unable to open '%s'\n\n", fname);
		exit(0);
	}
	name = file_name(fname);
	/* read problem line */
	read_designator(csa);
	if (strcmp(csa->field, "p") != 0)
		error(csa, "problem line missing or invalid");
	read_field(csa);
	if (strcmp(csa->field, "min") != 0)
		error(csa, "wrong problem designator; 'min' expected");
	read_field(csa);
	if (!(str2int(csa->field, &na) == 0 && na >= 0))
		error(csa, "number of arcs missing or invalid");
	read_field(csa);
	if (!(str2int(csa->field, &nv) == 0 && nv >= 0))
		error(csa, "number of nodes missing or invalid");
	read_field(csa);
	if (!(str2int(csa->field, &ns) == 0 && ns >= 0 && ns < nv))
		error(csa, "number of source nodes missing or invalid");
	read_field(csa);
	if (!(str2int(csa->field, &nt) == 0 && nt >= 0 && ns+nt < nv))
		error(csa, "number of terminal nodes missing or invalid");
	read_field(csa);
	if (!(str2int(csa->field, &nk) == 0 && nk >= 0))
		error(csa, "number of quality attributes missing or invalid");
	int ni = nv-ns-nt;
	fprintf(fo, "The pooling network '%s' ",        name.c_str());
	fprintf(fo, "has %d source%s, ",      ns, ns == 1 ? "" : "s");
	fprintf(fo, "%d pool%s,\n",           ni, ni == 1 ? "" : "s");
	fprintf(fo, "%d terminal%s and ",     nt, nt == 1 ? "" : "s");
	fprintf(fo, "%d arc%s ",              na, na == 1 ? "" : "s");
	fprintf(fo, "with %d qualit%s.\n", nk, nk == 1 ? "y" : "ies");
	end_of_line(csa);
	N.assign(nv+1, vertex());
	A = new vector<arc>[nv+1];
	/* read node descriptor lines */
	for (;;) {
		bl = 0.0;
		read_designator(csa);
		if (strcmp(csa->field, "n") != 0) break;
		read_field(csa);
		int i;
		if (str2int(csa->field, &i) != 0)
			error(csa, "node number missing or invalid");
		if (!(1 <= i && i <= nv))
			error(csa, "node number %d out of range", i);
		if (!(ns < i && i <= nv-nt)) {
			read_field(csa);
			if (str2num(csa->field, &bl) != 0)
				error(csa, "node bl capacity missing or invalid");
		}
		read_field(csa);
		if (str2num(csa->field, &bu) != 0)
			error(csa, "node ub capacity missing or invalid");
		N[i] = vertex(this, i, bl, bu);
		if (i <= ns || i > nv-nt) {
			(N[i].q).assign(nk+1, 0.0);
			for (int k = 1; k <= nk; k++) {
				read_field(csa);
				if (str2num(csa->field, &q) != 0)
					error(csa, "\nquality missing or invalid");
				N[i].q[k] = q;
			}
		}
		end_of_line(csa);
	}
	/* read arc descriptor lines */
	for (int a = 1; a <= na; a++) {
		if (a > 1) read_designator(csa);
		if (strcmp(csa->field, "a") != 0)
			error(csa, "wrong line designator; 'a' expected");
		read_field(csa);
		int i, j;
		if (str2int(csa->field, &i) != 0)
			error(csa, "starting node number missing or invalid");
		if (!(1 <= i && i <= nv))
			error(csa, "starting node number %d out of range", i);
		read_field(csa);
		if (str2int(csa->field, &j) != 0)
			error(csa, "ending node number missing or invalid");
		if (!(1 <= j && j <= nv))
			error(csa, "ending node number %d out of range", j);
		read_field(csa);
		if (str2num(csa->field, &c) != 0)
			error(csa, "per-unit cost of arc flow missing or invalid");
		if (isSpp && (ns < i && i <= nv-nt) && (ns < j && j <= nv-nt))
			isSpp = false;
		arc edge(this, i, j, c);
		A[i].push_back(edge);
		end_of_line(csa);
	}
	fprintf(fo, "%d lines were read.\n", csa->count);
	if (csa->fp != NULL)
		fclose(csa->fp);
	/**************************************************************************
	FILE *fp = fopen("test.dat", "w");
	fprintf(fp, "[[");
	for (int s = 1; s <= ns; s++) {
		for (int i = ns+1; i <= nv-nt; i++) {
			if (s < ns)
				fprintf(fp,"%g%s",is_arc(s,i)>=0?1.0:0.0, i==nv-nt?"],\n [":",");
			else
				fprintf(fp,"%g%s",is_arc(s,i)>=0?1.0:0.0,i==nv-nt?"]]\n":",");
		}
	}
	fprintf(fp, "[[");
	for (int i = ns+1; i <= nv-nt; i++) {
		for (int t = nv-nt+1; t <= nv; t++) {
			if (i < nv-nt)
				fprintf(fp,"%g%s",is_arc(i,t)>=0?1.0:0.0,t==nv?"],\n [":",");
			else
				fprintf(fp,"%g%s",is_arc(i,t)>=0?1.0:0.0,t==nv?"]]\n":",");
		}
	}
	fprintf(fp, "[[");
	for (int s = 1; s <= ns; s++) {
		for (int t = nv-nt+1; t <= nv; t++) {
			if (s < ns)
				fprintf(fp,"%g%s",is_arc(s,t)>=0?1.0:0.0,t==nv?"],\n [":",");
			else
				fprintf(fp,"%g%s",is_arc(s,t)>=0?1.0:0.0,t==nv?"]]\n":",");
		}
	}
	fclose(fp);
	 **************************************************************************/
	p_indeg = new vector<int>[ni+1];
	for (int i = ns+1; i <= nv-nt; i++) {
		for (int s = 1; s <= ns; s++) {
			if (is_path(s,i))
				p_indeg[i-ns].push_back(s);
		}
	}
	if (mdl == 5) {
		int nlt = 0;
		for (int l = ns+1; l <= nv-nt; l++) {
			nlt += A[l].size();
			for (int k = 1; k <= nk; k++) {
				list<double> qual;
				for (int s = 1; s <= ns; s++) {
					if (s_p(s,l)>=0)
						qual.push_back(N[s].q[k]);
				}
				qual.sort();
				double lb = qual.front(), ub = qual.back();
				N[l].wy[k-1].setLB(lb);
				N[l].wy[k-1].setUB(ub);
			}
		}
		nlt = nlt*nk;
		nlterms = nlt, nlconsts = nlt;
	}
	if (!isSpp && (1 <= mdl && mdl <= 3)) {
		fprintf(fo, "This model can't handle such network\n");
		exit(0);
	}
	return ret;
}

/* find random int less than a given fraction */
int randf(double p) {
	return drand48() < p ? 1 : 0;
}
/* find random double in a given interval */
double randi(double a, double b) {
	return (a + drand48()*(b - a));
}

/* generate node j randomly to the arc (i,j) */
int graph::rnd_node(int i, int net) {
	int j;
	if (i <= ns)
		j = (int)randi(ns+1, nv+1);
	else if (net == 0) {
		j = (int)randi(ns+1, nv+1);
		while (i == j)
			j = (int)randi(ns+1, nv+1);
	}
	else if (net == 1)
		j = (int)randi(i+1, nv+1);
	else
		j = (int)randi(nv-nt+1, nv+1);
	return j;
}

/******************************************************************************
 *  NAME
 *  random_graph - generate random pooling problem data to graph
 *
 *  SYNOPSIS
 *  int random_graph(int *val, double p, const char *fname);
 *
 *  DESCRIPTION
 *  The routine random_graph generates random pooling problem data by
 *  specifying the number of nodes, sources, terminals, and qualities
 *  given in the array val and the arc's probability density p.
 *
 *  RETURNS
 *  If the operation was successful, the routine returns zero. Otherwise
 *  it prints an error message and returns non-zero.
 *****************************************************************************/

int graph::random_graph(int *val, double p, const char *fname) {
	srand48(time(NULL));
	srand(time(NULL));
	mdl = 0;
	int net, ret = 0;
	double c, bl, bu, q;
	/* read pre-specified network configuration */
	nv = val[0], ns = val[1], nt = val[2], nk = val[3]; net = val[4];
	fprintf(fo, "\nGenerating random %spooling problem data...\n", (net < 2)? "generalized ":"standard ");
	if (net < 2)
		fprintf(fo, "Cycles are %sallowed.\n", (net < 1)? "":"not ");
	name = file_name(fname).c_str();
	if (!(nv >= 0)) {
		fprintf(fo, "number of nodes invalid\n");
		ret = 1;
		return ret;
	}
	if (!( ns >= 0 && ns < nv)) {
		fprintf(fo, "number of source nodes invalid\n");
		ret = 1;
		return ret;
	}
	if (!(nt >= 0 && ns+nt < nv)) {
		fprintf(fo, "number of terminal nodes invalid");
		ret = 1;
		return ret;
	}
	if (!(nk >= 0)) {
		fprintf(fo, "number of quality attributes invalid\n");
		ret = 1;
		return ret;
	}
	N.assign(nv+1, vertex());
	double *cij = new double[nv+1];
	for (int i = 1; i <= nv; i++) {
		bl = 0.0;                    /* [ 0,10 ] */
		bu = randi(20,60);           /* [100,200] */
		N[i] = vertex(this, i, bl, bu);
		cij[i] = 0.0;
		if (i <= ns) {
			cij[i] = randi(0,6);     /* [ 10,50 ] */
			(N[i].q).assign(nk+1, 0.0);
			for (int k = 1; k <= nk; k++) {
				q = randi(0,10);     /* [ 0,80] */
				N[i].q[k] = q;
			}
		}
		if (i > nv-nt) {
			cij[i] = randi(5,15);    /* [ 40,50 ] */
			(N[i].q).assign(nk+1, 0.0);
			for (int k = 1; k <= nk; k++) {
				q = randi(2,7);  /* [ 20,50 ] */;
				N[i].q[k] = q;
			}
		}
	}
	/* generate random network connectivity */
	A = new vector<arc>[nv+1];
	vector <int> nra;
	nra.assign(nv+1, 0);
	for (int i = 1; i <= nv-nt; i++) {
		int j = rnd_node(i, net);
		c = cij[i] - cij[j];
		arc edge(this, i, j, c);
		A[i].push_back(edge);
		nra[j]++;
		double pr = p;
		while(randf(pr) == 1) {
			j = rnd_node(i, net);
			if (is_arc(i,j) < 0) {
				c = cij[i] - cij[j];
				arc edge(this, i, j, c);
				A[i].push_back(edge);
				nra[j]++;
			}
			pr = pr/2;
		}
	}
	if (net < 2) {
		int i = (int)randi(ns+1,nv-nt+1);
		int j = rnd_node(i,net);
		c = cij[i] - cij[j];
		arc edge(this, i, j, c);
		A[i].push_back(edge);
		nra[j]++;
	}
	for (int i = ns+1; i <= nv; i++) {
		if (ns < i && i <= nv-nt) {
			while(nra[i] < 2) {
				int j = (int)randi(1, ns+1);
				if (is_arc(j,i) < 0) {
					c = cij[j] - cij[i];
					arc edge(this, j, i, c);
					A[j].push_back(edge);
					nra[i]++;
				}
			}
		}
		if (i > nv-nt && nra[i] < 1) {
			int j = (int)randi(1, nv-nt+1);
			c = cij[j] - cij[i];
			arc edge(this, j, i, c);
			A[j].push_back(edge);
			nra[i]++;
		}
	}
	delete[] cij;
	na = 0;
	for (int i = 1; i <= nv; i++)
		na += A[i].size();
	fprintf(fo, "Done.\n");
	return ret;
}

/******************************************************************************
 *  NAME
 *  write_graph - write pooling problem data in DIMACS format
 *
 *  SYNOPSIS
 *  int write_graph(const char *fname);
 *
 *  DESCRIPTION
 *  The routine write_graph writes pooling problem data in DIMACS
 *  format to a text file.
 *
 *  RETURNS
 *  If the operation was successful, the routine returns zero. Otherwise
 *  it prints na error message and returns non-zero.
 *****************************************************************************/

int graph::write_graph(const char *fname) {
	FILE *fp;
	int count = 0, ret = 0;
	fprintf(fo, "\nWriting the pooling problem in DIMACS format to '%s.dpp'...\n", fname);
	fp = fopen(get_str("%s.dpp", fname).c_str(), "w");
	if (fp == NULL)	{
		fprintf(fo, "Unable to create `%s'\n", fname);
		ret = 1;
		return ret;
	}
	fprintf(fp, "c %s.dpp, Date: ", name.empty() ? "unknown" : name.c_str());
	fprintf(fp, "%s", cdate()), count++;
	fprintf(fp, "c\n"), count++;
	fprintf(fp, "c The pooling problem data in DIMACS format.\n"), count++;
	fprintf(fp, "c\n"), count++;
	fprintf(fp, "p min %d %d %d %d %d\n", na, nv, ns, nt, nk), count++;
	fprintf(fp, "c\n"), count++;
	for (int i = 1; i <= nv; i++) {
		if (ns < i && i <= nv-nt)
			fprintf(fp, "n %d %.2g", N[i].v, N[i].bu);
		if (i <= ns) {
			fprintf(fp, "n %d %.2g %.2g", N[i].v, N[i].bl, N[i].bu);
			for (int k = 1; k <= nk; k++)
				fprintf(fp, " %.2g", N[i].q[k]);
		}
		if (i > nv-nt) {
			fprintf(fp, "n %d %.2g %.2g", N[i].v, N[i].bl, N[i].bu);
			for (int k = 1; k <= nk; k++)
				fprintf(fp, " %.2g", N[i].q[k]);
		}
		fprintf(fp, "\n"), count++;
	}
	fprintf(fp, "c\n"), count++;
	for (int i = 1; i <= nv; i++) {
		for (int j = 0; j < A[i].size(); j++)
			fprintf(fp, "a %d %d %.2g\n", i, A[i][j].j, A[i][j].c), count++;
	}
	fprintf(fp, "c\n"), count++;
	fprintf(fp, "c eof"), count++;
	fflush(fp);
	if (ferror(fp))	{
		printf("Write error on `%s'\n", fname);
		ret = 1;
		return ret;
	}
	fclose(fp);
	fprintf(fo, "%d lines were written\n", count);
	return ret;
}

/******************************************************************************
 *  NAME
 *  write_gams - write pooling problem data in GAMS file
 *
 *  SYNOPSIS
 *  int write_gams(const char *fname);
 *
 *  DESCRIPTION
 *  The routine write_gams writes pooling problem data in a GAMS
 *  format to a text file.
 *
 *  RETURNS
 *  If the operation was successful, the routine returns zero. Otherwise
 *  it prints an error message and returns non-zero.
 *****************************************************************************/

int graph::write_gams(const char *fname) {
	FILE *fp;
	int nd, count = 0, ret = 0;
	fprintf(fo, "\nWriting the pooling problem in GAMS format to '%s.gms'...\n", fname);
	fp = fopen(get_str("%s.gms", fname).c_str(), "w");
	if (fp == NULL)	{
		printf("Unable to create `%s'\n", fname);
		ret = 1;
		return ret;
	}
	nd = (int)log10(nv)+1;
	fprintf(fp, "$ontext\n"), count++;
	fprintf(fp, "    %s ", name.empty() ? "unknown" : name.c_str());
	fprintf(fp, "pooling problem data.\n"), count++;
	fprintf(fp, "    Author: Mohammed Alfaki, %s", cdate()), count++;
	fprintf(fp, "$offtext\n\n"), count++, count++;
	fprintf(fp, "$eolcom #\n\n"), count++, count++;
	fprintf(fp, "# Declare sets\n"); count++;
	fprintf(fp, "    set i    /1*%d/;\n", nv), count++;
	fprintf(fp, "    set s(i) /1*%d/;\n", ns), count++;
	fprintf(fp, "    set t(i) /%d*%d/;\n", nv-nt+1, nv), count++;
	fprintf(fp, "    set k    /1*%d/;\n\n", nk), count++, count++;
	fprintf(fp, "alias (i,j);\n\n"), count++, count++;
	fprintf(fp, "# The arc unit cost c_{ij}\n"), count++;
	fprintf(fp, "table c(i,j)\n %*s", nd, ""), count++;
	for (int j = ns+1; j <= nv; j++)
		fprintf(fp, " %7d", j);
	for (int i = 1; i <= nv-nt; i++) {
		fprintf(fp, "\n"), count++;
		for (int j = ns+1; j <= nv; j++) {
			double c = 0;
			for (int it = 0; it < A[i].size(); it++) {
				if (A[i][it].j == j)
					c = A[i][it].c;
			}
			if (j == ns+1)
				fprintf(fp, " %*d", nd, i);
			fprintf(fp, " %7.2f", c);
		}
	}
	fprintf(fp, " ;\n\n"), count++, count++;
	fprintf(fp, "# The adjacency matrix (the arcs set A)\n"), count++;
	fprintf(fp, "table a(i,j)\n %*s", nd, ""), count++;
	for (int j = ns+1; j <= nv; j++)
		fprintf(fp, " %*d", nd, j);
	for (int i = 1; i <= nv-nt; i++) {
		fprintf(fp, "\n"), count++;
		for (int j = ns+1; j <= nv; j++) {
			if (j == ns+1)
				fprintf(fp, " %*d", nd, i);
			fprintf(fp, " %*d", nd, is_arc(i,j)>= 0 ? 1 : 0);
		}
	}
	fprintf(fp, " ;\n\n"), count++, count++;
	fprintf(fp, "# Source qualities/terminal quality upper bounds\n"), count++;
	fprintf(fp, "table q(i,k)\n %*s", nd, ""), count++;
	for (int k = 1; k <= nk; k++)
		fprintf(fp, " %7d", k);
	for (int i = 1; i <= nv; i++) {
		if (i <= ns || i > nv-nt) {
			fprintf(fp, "\n"), count++;
			for (int k = 1; k <= nk; k++) {
				if (k == 1)
					fprintf(fp, " %*d", nd, i);
				fprintf(fp, " %7.2f", N[i].q[k]);
			}
		}
	}
	fprintf(fp, " ;\n\n"), count++, count++;
	fprintf(fp, "# Node capacity lower bound\n"), count++;
	fprintf(fp, "parameter bl(i) /%*d  %7.2f\n", nd, 1, N[1].bl), count++;
	for (int i = 2; i < nv; i++) {
		if (ns < i && i <= nv-nt) continue;
		fprintf(fp, "                 %*d  %7.2f\n", nd, i, N[i].bl), count++;
	}
	fprintf(fp, "                 %*d  %7.2f / ;\n\n", nd, nv, N[nv].bl);
	count++, count++;
	fprintf(fp, "# Node capacity upper bound\n"), count++;
	fprintf(fp, "parameter bu(i) /%*d  %7.2f\n", nd, 1, N[1].bu), count++;
	for (int i = 2; i < nv; i++) {
		fprintf(fp, "                 %*d  %7.2f\n", nd, i, N[i].bu);
		count++;
	}
	fprintf(fp, "                 %*d  %7.2f / ;\n\n", nd, nv, N[nv].bu);
	count++, count++;
	fprintf(fp, "$include xmodel.gms"), count++;
	fflush(fp);
	if (ferror(fp))	{
		printf("Write error on `%s'\n", fname);
		ret = 1;
		return ret;
	}
	fprintf(fo, "%d lines were written\n", count);
	if (fp != NULL)
		fclose(fp);
	return ret;
}

class var {
public:
	int id;
	string name;
	var(int i) { id = i, name = get_str("x(%d)", i); };
};

int graph::write_matlab(const char *fname, int flag) {
	int ret = 0;
	string file = fname;
	file += (flag == 1)?"p":"";
	fprintf(fo, "\nExpress the pooling problem in Matlab format:");
	fprintf(fo, "\nWriting the objective function in obj%s.m...", file.c_str());
	FILE *fp = fopen(("mat/obj"+file+".m").c_str(), "w");
	if (fp == NULL) {
		printf("Unable to create `obj%s'\n", file.c_str());
		ret = 1;
		return ret;
	}
	/**************************************************************************/
	fprintf(fp, "function f = obj%s(x)\n", file.c_str());
	fprintf(fp, "if strcmp(x,'init')\n");
	int numOfVar = 0;
	string str, cExpr, LB, UB;
	vector <var> *wy;
	if (flag == 1) {
		LB = UB = "[";
		wy = new vector<var>[nv-ns-nt+1];
		for (int i = ns+1; i <= nv-nt; i++) {
			for (int k = 1; k <= nk; k++) {
				list<double> qual;
				for (int s = 1; s <= ns; s++) {
					if (is_path(s,i))
						qual.push_back(N[s].q[k]);
				}
				qual.sort();
				double lb = qual.front(), ub = qual.back();
				numOfVar++;
				wy[i-ns].push_back(var(numOfVar));
				LB += get_str("%s%g", numOfVar == 1?"":" ", lb);
				UB += get_str("%s%g", numOfVar == 1?"":" ", ub);
			}
		}
	}
	else {
		UB = "[";
		wy = new vector<var>[ns+1];
		for (int s = 1; s <= ns; s++) {
			wy[s].assign(nv-nt+1, var(-1));
			for (int i = ns+1; i <= nv-nt; i++) {
				if (is_path(s,i)) {
					numOfVar++;
					wy[s][i] = var(numOfVar);
					UB += get_str("%s%g", numOfVar == 1?"":" ", 1.0);
				}
			}
		}
	}
	vector <var> *f;
	f = new vector<var>[nv-nt+1];
	for (int i = 1; i <= nv-nt; i++) {
		for (int jp = 0; jp < A[i].size(); jp++) {
			int j = A[i][jp].j;
			double ub = fmin(N[i].bu, N[j].bu);
			numOfVar++;
			f[i].push_back(var(numOfVar));
			LB += get_str(" %g", 0.0);
			UB += get_str(" %g", ub);
		}
	}
	if (flag == 1) {
		LB += "];";
		UB += "];";
	}
	else {
		LB = get_str("zeros(1,%d);", numOfVar);
		UB += "];";
	}
	fprintf(fp, "    f.LB = %s\n", LB.c_str());
	fprintf(fp, "    f.UB = %s\n", UB.c_str());
	int c = 0;
	vector<int> I, J;
	vector<double> S, b;
	for (int i = 1; i <= nv; i++) {
		if (N[i].bl > 0 && !(ns < i && i <= nv-nt)) {
			c++;
			if (i <= ns) {
				for (int j = 0; j < A[i].size(); j++) {
					I.push_back(c), J.push_back(f[i][j].id);
					S.push_back(-1.0);
				}
			}
			else {
				for (int j = 1; j <= nv; j++) {
					int ip = is_arc(j,i);
					if (ip>=0) {
						I.push_back(c), J.push_back(f[j][ip].id);
						S.push_back(-1.0);
					}
				}
			}
			b.push_back(-N[i].bl);
		}
	}
	for (int i = 1; i <= nv; i++) {
		if (N[i].bu < inf) {
			c++;
			if (i <= ns) {
				for (int j = 0; j < A[i].size(); j++) {
					I.push_back(c), J.push_back(f[i][j].id);
					S.push_back(1.0);
				}
			}
			else {
				for (int j = 1; j <= nv; j++) {
					int ip = is_arc(j,i);
					if (ip>=0) {
						I.push_back(c), J.push_back(f[j][ip].id);
						S.push_back(1.0);
					}
				}
			}
			b.push_back(N[i].bu);
		}
	}
	fprintf(fp, "    i = [");
	for (int i = 0; i < I.size(); i++)
		fprintf(fp, "%s%d", i==0?"":" ", I[i]);
	fprintf(fp, "];\n");
	fprintf(fp, "    j = [");
	for (int i = 0; i < I.size(); i++)
		fprintf(fp, "%s%d", i==0?"":" ", J[i]);
	fprintf(fp, "];\n");
	fprintf(fp, "    s = [");
	for (int i = 0; i < I.size(); i++)
		fprintf(fp, "%s%g", i==0?"":" ", S[i]);
	fprintf(fp, "];\n");
	fprintf(fp, "    f.Aineq = full(sparse(i,j,s,%d,%d));\n", c, numOfVar);
	fprintf(fp, "    f.bineq = [");
	for (int i = 0; i < b.size(); i++)
		fprintf(fp, "%s%g", i==0?"":" ", b[i]);
	fprintf(fp, "]';\n");
	I.clear(), J.clear(), S.clear(), b.clear(); c = 0;
	/**************************************************************************/
	if (flag == 1) {
		for (int i = ns+1; i <= nv-nt; i++) {
			c++;
			for (int j = 1; j <= nv-nt; j++) {
				int ip = is_arc(j,i);
				if (ip>=0) {
					I.push_back(c), J.push_back(f[j][ip].id);
					S.push_back(1.0);
				}
			}
			for (int j = 0; j < A[i].size(); j++) {
				I.push_back(c), J.push_back(f[i][j].id);
				S.push_back(-1.0);
			}
			b.push_back(0.0);
		}
	}
	else {
		for (int i = ns+1; i <= nv-nt; i++) {
			c++;
			for (int s = 1; s <= ns; s++) {
				if (is_path(s,i)) {
					I.push_back(c), J.push_back(wy[s][i].id);
					S.push_back(1.0);
				}
			}
			b.push_back(1.0);
		}
	}
	fprintf(fp, "    i = [");
	for (int i = 0; i < I.size(); i++)
		fprintf(fp, "%s%d", i==0?"":" ", I[i]);
	fprintf(fp, "];\n");
	fprintf(fp, "    j = [");
	for (int i = 0; i < I.size(); i++)
		fprintf(fp, "%s%d", i==0?"":" ", J[i]);
	fprintf(fp, "];\n");
	fprintf(fp, "    s = [");
	for (int i = 0; i < I.size(); i++)
		fprintf(fp, "%s%g", i==0?"":" ", S[i]);
	fprintf(fp, "];\n");
	fprintf(fp, "    f.Aeq = full(sparse(i,j,s,%d,%d));\n", c, numOfVar);
	fprintf(fp, "    f.beq = [");
	for (int i = 0; i < b.size(); i++)
		fprintf(fp, "%s%g", i==0?"":" ", b[i]);
	fprintf(fp, "]';\n");
	I.clear(), J.clear(), S.clear(), b.clear();
	/**************************************************************************/
	fprintf(fp, "    f.fitnessfcn = @obj%s;\n", file.c_str());
	fprintf(fp, "    f.nonlcon = 'nlc%s';\n", file.c_str());
	fprintf(fp, "    f.nvars = %d;\n", numOfVar);
	fprintf(fp, "    f.options.PopulationSize = 50;\n");
	fprintf(fp, "    f.options.Generations = 10;\n");
	fprintf(fp, "    f.options.ConstrBoundary = 'absorb';\n");
	fprintf(fp, "else\n");
	/**************************************************************************/
	cExpr = "";
	for (int i = 1; i <= nv-nt; i++) {
		for (int j = 0; j < A[i].size(); j++)
			cExpr += A[i][j].c == 0?"":(get_str("%+g", A[i][j].c)+"*"+f[i][j].name);
	}
	fprintf(fp, "    f = %s;\n", cExpr.c_str());
	fprintf(fp, "end\n");
	fflush(fp);
	if (ferror(fp)) {
		printf("Write error on `obj%s'\n", file.c_str());
		ret = 1;
		return ret;
	}
	if (fp != NULL)
		fclose(fp);
	/**************************************************************************/
	fprintf(fo, "\nWriting the nonlinear constraints in nlc%s.m...\n", file.c_str());
	fp = fopen(("mat/nlc"+file+".m").c_str(), "w");
	if (fp == NULL) {
		printf("Unable to create `nlc%s'\n", file.c_str());
		ret = 1;
		return ret;
	}
	fprintf(fp, "function [c, ceq] = nlc%s(x)", file.c_str());
	fprintf(fp, "\nc   = [");
	if (flag == 1) {
		for (int i = nv-nt+1; i <= nv; i++) {
			for (int k = 1; k <= nk; k++) {
				string cExpr = "";
				for (int s = 1; s <= ns; s++) {
					int ip = is_arc(s,i);
					if (ip>=0) {
						double coef = (N[s].q[k]-N[i].q[k]);
						cExpr += coef==0?"":(get_str("%+g", coef)+"*"+f[s][ip].name);
					}
				}
				for (int l = ns+1; l <= nv-nt; l++) {
					int ip = is_arc(l,i);
					if (ip>=0)
						cExpr += "+"+f[l][ip].name+"*"+wy[l-ns][k-1].name;
				}
				for (int j = ns+1; j <= nv-nt; j++) {
					int ip = is_arc(j,i);
					if (ip>=0) {
						double coef = -1*N[i].q[k];
						cExpr += coef==0?"":(get_str("%+g", coef)+"*"+f[j][ip].name);
					}
				}
				if (0 < fabs(N[i].q[k]) && fabs(N[i].q[k]) < inf) {
					str = (i==nv-nt+1 && k==1)?"":"\n       ";
					fprintf(fp, "%s%s", str.c_str(), cExpr.c_str());
				}
			}
		}
	}
	else {
		for (int i = nv-nt+1; i <= nv; i++) {
			for (int k = 1; k <= nk; k++) {
				string cExpr = "";
				for (int s = 1; s <= ns; s++) {
					for (int j = ns+1; j <= nv-nt; j++) {
						int ip = is_arc(j,i);
						if (is_path(s,j) && ip>=0) {
							double coef = (N[s].q[k]-N[i].q[k]);
							cExpr += coef==0?"":(get_str("%+g", coef)+"*"+wy[s][j].name+"*"+f[j][ip].name);
						}
					}
					int ip = is_arc(s,i);
					if (ip>=0) {
						double coef = (N[s].q[k]-N[i].q[k]);
						cExpr += coef==0?"":(get_str("%+g", coef)+"*"+f[s][ip].name);
					}
				}
				if (0 < fabs(N[i].q[k]) && fabs(N[i].q[k]) < inf) {
					str = (i==nv-nt+1 && k==1)?"":"\n       ";
					fprintf(fp, "%s%s", str.c_str(), cExpr.c_str());
				}
			}
		}
	}
	fprintf(fp, "];\n");
	fprintf(fp, "\nceq = [");
	if (flag == 1) {
		for (int i = ns+1; i <= nv-nt; i++) {
			for (int k = 1; k <= nk; k++) {
				string cExpr = "";
				for (int s = 1; s <= ns; s++) {
					int ip = is_arc(s,i);
					if (ip>=0) {
						double coef = N[s].q[k];
						cExpr += coef==0?"":(get_str("%+g", coef)+"*"+f[s][ip].name);
					}
				}
				for (int l = ns+1; l <= nv-nt; l++) {
					int ip = is_arc(l,i);
					if (ip>=0)
						cExpr += "+"+f[l][ip].name+"*"+wy[l-ns][k-1].name;
				}
				for (int j = 0; j < A[i].size(); j++)
					cExpr += "-"+f[i][j].name+"*"+wy[i-ns][k-1].name;
				str = (i==ns+1 && k==1)?"":"\n       ";
				fprintf(fp, "%s%s", str.c_str(), cExpr.c_str());
			}
		}
	}
	else {
		int l = 1;
		for (int i = ns+1; i <= nv-nt; i++) {
			for (int s = 1; s <= ns; s++) {
				if (is_path(s,i)) {
					string cExpr = "";
					for (int j = 0; j < A[i].size(); j++)
						cExpr += "+"+wy[s][i].name+"*"+f[i][j].name;
					int ip = is_arc(s,i);
					if (ip>=0)
						cExpr += "-"+f[s][ip].name;
					for (int j = ns+1; j <= nv-nt; j++) {
						int ip = is_arc(j,i);
						if (is_path(s,j) && ip>=0)
							cExpr +="-"+ wy[s][j].name+"*"+f[j][ip].name;
					}
					str = (l == 1)?"":"\n       ";
					fprintf(fp, "%s%s", str.c_str(), cExpr.c_str());
					l++;
				}
			}
		}
	}
	fprintf(fp, "];\n");
	delete[] wy;
	delete[] f;
	fflush(fp);
	if (ferror(fp)) {
		printf("Write error on `nlc%s'\n", file.c_str());
		ret = 1;
		return ret;
	}
	if (fp != NULL)
		fclose(fp);
	return ret;
}

/* eof */
