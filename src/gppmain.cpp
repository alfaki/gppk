/* gppmain.c */

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

#define p(str) (strcmp(argv[k], str) == 0)

/* global object */
graph *G;

/* catch interrupt signal from the keyboard */
void signal_handler(int sgnl) {
	G->print_summary(sgnl);
	exit(1);
}

/* check if str is number or not */
bool strIsNum(char *str) {
	bool ret = true;
	int i = 0;
	while (str[i]) {
		if (!isdigit(str[i])) {
			ret = false;
			break;
		}
		i++;
	}
	return ret;
}

bool fexists(const char *filename) {
	ifstream ifile(filename);
	return ifile;
}

/* the main function */
int main(int argc, char **argv) {
	FILE *fo = stdout;
	const char *in_file, *out_file;
	int mdl = 6, val[5], n = 1, flag = 1, o = 0, k = 1;
	double p;
	if (argc <= 2 || argv[k][0] == '\0' || p("-h") || p("-help")) {
		print_usage (argv[0]);
		return 0;
	}
	if p("-o") {
		fo = fopen(get_str("%s.log", argv[0]).c_str(), "w");
		o++;
		k++;
	}
	else if (argc < 3) {
		print_usage (argv[0]);
		return 0;
	}
	if      p("rnd") mdl = 0;
	else if p("pq")  mdl = 1;
	else if p("tp")  mdl = 2;
	else if p("stp") mdl = 3;
	else if p("mcf") mdl = 4;
	else if p("p")   mdl = 5;
	else if p("dsc") mdl = 6;
	else if p("hrs") mdl = 7;
	else if p("cdg") mdl = 8;
	else if p("cdd") mdl = 9;
	else if p("mat") mdl = 10;
	else mdl = 4; // default
	k++;
	in_file = argv[k++];
	if (1 <= mdl && mdl <= 5 && argc-o > 2)
		fprintf(fo, "Starting Branch-and-bound...\n");
	else if (mdl == 0 && argc-o > 8) {
		val[0] = strtol(argv[k++], NULL, 10);
		val[1] = strtol(argv[k++], NULL, 10);
		val[2] = strtol(argv[k++], NULL, 10);
		val[3] = strtol(argv[k++], NULL, 10);
		p = strtod(argv[k++], NULL);
		val[4] = strtol(argv[k], NULL, 10);
	}
	else if (mdl == 6 && argc-o > 3) {
		if (!strIsNum(argv[k])) {
			fprintf(fo, "The number of discretized intervals is not valid!\n");
			print_usage(argv[0]);
			return 0;
		}
		n = strtol(argv[k], NULL, 10);
	}
	else if (mdl == 7 && argc-o > 3) {
		if (!fexists(get_str("%shrs.gms", argv[k]).c_str())) {
			fprintf(fo, "The GAMS model file to solve the bilinear problem is not exists!\n");
			print_usage(argv[0]);
			return 0;
		}
		fprintf(fo, "Starting Heuristic with %s-model...\n", argv[k]);
		FILE *fmodel = fopen("xmodel.gms", "w");
		fprintf(fmodel, "$include %shrs.gms\n", argv[k++]);
		fclose(fmodel);
	}
	else if ((8 <= mdl && mdl <= 9) && argc-o > 3)
		out_file = argv[k];
	else if (mdl == 10 && argc-o > 4) {
		out_file = argv[k++];
		if (!strIsNum(argv[k])) {
			fprintf(fo, "Enter either 0 for P-model or 1 for MCF-model!\n");
			print_usage(argv[0]);
			return 0;
		}
		flag = strtol(argv[k], NULL, 10);
	}
	else {
		print_usage (argv[0]);
		return 0;
	}
	clock_t start;
	double ttime = 0.0;
	G = new graph();
	G->fo = fo;
	if (mdl == 0) {
		G->random_graph(val, p, in_file);
		G->write_graph(in_file);
		G->write_gams(in_file);
	}
	else {
		G->read_graph(in_file, mdl);
		struct sigaction sigIntHandler;
		sigIntHandler.sa_handler = signal_handler;
		sigemptyset(&sigIntHandler.sa_mask);
		sigIntHandler.sa_flags = 0;
		sigaction(SIGINT, &sigIntHandler, NULL);
		if (1 <= mdl && mdl <= 5) {
			try {
				G->build_relaxation();
				start = clock();
				G->check_feasibility();
				G->branch_and_bound(start);
			}
			catch (IloException& ex) {
				cerr << "Error: " << ex << endl;
			}
		}
		if (mdl == 6) G->discrete_model(n, 1, TIME);
		if (mdl == 7) G->heuristic_alg();
		if (mdl == 8) G->write_gams(out_file);
		if (mdl == 9) G->write_graph(out_file);
		if (mdl == 10) G->write_matlab(out_file, flag);
	}
	G->print_summary(-1);
	delete G;
	fclose(fo);
	return 0;
}

/* print problem information */
void graph::problem_info() {
	const char *modl[] = {"Random data", "PQ", "TP", "STP", "MCF", "P", "Discretized", "Heuristic"};
	int sum, nline = 1;
	fprintf(fo, "\nReading problem...\n");
	fprintf(fo, "The %s-model of '%s' has %d variable%s,\n", modl[mdl], name.c_str(), vars, vars == 1 ? "" : "s");
	if(intvars > 0) {
		sum = nlterms+lconsts+nlconsts;
		fprintf(fo, "%s%d binary variable%s%s%s", sum == 0 ? "and " : "", intvars, intvars == 1 ? "" : "s", sum == 0 ? ".\n" : ", ", nline%2 == 1 ? "" : "\n");
		nline++;
	}
	if(nlterms > 0) {
		sum = lconsts+nlconsts;
		fprintf(fo, "%s%d nonlinear term%s%s%s", sum == 0 ? "and " : "", nlterms, nlterms == 1 ? "" : "s", sum == 0 ? ".\n" : ", ", nline%2 == 1 ? "" : "\n");
		nline++;
	}
	if(lconsts > 0) {
		sum = nlconsts;
		fprintf(fo, "%s%d linear constraint%s%s%s", sum == 0 ? "and " : "", lconsts, lconsts == 1 ? "" : "s", sum == 0 ? "." : ", ", nline%2 == 1 ? "" : "\n");
		nline++;
	}
	if(nlconsts > 0)
		fprintf(fo, "and %d nonlinear constraint%s.\n", nlconsts, nlconsts == 1 ? "" : "s");
}


/* print solution summary */
void graph::print_summary(int sgnl) {
	if (sgnl >= 0) fprintf(fo, "\nThe solver is interrupted!\n");
	else {
		if (mdl == 0)
			fprintf(fo, "Random pooling problem generated successfully.\n");
		if (mdl >= 8)
			fprintf(fo, "The file converted successfully.\n");
	}
	if (1 <= mdl && mdl <= 7) {
		double agap = fabs(gpp_ub-gpp_lb);
		double gap = agap/fabs(gpp_ub+1e-10);
		FILE *fp = fopen(("out/"+(name.empty()?"unknown":name)+".out").c_str(), "w");
		fprintf(fp, "%g %g %g %g %g\n", gpp_ub, gpp_lb, ttime, gap, agap);
		fprintf(fp, "%d %d %d %d %d %d %d %d\n", solstat, nnode, vars, nlterms, lconsts, nlconsts, intvars, mdl);
		fprintf(fp, "/* First row: ub, lb, time, rgap, agap */\n");
		fprintf(fp, "/* Second row: solstat, #node, vars, nlterms, lconsts, nlconsts, intvars, mdl */\n");
		fclose(fp);
		string line = "===============================================";
		fprintf(fo, "\n%s\n", line.c_str());
		fprintf(fo, "Optimal solution found     = %16.6f\n", gpp_ub);
		fprintf(fo, "Best objective value       = %16.6f\n", gpp_lb);
		fprintf(fo, "Relative gap               = %16.6f%%\n", gap);
		fprintf(fo, "Elapsed time (hh:mm:ss:ms) = %*s\n", 16, clocktime(ttime).c_str());
		fprintf(fo, "%s\n", line.c_str());
	}
}

/* print the package usage instructions */
void print_usage(const char* pkg) {
	printf("\nUsage: %s [option] [file] <values>\n\n", pkg);
	printf("for help use -h option for information and usage.\n");
	printf("option:\n");
	printf("      p   P-formulation.\n");
	printf("      pq  TP-formulation.\n");
	printf("      tp  TP-formulation.\n");
	printf("      stp STP-formulation.\n");
	printf("      mcf MCF-formulation.\n");
	printf("      cdg Convert dpp to gms.\n");
	printf("      cdd Convert dpp to dpp.\n");
	printf("      rnd Generate random pooling network.\n");
	printf("      dsc Discretized model for the pooling problem.\n");
	printf("      hrs Heuristic algorithm for the pooling problem.\n");
	printf("      mat Express the pooling problem in Matlab files.\n");
	printf("file: file containing the pooling problem data in DIMACS format.");
	printf("\nvalues (optional) if:\n");
	printf("      'cdg', 'cdd', or 'mat': The converted file name.\n");
	printf("      'dsc': Number of discretized intervals.\n");
	printf("      'hrs': The GAMS model file to solve the sub bilinear problem.\n");
	printf("      'rnd': The values should be in order: number of nodes,\n");
	printf("             sources, terminals, quality attributes, the expected\n");
	printf("             network density and flag that can take the values 0,\n");
	printf("             1, and 2 to choose between generalized with and without\n");
	printf("             cycles and standard pooling problem, respectively.\n\n");
	exit(0);
}

/* eof */
