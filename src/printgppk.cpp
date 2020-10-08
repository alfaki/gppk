#include <cstring>
#include <stdarg.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

#define p(str) (strcmp(shrt, str) == 0)
#define STRING_LENGTH 1000

void usage(const char* name) {
	printf("usage: %s [model_name] <val>\n",name);
	printf("For help: -h help information and usage.\n");
	printf("model_name: p, pq, tp, stp, and mcf-formulation\n");
	printf("            cdg convert dpp to gms\n");
	printf("            cdd convert dpp to dpp\n");
	printf("            dsc the discrete model for the pooling problem\n");
	printf("            hrs Heuristic algorithm for the pooling problem\n\n");
	exit(0);
}

/* Write formatted output to char string */
char *get_cstr(const char *fmt, ...) {
	va_list arg;
	va_start(arg, fmt);
	char *cstr;
	cstr = new char[STRING_LENGTH];
	vsprintf(cstr, fmt, arg);
	va_end(arg);
	return cstr;
}

int main (int argc, char **argv) {
	if ( argc < 2 ) {
		usage(argv[0]);
		return 0;
	}
	int len = strlen(argv[1]), n, ndsc, i = 0, flg;
	bool dpp2gms = false, dpp2dpp = false;
	char shrt[len], str[STRING_LENGTH];
	strcpy(shrt, argv[1]);
	char *problems[STRING_LENGTH], *pch;
	const char *val;
	if p("-h") {
		usage(argv[0]);
		return 0;
	}
	if p("cdg") dpp2gms = true;
	if p("cdd") dpp2dpp = true;
	if ((p("dsc") || p("hrs")) && argc < 3) {
		printf("%d here\n", argc);
		usage(argv[0]);
		return 0;
	}
	val = (p("dsc") || p("hrs")) ? argv[2]:"";
	ifstream file("inst.txt");
	if(!file) {
		printf("\nNo file or the file is impty!\n");
		exit(0);
	}
	file.getline(str,STRING_LENGTH);
	pch = strtok(str, " ");
	while ( pch != NULL ) {
		problems[i] = pch;
		pch = strtok (NULL, " ");
		i = i+1;
	}
	n = i;
	file.close();
	FILE *f1, *f2;
	char ln[] = "===================================";
	if (!dpp2gms && !dpp2dpp) {
		f1 = fopen(get_cstr("%s%s_inst_sol.tab", shrt, val), "w");
		//f1 = stdout;
		fprintf(f1, "Model name: %smodel %s\n", shrt, val);
		fprintf(f1, "%s%s\n", ln, ln);
		fprintf(f1, "%*s %*s %*s %*s %*s\n", 9, "instance", 8, "nods", 15,
						   "lb (time)", 15, "ub", 15, "rgap");
		fprintf(f1, "%s%s\n", ln, ln);
		f2 = fopen(get_cstr("%s%s_inst_prob.tab", shrt, val), "w");
		fprintf(f2, "Model name: %smodel %s\n", shrt, val);
		fprintf(f2, "%s========\n", ln);
		fprintf(f2, "%*s %*s %*s %*s %*s\n", 9, "instance", 7, "vars", 7,
						    p("dsc")?"invrs":"nlts", 7, "lcs", 7, "nlcs");
		fprintf(f2, "%s========\n", ln);
	}
	for (int i = 0; i < n; i++) {
		if (dpp2dpp || dpp2gms) {
			flg = system(get_cstr("./gppk -o %s dpp/%s.dpp inst/%s", shrt,
					problems[i], problems[i]));
			continue;
		}
		else {
			flg = system(get_cstr("./gppk -o %s dpp/%s.dpp %s", shrt,
														problems[i], val));
			if p("hrs") 
				flg = system(get_cstr("mv gppk.log %s/%s.log", val, problems[i]));
			ifstream file(get_cstr("out/%s%s.out", problems[i], p("dsc")?val:""));
			string line1, line2;
			getline(file, line1);
			getline(file, line2);
			file.close();
			char line[STRING_LENGTH];
			strcpy(line, line1.c_str());
			double ub, lb, time, rgap, agap;
			ub = strtod(line, &pch);
			lb = strtod(pch, &pch);
			time = strtod(pch, &pch);
			rgap = strtod(pch, &pch);
			agap = strtod(pch, NULL);
			strcpy(line, line2.c_str());
			int solstat, nods, vars, nlts, lcs, nlcs, intvars, mdl;
			solstat = strtol(line, &pch,10);
			nods = strtol(pch, &pch,10);
			vars = strtol(pch, &pch,10);
			nlts = strtol(pch, &pch,10);
			lcs = strtol(pch, &pch,10);
			nlcs = strtol(pch, &pch,10);
			intvars = strtol(pch, &pch,10);
			mdl = strtol(pch, NULL,10);
			fprintf(f2, "%*s %7d %7d %7d %7d\n", 9, problems[i], vars,
							    p("dsc")? intvars:nlts, lcs, nlcs);
			if (solstat == 1) {
				fprintf(f1, "%*s %8d %*s) %15.2f %15.6f\n", 9, problems[i], nods,
						14, get_cstr("(%.2f", time), ub, rgap);
			}
			else {
				fprintf(f1, "%*s %8d %15.2f %15.2f %15.6f   *\n", 9, problems[i],
						nods, lb, ub, rgap);
			}
		}
	}
	if (f1 != NULL) {
		fprintf(f1, "%s%s\n", ln, ln);
		fclose(f1);
	}
	if (f2 != NULL) {
		fprintf(f2, "%s========\n", ln);
		fclose(f2);
	}
	return 0;
}
