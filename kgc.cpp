#include "kgc.h"

#define IJ(BUF, I, J) ((BUF)[I * N + J])

template <typename type>
__attribute__((always_inline)) inline
void printvararray(type &ia, IloCplex cplex) {

	for (id i = 0; i < ia.getSize(); i++) {
		try {
			if (cplex.getValue(ia[i]))
				cout << ia[i].getName() << " = " << cplex.getValue(ia[i]) << endl;
		}
		catch (IloException& e) { e.end(); }
	}
	puts("");
}

template <typename type>
__attribute__((always_inline)) inline
void printvarmatrix(type &ia, IloCplex cplex, const char *name = NULL, const char *format = NULL) {

	if (name) printf("%s =\n", name);

	for (id i = 0; i < N; i++) {
		printf("[ ");
		for (id j = 0; j < N; j++) {
			try {
				if (format) { printf(format, cplex.getValue(IJ(ia, i, j))); printf(" "); }
				else std::cout << cplex.getValue(IJ(ia, i, j)) << " ";
			}
			catch (IloException& e) {
				if (format) { printf(format, 0); printf(" "); }
				else std::cout << 0 << " ";
			e.end(); }
		}
		puts("]");
	}

	puts("");
}

int main(int argc, char *argv[]) {

	#ifdef DEBUG
	printf("# Vertices = %u\n# Edges = %u\nMaximum cardinality = %u\n\n", N, E, K);
	#endif

	// Allocate data structures

	value *adj = (value *)calloc(N * N, sizeof(value));
	uint8_t *la = (uint8_t *)calloc(N, sizeof(uint8_t));

	// Read input file

	#define MAXLINE 1000
	char line[MAXLINE];
	FILE *f = fopen(argv[1], "r");

	for (id i = 0; i < N; i++) {
		fgets(line, MAXLINE, f);
		char *pch = line;
		if (*pch == '*') {
			pch++;
			la[i] = 1;
		}
		IJ(adj, i, i) = atof(pch);
	}

	for (id i = 0; i < E; i++) {
		id v1, v2;
		value w;
		fscanf(f, "%u %u %f", &v1, &v2, &w);
		IJ(adj, v1, v2) = IJ(adj, v2, v1) = w;
	}

	fclose(f);

	#ifdef DEBUG
	puts("Weights:");
	for (id i = 0; i < N; i++)
		printbuf(adj + i * N, N, NULL, "%+.4f");
	puts("");
	printbuf(la, N, "Leaders", "%hu", "\n");
	puts("");
	#endif

	// Create CPLEX model

	IloEnv env;
	IloIntVarArray xa(env, N * N);
	IloModel model(env);
	ostringstream ostr;

	// Create decision variables

	for (id i = 0; i < N; i++)
		for (id j = 0; j < N; j++) {
			ostr << "X_" << i << "," << j;
			IJ(xa, i, j) = IloIntVar(env, 0, 1, ostr.str().c_str());
			ostr.str("");
		}

	// Reflexivity constraints

	for (id i = 0; i < N; i++)
		model.add(IJ(xa, i, i) == 1);

	// Simmetry constraints

	for (id i = 0; i < N; i++)
		for (id j = i + 1; j < N; j++)
			model.add(IJ(xa, i, j) == IJ(xa, j, i));

	// Transitivity constraints

	for (id i = 0; i < N; i++)
		for (id j = i + 1; j < N; j++)
			for (id k = j + 1; k < N; k++) {
				model.add(IJ(xa, i, j) + IJ(xa, j, k) - 2 * IJ(xa, i, k) <= 1);
				model.add(IJ(xa, i, k) + IJ(xa, i, j) - 2 * IJ(xa, j, k) <= 1);
				model.add(IJ(xa, j, k) + IJ(xa, i, k) - 2 * IJ(xa, i, j) <= 1);
			}

	// Cardinality constraints

	IloIntVarArray ka(env, N); // cardinality variables

	#ifdef DEBUG
	puts("Cardinality constraints:");
	#endif

	for (id i = 0; i < N; i++) {

		ostr << "K_" << i;
		ka[i] = IloIntVar(env, 0, K, ostr.str().c_str());
		ostr.str("");
		IloExpr cardexpr(env);

		for (id j = 0; j < N; j++)
			cardexpr += IJ(xa, i, j);

		#ifdef DEBUG
		cout << (cardexpr <= K) << endl;
		#endif

		model.add(ka[i] == cardexpr);
		model.add(cardexpr <= K);
		cardexpr.end();
	}

	#ifdef DEBUG
	puts("");
	#endif

	/*
	// Alone variables

	IloIntVarArray na(env, N);
	IloIntVarArray a(env, N); // is i alone?

	for (id i = 0; i < N; i++) {

		ostr << "NA_" << i;
		na[i] = IloIntVar(env, 0, 1, ostr.str().c_str());
		ostr.str("");

		ostr << "A_" << i;
		a[i] = IloIntVar(env, 0, 1, ostr.str().c_str());
		ostr.str("");

		IloExpr aexpr(env);

		for (id j = 0; j < N; j++)
			if (i != j) aexpr += IJ(xa, i, j);

		#ifdef DEBUG
		cout << (na[i] <= lexpr) << endl;
		cout << (lexpr <= (N - 1) * na[i]) << endl;
		#endif

		model.add(na[i] <= aexpr);
		model.add(aexpr <= N * na[i]);
		model.add(a[i] == 1 - na[i]);
		aexpr.end();
	}

	// Leaders variables

	IloIntVarArray nla(env, N); // number of leaders in the cluster of i
	IloIntVarArray hla(env, N); // does the cluster of i have a leader?

	for (id i = 0; i < N; i++) {

		ostr << "NL_" << i;
		nla[i] = IloIntVar(env, 0, N, ostr.str().c_str());
		ostr.str("");

		ostr << "HL_" << i;
		hla[i] = IloIntVar(env, 0, 1, ostr.str().c_str());
		ostr.str("");

		IloExpr lexpr(env);
		lexpr += la[i];

		for (id j = 0; j < N; j++)
			if (i != j) lexpr += IJ(xa, i, j) * la[j];

		#ifdef DEBUG
		cout << (nla[i] == lexpr) << endl;
		cout << (nla[i] <= MAXLEADERS) << endl;
		#endif

		model.add(nla[i] == lexpr);
		model.add(nla[i] <= MAXLEADERS);
		model.add(hla[i] <= nla[i]);
		model.add(nla[i] <= N * hla[i]);
		lexpr.end();
	}
	*/

	// Create objective expression

	IloExpr objexpr(env);

	for (id i = 0; i < N; i++)
		for (id j = i; j < N; j++)
			objexpr += IJ(adj, i, j) * IJ(xa, i, j);

	#ifdef DEBUG
	cout << "Objective function:" << endl << objexpr << endl << endl;
	#endif

	model.add(IloMaximize(env, objexpr));
	objexpr.end();

	// Solve

	IloCplex cplex(model);
	struct timeval t1, t2;
	gettimeofday(&t1, NULL);

	#ifdef CSV
	cplex.setOut(env.getNullStream());
	#endif

	#ifndef PARALLEL
	cplex.setParam(IloCplex::Threads, 1);
	#endif

	#ifdef TOLERANCE
	cplex.setParam(IloCplex::EpGap, TOLERANCE);
	#endif

	try {
		if (!cplex.solve()) {
			#ifndef CSV
			env.out() << "Unable to find a solution" << endl;
			#endif
			exit(EXIT_FAILURE);
		}
	}
	catch (IloCplex::Exception e) {
		#ifndef CSV
		env.out() << "An exception occurred" << endl;
		#endif
		exit(EXIT_FAILURE);
	}

	gettimeofday(&t2, NULL);

	// Print solution

	#ifdef CSV
	printf("%f,%f\n", cplex.getObjValue(), (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec);
	#else
	puts("");
	printvarmatrix(xa, cplex, "X");
	printvarmatrix(fa, cplex, "F");
	//printvars(a, cplex);
	//printvars(nla, cplex);
	//printvars(hla, cplex);
	env.out() << "Optimal solution = " << cplex.getObjValue() << endl;
	printf("Clock elapsed time = %f\n", (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec);
	#endif

	// Free data structures

	free(adj);
	free(la);

	return 0;
}
