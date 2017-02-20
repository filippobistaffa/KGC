#include "kgc.h"

#define IJ(BUF, I, J) ((BUF)[(I) * N + (J)])

template <typename type>
__attribute__((always_inline)) inline
void printvararray(type &ia, IloCplex &cplex) {

	for (id i = 0; i < ia.getSize(); i++) {
		try {
			if (abs(cplex.getValue(ia[i])) > EPSILON)
				cout << ia[i].getName() << " = " << cplex.getValue(ia[i]) << endl;
		}
		catch (IloException& e) { e.end(); }
	}
	puts("");
}

template <typename type>
__attribute__((always_inline)) inline
void printvarmatrix(type &ia, IloCplex &cplex, const char *name = NULL, const char *format = NULL) {

	if (name) printf("%s =\n", name);

	for (id i = 0; i < N; i++) {
		printf("[ ");
		for (id j = 0; j < N; j++) {
			try {
				if (format) {
					printf(format, abs(cplex.getValue(IJ(ia, i, j))) > EPSILON ?
						       cplex.getValue(IJ(ia, i, j)) : 0);
					printf(" ");
				} else std::cout << (abs(cplex.getValue(IJ(ia, i, j))) > EPSILON ?
						     cplex.getValue(IJ(ia, i, j)) : 0) << " ";
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

bool checksimmetry(IloIntVarArray &xa, IloCplex &cplex) {

	for (id i = 0; i < N; i++)
		for (id j = 0; j < N; j++)
			if (cplex.getValue(IJ(xa, i, j)) != cplex.getValue(IJ(xa, j, i)))
				return false;

	return true;
}

bool checktransitivity(IloIntVarArray &xa, IloCplex &cplex) {

	for (id i = 0; i < N; i++)
		for (id j = 0; j < N; j++)
			for (id k = 0; k < N; k++) {
				if (!(cplex.getValue(IJ(xa, i, j)) +
				      cplex.getValue(IJ(xa, j, k)) -
				      2 * cplex.getValue(IJ(xa, i, k)) <= 1))
					return false;
				if (!(cplex.getValue(IJ(xa, i, k)) +
				      cplex.getValue(IJ(xa, i, j)) -
				      2 * cplex.getValue(IJ(xa, j, k)) <= 1))
					return false;
				if (!(cplex.getValue(IJ(xa, j, k)) +
				      cplex.getValue(IJ(xa, i, k)) -
				      2 * cplex.getValue(IJ(xa, i, j)) <= 1))
					return false;
			}

	return true;
}

bool checkflow(IloIntVarArray &xa, IloIntVarArray &sfa, IloFloatVarArray &fa, IloCplex cplex, const uint8_t *adj) {

	for (id i = 0; i < N; i++) {

		value f = cplex.getValue(sfa[i]);

		for (id j = 0; j < N; j++) {

			if (!(cplex.getValue(IJ(fa, i, j)) <= IJ(adj, i, j) * cplex.getValue(IJ(xa, i, j)) * K)) {
				//printf("F_%u,%u = %f </= %f\n", i, j, cplex.getValue(IJ(fa, i, j)),
				//					IJ(adj, i, j) * cplex.getValue(IJ(xa, i, j) * K));
				return false;
			}

			if (i != j && IJ(adj, i, j))
				f += cplex.getValue(IJ(fa, j, i)) - cplex.getValue(IJ(fa, i, j));
		}

		if (f != 1) return false;
	}

	return true;
}

int main(int argc, char *argv[]) {

	#ifndef CSV
	printf("# Vertices = %u\n# Edges = %u\nMaximum cardinality = %u\n\n", N, E, K);
	#endif

	// Allocate data structures

	value *wm = (value *)calloc(N * N, sizeof(value));
	uint8_t *adj = (uint8_t *)calloc(N * N, sizeof(uint8_t));
	uint8_t *la = (uint8_t *)calloc(N, sizeof(uint8_t));

	// Read input file

	#define MAXLINE 1000
	char line[MAXLINE];
	FILE *f = fopen(argv[1], "r");

	#ifdef DOT
	puts("graph G {");
	#endif

	for (id i = 0; i < N; i++) {
		fgets(line, MAXLINE, f);
		char *pch = line;
		if (*pch == '*') {
			pch++;
			la[i] = 1;
		}
		IJ(wm, i, i) = atof(pch);
		#ifdef DOT
		printf("\t%u -- %u [label=\"%f\"];\n", i, i, IJ(wm, i, i));
		#endif
	}

	for (id i = 0; i < E; i++) {
		id v1, v2;
		value w;
		fscanf(f, "%u %u %f", &v1, &v2, &w);
		IJ(wm, v1, v2) = IJ(wm, v2, v1) = w;
		IJ(adj, v1, v2) = IJ(adj, v2, v1) = 1;
		#ifdef DOT
		printf("\t%u -- %u [label=\"%f\"];\n", v1, v2, w);
		#endif
	}

	#ifdef DOT
	puts("}\n");
	#endif
	fclose(f);

	#ifdef DEBUG
	puts("Weights:");
	for (id i = 0; i < N; i++)
		printbuf(wm + i * N, N, NULL, "%+08.4f");
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
			//IJ(xa, i, j) = IloIntVar(env, IJ(xk, i, j), IJ(xk, i, j), ostr.str().c_str());
			ostr.str("");
		}

	// Reflexivity constraints

	for (id i = 0; i < N; i++)
		model.add(IJ(xa, i, i) == 1);

	// Simmetry constraints

	for (id i = 0; i < N; i++)
		for (id j = 0; j < N; j++) {
			model.add(xa[i * N + j] == xa[j * N + i]);
			model.add(xa[j * N + i] == xa[i * N + j]);
		}

	// Transitivity constraints

	for (id i = 0; i < N; i++)
		for (id j = 0; j < N; j++)
			for (id k = 0; k < N; k++) {
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
		for (id j = 0; j < N; j++) cardexpr += IJ(xa, i, j);
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

	IloIntVarArray kia(env, N); // cardinality-up-to-i variables
	IloIntVarArray sfa(env, N); // source flow variables

	for (id i = 0; i < N; i++) {
		ostr << "KI_" << i;
		kia[i] = IloIntVar(env, 0, K, ostr.str().c_str());
		ostr.str("");
		ostr << "SF_" << i;
		sfa[i] = IloIntVar(env, 0, K, ostr.str().c_str());
		ostr.str("");
		IloExpr cardexpr(env);
		for (id j = 0; j < i; j++) cardexpr += IJ(xa, i, j);
		model.add(kia[i] == cardexpr);
		model.add(IloIfThen(env, kia[i] == 0, sfa[i] == ka[i]));
		model.add(IloIfThen(env, kia[i] > 0, sfa[i] == 0));
		cardexpr.end();
	}

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

	// Connectivity constraints

	IloFloatVarArray fa(env, N * N); // flow variables

	for (id i = 0; i < N; i++)
		for (id j = 0; j < N; j++) {
			ostr << "F_" << i << "," << j;
			IJ(fa, i, j) = IloFloatVar(env, 0, K, ostr.str().c_str());
			ostr.str("");
		}

	#ifdef DEBUG
	puts("Flow constraints:");
	#endif

	for (id i = 0; i < N; i++) {
		IloExpr fexpr(env);
		fexpr += sfa[i];
		for (id j = 0; j < N; j++) {
			#ifdef DEBUG
			cout << (IJ(fa, i, j) <= IJ(adj, i, j) * IJ(xa, i, j) * K) << endl;
			#endif
			model.add((IJ(fa, i, j) <= IJ(adj, i, j) * IJ(xa, i, j) * K));
			if (i != j && IJ(adj, i, j)) fexpr += IJ(fa, j, i) - IJ(fa, i, j);
		}
		#ifdef DEBUG
		cout << (fexpr == 1) << endl;
		#endif
		model.add(fexpr == 1);
		fexpr.end();
	}

	#ifdef DEBUG
	puts("");
	#endif

	// Create objective expression

	IloExpr objexpr(env);

	for (id i = 0; i < N; i++)
		for (id j = i; j < N; j++)
			objexpr += IJ(wm, i, j) * IJ(xa, i, j);

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

	//cplex.setParam(IloCplex::EpInt, 0);

	#ifdef DEBUG
	cplex.setParam(IloCplex::MIPKappaStats, CPX_MIPKAPPA_FULL);
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

	#ifdef DEBUG
	assert(checksimmetry(xa, cplex));
	assert(checktransitivity(xa, cplex));
	assert(checkflow(xa, sfa, fa, cplex, adj));
	#endif

	gettimeofday(&t2, NULL);

	#ifdef DEBUG
	env.out() << "\nSolution status = " << cplex.getStatus() << endl;
	env.out() << "Max kappa value = " << cplex.getQuality(IloCplex::KappaMax) << endl;
	env.out() << "Percentage of numerically stable simplex bases = " << cplex.getQuality(IloCplex::KappaStable) << endl;
	env.out() << "Percentage of numerically unstable simplex bases = " << cplex.getQuality(IloCplex::KappaUnstable) << endl;
	env.out() << "Kappa attention value = " << cplex.getQuality(IloCplex::KappaAttention) << endl;
	#endif

	// Print solution

	#ifdef CSV
	printf("%f,%f\n", cplex.getObjValue(), (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec);
	#else
	puts("");
	printvarmatrix(xa, cplex, "X");
	printvarmatrix(fa, cplex, "F");
	//printvararray(ka, cplex);
	//printvararray(kia, cplex);
	//printvararray(kiba, cplex);
	printvararray(sfa, cplex);
	//printvars(nla, cplex);
	//printvars(hla, cplex);
	env.out() << "Optimal solution = " << cplex.getObjValue() << endl;
	printf("Clock elapsed time = %f\n", (double)(t2.tv_usec - t1.tv_usec) / 1e6 + t2.tv_sec - t1.tv_sec);
	#endif

	// Free data structures

	free(adj);
	free(wm);
	free(la);

	return 0;
}
