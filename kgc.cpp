#include "kgc.h"

int main(int argc, char *argv[]) {

	#ifdef DEBUG
	printf("# Vertices = %u\n# Edges = %u\nMaximum cardinality = %u\n\n", N, E, K);
	#endif

	// Allocate data structures

	value *adj = (value *)calloc(N * N, sizeof(value));

	// Read input file

	#define MAXLINE 1000
	char line[MAXLINE];
	FILE *f = fopen(argv[1], "r");

	for (id i = 0; i < N; i++) {
		fgets(line, MAXLINE, f);
		char *pch = line;
		if (*pch == '*') {
			pch++;
		}
		adj[i * N + i] = atof(pch);
	}

	for (id i = 0; i < E; i++) {
		id v1, v2;
		value w;
		fscanf(f, "%u %u %f", &v1, &v2, &w);
		adj[v1 * N + v2] = adj[v2 * N + v1] = w;
	}

	fclose(f);

	#ifdef DEBUG
	puts("Weights:");
	for (id i = 0; i < N; i++)
		printbuf(adj + i * N, N, NULL, "%+.4f");
	puts("");
	#endif

	// Free data structures

	free(adj);

	return 0;
}
