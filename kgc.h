#ifndef KGC_H_
#define KGC_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>

// CPLEX
#include <ilcplex/ilocplex.h>

#include "instance.h"
#include "params.h"
#include "macros.h"
#include "types.h"

using namespace std;

#ifdef CSV
#undef DEBUG
#undef DOT
#endif

#ifdef DEBUG
#include <assert.h>
#endif

#endif /* KGC_H_ */
