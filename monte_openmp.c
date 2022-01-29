#include "omp.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define RADIUS 1000

int main(int argc, char **argv) {
	long n, i, count = 0;
	double x, y, pi, startTime, endTime;
	n = strtol(argv[1], NULL, 10);

	srandom(time(NULL));
	startTime = omp_get_wtime();
	#pragma omp parallel for private(i,x,y) shared(n) reduction(+:count)
	for (i = 0; i < 1000 * n; i++) {
		x = (double)(random() % (RADIUS*1000))/1000;
		y = (double)(random() % (RADIUS*1000))/1000;
		if (x*x + y*y < RADIUS*RADIUS) {
			count++;
		}
	}

	pi = (double) 4 * count / (1000 * n);
	endTime = omp_get_wtime();
	printf("%lf\n", endTime - startTime);
	return EXIT_SUCCESS;
}
