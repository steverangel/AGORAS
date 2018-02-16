/*
Author:	Esteban Rangel steverangel@gmail.com
Date:	6-19-2012
Description: This is an implementation of the PAM (Partitioning Around Medoids) algorithm.
*/
#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>

#define INLINE_UTIL
#include "cluster_utils.h"
#include "pam.h"

float swapmedoids(int nclusters, int nobjects, int ncoords, float **objects, int *objMembership, int *clusterMedoids, float(*distcalc)(int,float *, float *),float oldcost,float** clusterMedoidsDistMatrix) {

	int iter1, iter2, m, n, i, index;
	int *testMembership;
	int *tempMembership;
	float oCost;
	float existingCost;
	float configCost;
	float bestCost;
	float diffCost;
	float *tempDist;
	float *testDist;
	int swap_medoid, swap_object;

	/* allocate memory for temp Dist array */
	tempDist = (float *) malloc(nobjects * sizeof(float));
	assert(tempDist != NULL);

	/* allocate memory for test Dist array */
	testDist = (float *) malloc(nobjects * sizeof(float));
	assert(testDist != NULL);
	
	/* allocate memory for test membership array */
	testMembership = (int *) malloc(nobjects * sizeof(int));
	assert(testMembership != NULL);

	/* allocate memory for temp membership array */
	tempMembership = (int *) malloc(nobjects * sizeof(int));
	assert(tempMembership != NULL);
	
	bestCost = 0.0;
	for (iter1 = 0; iter1 < nclusters; iter1++) {		// for each medoid	
		for (iter2 = 0; iter2 < nobjects; iter2++) {	// swap medoid with each non-medoid
			configCost = 0.0;
			for (n = 0; n < nobjects; n++) {	// reset these values on each swap
				testMembership[n] = objMembership[n];
			}
			configCost = partition_swap(nclusters, nobjects, ncoords, objects, clusterMedoids, testMembership, distcalc, clusterMedoidsDistMatrix, iter1, iter2, testDist);
			if (configCost < bestCost) {		// if it is better save the config
				bestCost = configCost;
				swap_medoid = iter1; swap_object = iter2;
				// update membership array and dist array
				for (m = 0; m < nobjects; m++) {
					tempMembership[m] = testMembership[m];
					tempDist[m] = testDist[m];
				}
			}
		}
	}
	if (bestCost < 0.0) {	// if there was a better configuration found
		clusterMedoids[swap_medoid] = swap_object;
		for (m = 0; m < nobjects; m++) {
			objMembership[m] = tempMembership[m];
			clusterMedoidsDistMatrix[swap_medoid][m]=tempDist[m]; // update the row for the new medoid object
		}
	}
	free(testMembership);
	free(tempMembership);
	return bestCost;
}
	
int* seq_pam(float **objects, int numCoords, int numObjs, int numClusters, int *objMembership, float(*distcalc)(int,float *, float *), int maxiter, unsigned int randseed, float *newcost) {

	int i, j, r;
	int *clusterMedoids;
	float **clusterMedoidsDistMatrix;
	float objCost;
	float diffCost;

	srand(randseed);
	
	/* initialize object membership array */
	for (i = 0; i<numObjs; i++) objMembership[i] = 0;
	
	/* allocate space for cluster medoids dist matrix */
	// this is the distance from each object to the medoids
	clusterMedoidsDistMatrix = (float **)malloc(numClusters * sizeof(float*));
	assert(clusterMedoidsDistMatrix != NULL);
	// the size is numClusters by numObjs
	clusterMedoidsDistMatrix[0] = (float*)malloc(numClusters * numObjs * sizeof(float));
	assert(clusterMedoidsDistMatrix[0] != NULL);
	// fix the indicies for 2d refrence
	for (i=1; i<numClusters; i++)
		clusterMedoidsDistMatrix[i] = clusterMedoidsDistMatrix[i-1] + numObjs;
	
	/* allocate space for cluster medoids array */
	// holds indexes of objects that are treated as medoids
	clusterMedoids = (int*)malloc(numClusters * sizeof(int));
	assert(clusterMedoids != NULL);

	/* pick random data objects as initial set of medoids */
	for (i=0; i<numClusters; i++) {
		r = rand() % numObjs;
		clusterMedoids[i] = r;	// r is a random index of an object
	}
	
	float oldcost;
	(*newcost) = -1;
	// initial partition fills out clusterMedoidsDistMatrix and objMembership
	partition_init(numClusters, numObjs, numCoords, objects, clusterMedoids, objMembership, distcalc, clusterMedoidsDistMatrix);
	while ((*newcost)<0) {
		oldcost = (*newcost);
		(*newcost) = swapmedoids(numClusters, numObjs, numCoords, objects, objMembership, clusterMedoids, distcalc, oldcost, clusterMedoidsDistMatrix);
	}
	return(clusterMedoids);
}
