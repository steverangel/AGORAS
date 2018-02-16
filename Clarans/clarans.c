#include <stdlib.h>
#include <unistd.h>
#include <assert.h>
#include <stdio.h>

#define INLINE_UTIL
#include "cluster_utils.h"
#include "clarans.h"

// Returns bestCost: the cost of partitioning the best configuration found.
// Modifies: objMembership, currentClusterSize, currentClusterCost, currentClusterMedoids
float searchNeighborhood(int nclusters, int nobjects, int ncoords,
		 float **objects, int *objMembership, int *currentClusterSize,
		 float *currentClusterCost, int *currentMedoids,
		 float(*distcalc)(int,float *, float *), int maxneighbor, float ** clusterMedoidsDistMatrix)
{
	int neighborsChecked, i, j;
	int *randneighborMembership;
	float randneighborCost;
	float currentCost;
	float diffCost;
	int randMedoid;
	int randObject;
	
	/* allocate memory for randneighbor Membership array */
	randneighborMembership = (int *) malloc(nobjects * sizeof(int));
	assert(randneighborMembership != NULL);

	/* allocate memory for test Dist array */
	float* testDist = (float *) malloc(nobjects * sizeof(float));
	assert(testDist != NULL);

	currentCost = 0.0;

	neighborsChecked = 1;
	while(neighborsChecked<maxneighbor)
	{
		for (i=0;i<nobjects;i++) { // reset the membership
			randneighborMembership[i]=objMembership[i];
		}

		// Consider a random neighbor randneighborMedoids of currentClusterMedoids //

		randMedoid = rand() % nclusters; 	// random medoid of current
		randObject = rand() % nobjects;		// random object in the data
				
		randneighborCost = partition_swap(nclusters, nobjects, ncoords, objects, currentMedoids, 
					randneighborMembership, distcalc,clusterMedoidsDistMatrix,randMedoid,randObject,testDist);

		// If randneighborMedoids has a lower cost, set currentClusterMedoids to randneighborMedoids
		if (randneighborCost < 0.0)
		{
			neighborsChecked = 0;
			//currentCost = randneighborCost;
			currentMedoids[randMedoid]=randObject;
			for (i = 0; i < nobjects; i++) {	// update membership array
				objMembership[i] = randneighborMembership[i];
				clusterMedoidsDistMatrix[randMedoid][i]=testDist[i];
			}
			
		}
		else 
			neighborsChecked++;
	}
	free(randneighborMembership);
	
	return currentCost;
}

	
int *seq_clarans(float **objects,
			  int numCoords,
			  int numObjs,
			  int numClusters,
			  int *objMembership,
			  float(*distcalc)(int,float *, float *),
			  int maxneighbor,
			  int numlocal,
			  float *mincost,
			  unsigned int randseed)
{
	int i, j, index, local;
	int *newClusterSize;
	int *currentClusterMedoids;
	int* bestClusterMedoids;
	float *newClusterCost;
	float objCost;
	float diffCost;
	float currentcost;
	int *currentMembership;

	srand(randseed);
	
	/* initialize object membership array */
	for (i = 0; i<numObjs; i++) 
		objMembership[i] = -1;
	
	/* allocate space for currentClusterMedoids array */
	currentClusterMedoids = (int*)malloc(numClusters * sizeof(int));
	assert(currentClusterMedoids != NULL);

	/* allocate space for bestClusterMedoids array */
	bestClusterMedoids = (int*)malloc(numClusters * sizeof(int));
	assert(bestClusterMedoids != NULL);

	/* allocate space for currentMembership array and clusterSize array */
	currentMembership = (int *)calloc(numObjs,sizeof(int));
	assert(currentMembership != NULL);
	
	/* allocate space for newClusters array and clusterSize array */
	newClusterSize = (int *)calloc(numClusters,sizeof(int));
	assert(newClusterSize != NULL);
	
	newClusterCost = (float *)calloc(numClusters,sizeof(float));	
	assert(newClusterCost != NULL);
	
	local = 1;
	(*mincost) = 1.0; // set mincost to a large value
	float **clusterMedoidsDistMatrix;

	do
	{
		/* allocate space for cluster medoids dist matrix */
		clusterMedoidsDistMatrix = (float **)malloc(numClusters * sizeof(float*));
		assert(clusterMedoidsDistMatrix != NULL);
		clusterMedoidsDistMatrix[0] = (float*)malloc(numClusters * numObjs * sizeof(float));
		assert(clusterMedoidsDistMatrix[0] != NULL);
		for (i=1; i<numClusters; i++)
			clusterMedoidsDistMatrix[i] = clusterMedoidsDistMatrix[i-1] + numObjs;

		int r;
		/* pick numClusters random data objects as initial set of medoids */
		// an arbitrary node of G
		for (i=0; i<numClusters; i++)
		{
			r = rand() % numObjs;
			currentClusterMedoids[i] = r;
		}
		// Compute the cost and partitioning of the arbitrary node
		partition_init(numClusters, numObjs, numCoords, objects, currentClusterMedoids, 
					currentMembership, distcalc, clusterMedoidsDistMatrix);

		// Search the local neighborhood of the arbitrary node
		currentcost = searchNeighborhood(numClusters, numObjs, numCoords, objects, 
			currentMembership, newClusterSize, newClusterCost, currentClusterMedoids, distcalc, maxneighbor, clusterMedoidsDistMatrix);
		if (currentcost<(*mincost))
		{
			//(*mincost) = currentcost;	// set mincost to currentcost
			// set bestClusterMedoids to currentClusterMedoids
			for (i=0; i<numClusters; i++)
			{
				bestClusterMedoids[i] = currentClusterMedoids[i];
			}
			for (i=0;i<numObjs;i++)	// copy the membership
				objMembership[i]=currentMembership[i];
		}
		local++;
		free(clusterMedoidsDistMatrix[0]);
	} while (local < numlocal);
	
	free(clusterMedoidsDistMatrix);
	free(newClusterSize);
	free(newClusterCost);
	free(currentClusterMedoids);
	free(currentMembership);
	
	return(bestClusterMedoids);
}
