#include <math.h>
#include <string.h>
#ifndef INLINE_UTIL
# define INLINE_UTIL extern inline
#endif

INLINE_UTIL int find_nearest_cluster(int numClusters, int numCoords, float *object, float **clusters, float *cost, float(*distcalc)(int,float *, float *)) {

	int closest, i;
	float dist;
	
	/* find the closest cluster */
	closest = 0;
	*cost = (*distcalc)(numCoords, object, clusters[0]);
	for (i=1; i<numClusters; i++) {
		dist = (*distcalc)(numCoords, object, clusters[i]);
		if (dist < *cost) {
			*cost = dist;
			closest = i;
		}
	}
	return(closest);
}

INLINE_UTIL float partition(int nclusters, int nobjects, int ncoords, float** objects, float** medoids, int* membership, int* clusterSize, float* clusterCost, float(*distcalc)(int,float *, float *)) {

	float cost = 0.0;
	int i;
	for (i=0;i<nclusters;i++)
		clusterSize[i] = 0;
	for (i=0; i<nobjects; i++) {
		int index = 0;
		float oCost = 0.0;
		index = find_nearest_cluster(nclusters, ncoords, objects[i], 
					 medoids, &oCost, distcalc);
		membership[i] = index;
		clusterSize[index]++;
		clusterCost[index] += oCost;
		cost +=oCost;
	}
	return cost;
}

INLINE_UTIL int find_nearest_cluster_v2(int numClusters, int numCoords, float *object, float **objects, int* clusters, float *cost, float *dist, float(*distcalc)(int,float *, float *)) {

	int closest, i;	
	/* find the closest cluster */
	closest = 0;
	*dist = (*distcalc)(numCoords, object, objects[clusters[0]]);
	cost[0] = *dist;
	for (i=1; i<numClusters; i++) {
		cost[i] = (*distcalc)(numCoords, object, objects[clusters[i]]);
		if (cost[i] < *dist) {
			*dist = cost[i];
			closest = i;
		}
	}
	return(closest);
}
// This function produces the initial partition. 
// membership and clusterMedoidsDistMatrix are altered after calling this function
// Note: medoid-to-medoid costs are computed to save the dist when they are swapped (i.e. considered as non-medoids)
INLINE_UTIL void partition_init(int nclusters, int nobjects, int ncoords, float** objects, int* medoids, int* membership, float(*distcalc)(int,float *, float *),float** clusterMedoidsDistMatrix) {
	int i;
	for (i=0; i<nobjects; i++) {
		int index = 0;
		float cCost[nclusters];	// the cost from objects[i] to all the medoids
		float oCost;
		index = find_nearest_cluster_v2(nclusters, ncoords, objects[i],objects, 
					 medoids, cCost, &oCost, distcalc);
		int q;
		for (q=0;q<nclusters;q++) {
			clusterMedoidsDistMatrix[q][i]=cCost[q];
		}
		membership[i] = index;
	}
}
// M is an index of a medoid. medoid[M] is an index of objects
// A is an index of an object that is swapping with medoid[M]
INLINE_UTIL float partition_swap(int nclusters, int nobjects, int ncoords, float** objects, int* medoids, int* membership, float(*distcalc)(int,float *, float *),float** clusterMedoidsDistMatrix, int M, int A, float* A_dist) {

	float cost = 0.0;
	int i,j;	

	for (i=0; i<nobjects; i++) {
		A_dist[i] = (*distcalc)(ncoords,objects[i],objects[A]);
		float doj2,tmp; // the nearest medoid to oi without M or A
		int doj2_i;
		doj2 = 1000000.0;
		if (membership[i]==M) {
			for (j=0;j<nclusters;j++) {
				if (j!=M) {
					tmp = clusterMedoidsDistMatrix[j][i];
					if (tmp<doj2) {
						doj2 = tmp;
						doj2_i = j;
					}
				}
			}
			if (A_dist[i]>=doj2) {
				cost += doj2 - clusterMedoidsDistMatrix[M][i];
				membership[i] = doj2_i;
			}
			else if(A_dist[i]<doj2) {
				cost += A_dist[i] - clusterMedoidsDistMatrix[M][i];
			}
		}
		else {
			if (A_dist[i]>=clusterMedoidsDistMatrix[membership[i]][i])
				cost += 0;
			else if (A_dist[i]<clusterMedoidsDistMatrix[membership[i]][i]) {
				cost += A_dist[i]-clusterMedoidsDistMatrix[membership[i]][i];
				membership[i] = M; // M will be replaced by A later
			}
		}
	}
	return cost;
}
