#ifndef INLINE_DIST
# define INLINE_DIST extern inline
#endif
INLINE_DIST float DTW(int numVars, float *coord1, float *coord2) {
	int i,j;
	int lenC1 = 0;
	int lenC2 = 0;

	//count the number of non-zero pixels
	for (i=0;i<numVars;i++) {
		if (coord1[i]>0)
			lenC1++;
		if (coord2[i]>0)
			lenC2++;
	}

	// arrays for the non-zero indices
	int *coord1_idx = (int*)malloc(lenC1*sizeof(int));
	int *coord2_idx = (int*)malloc(lenC2*sizeof(int));

	lenC1 = 0;
	lenC2 = 0;
	for (i=0;i<numVars;i++) {
		if (coord1[i]>0)
			coord1_idx[lenC1++]=i;
		if (coord2[i]>0)
			coord2_idx[lenC2++]=i;
	}
	
	float **dist = (float**)malloc(lenC1*sizeof(float*));
	dist[0] = (float*)malloc(lenC1*lenC2*sizeof(float));
	for (i=1;i<lenC1;i++) {
		dist[i] = dist[i-1] + lenC2;
	}
	for (i=0;i<lenC1;i++) {
		for (j=0;j<lenC2;j++) {
			dist[i][j] = pow((float)coord1_idx[i]-(float)coord2_idx[j],2.0);
		}
	}
	float **cost = (float**)malloc(lenC1*sizeof(float*));
	cost[0] = (float*)malloc(lenC1*lenC2*sizeof(float));
	for (i=1;i<lenC1;i++) {
		cost[i] = cost[i-1] + lenC2;
	}
	cost[0][0] = dist[0][0];
	for (i=1;i<lenC1;i++) {
		cost[i][0] = dist[i][0]+cost[i-1][0];
	}
	for (i=1;i<lenC2;i++) {
		cost[0][i] = dist[0][i]+cost[0][i-1];
	}
	for (i=1;i<lenC1;i++) {
		for (j=1;j<lenC2;j++) {
			float temp;
			if (cost[i-1][j]<cost[i][j-1])
				temp = cost[i-1][j];
			else
				temp = cost[i][j-1];
			
			if (temp<cost[i-1][j-1])
				cost[i][j] = temp;
			else
				cost[i][j] = cost[i-1][j-1];
		}
	}
	float total = cost[lenC1-1][lenC2-1];
	free(dist[0]);
	free(cost[0]);
	free(dist);
	free(cost);
	//printf("DTW %f\n",total);
	return total;
}

INLINE_DIST float cosine_sim(int numVars, float *coord1, float *coord2) {
	float mag_coord1 = 0.0;
	float mag_coord2 = 0.0;
	float dot_prod = 0.0;
	int i;
	for (i=0;i<numVars;i++) {
		mag_coord1 += coord1[i]*coord1[i];
		mag_coord2 += coord2[i]*coord2[i];
		dot_prod += coord1[i]*coord2[i];
	}
	float cos_sim = dot_prod/(sqrt(mag_coord1)*sqrt(mag_coord2));
	return -cos_sim;
}
INLINE_DIST float jaccard(int numVars, float *coord1, float *coord2) {
	int i, vect_union, vect_intersection;
	vect_union = 0;
	vect_intersection = 0;
	for (i=0;i<numVars;i++) {
		if ((coord1[i]>0)||(coord2[i]>0)) {
			vect_union++;
		}
		
		if ((coord1[i]>0)&&(coord2[i]>0)) {
			vect_intersection++;
		}
	}
//	printf("%d,%d\n",vect_intersection,vect_union);
	float dist = -log(((float)vect_intersection/(float)vect_union))/log(2.0);
	assert(dist>=0);
//	printf("dist %f\n",dist);
	return dist;
}

INLINE_DIST float L1_norm(int numVars, float *coord1, float *coord2) {
	int j;
	float distance = 0.0;
	
	for (j=0; j<numVars; j++)
	{
		distance += fabsf(coord1[j] - coord2[j]);
	}
	return(distance);
}

INLINE_DIST float squared_euclid_dist(int numVars, float *coord1, float *coord2) {
	int j;
	float distance = 0.0;
	
	for (j=0; j<numVars; j++)
	{
		distance += (coord1[j] - coord2[j]) * (coord1[j] - coord2[j]);
	}
	assert(distance>=0.0);
	return(distance);
}

INLINE_DIST float euclid_dist(int numVars, float *coord1, float *coord2) {
	float distance = 0.0;
	distance = sqrt(squared_euclid_dist(numVars,coord1,coord2));
	assert(distance>=0.0);
	return(distance);
}

INLINE_DIST float simple_partition(int nclusters, int nobjects, int ncoords, float** objects, float** medoids, 
					int* membership, float(*distcalc)(int,float *, float *)) {
	float cost = 0.0;
	int i;
	for (i=0; i<nobjects; i++)
	{
		int index = 0;
		float oCost = 0.0;
		index = find_nearest_cluster(nclusters, ncoords, objects[i], 
					 medoids, &oCost, distcalc);
		membership[i] = index;
		cost +=oCost;
	//	assert(cost>=0.0);
	}
	return cost;
}
INLINE_DIST float simple_partition_v2(int nclusters, int nobjects, int ncoords, float** objects, int* medoids, 
					int* membership, float(*distcalc)(int,float *, float *))
{
	float cost = 0.0;
	int i;
	for (i=0; i<nobjects; i++)
	{
		int index = 0;
		float oCost = 0.0;
		float cCost[nclusters];
		index = find_nearest_cluster_v2(nclusters, ncoords, objects[i], objects, 
					 medoids, cCost, &oCost, distcalc);
		membership[i] = index;
		cost +=oCost;
		assert(cost>=0.0);
	}
	return cost;
}
