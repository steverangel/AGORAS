/*!
	AGORAS - an algorithm for approximating the k-medoids in clustered data
	Steve Rangel 2012
	Northwestern University CUCIS
	steverangel@u.northwestern.edu
*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <stdbool.h>
#include <string.h>
#include "agoras.h"

static const float GAMMA = 0.577215665;	// Euler-Mascheroni constant
static const float GC = 1.0;		// the coefficient for the sample size adaptive growth formula
static const float RC = 0.95;		// the sample size reduction factor

struct sample_item {			// struct for sample set objects
	float *data_item;		// pointer to a row in the data_set	
	struct sample_item *next;	// used in phase 2 for arranging objects in a ...
	struct sample_item *tail;	// mapping tree into a linked list of objects
	float cost;			// cost is the sum dist to all other objects in the list
};
typedef struct sample_item item;

static item *magic_stop = (item*)0xFFFFDEAD;

// return the index of the nearest neighbor using the provided distance function, object, and neighbors
int nn(float (*dist)(int, float *, float *), int dim, item obj, item *neighbor, int len) {
	int m,i;
	float d,tmp;
	d = (*dist)(dim,obj.data_item,neighbor[0].data_item);
	m = 0;
	for (i=1;i<len;i++) {
		tmp = (*dist)(dim,obj.data_item,neighbor[i].data_item);
		if (tmp<d) {
			d=tmp;
			m=i;
		}
	}
	return m;
}

// merge lists whose heads are objects in A with a new head that is an object in B, B's head is the nn of A's head
int merge_lists(item *A, item *B,float (*dist)(int, float *, float *), int dim, int size_sample) {
	int j;
	for(j=0;j<size_sample;j++) {		// for each sample item
		if (A[j].next) {		// if it has a mapping, or in the first sample set
			int m = nn(dist,dim,A[j],B,size_sample);	// m is the nn in B to A[j]
			B[m].tail->next=&A[j];	// add A[j] to the end of B[m]'s list
			B[m].tail=A[j].tail;	// Update B[m]'s end of list
		}
	}
	int list_count=0;			// reset the counter between each set mapping
	for(j=0;j<size_sample;j++) {		// count the items where is_mapped is true
		if(B[j].next) {			// more accurate than just counting the mappings
			list_count++;
		}
	}
	return list_count;
}

void **draw_samples(item **sample, int num_sets, int sample_size, float **data_set, int len) {
	int i,j,r;
	sample[0]=(item*)malloc(sample_size*sizeof(item)); // the first sample set
	for(j=0;j<sample_size;j++) {
		r = rand() % len;
		sample[0][j].data_item = data_set[r];
		sample[0][j].next = magic_stop;
		sample[0][j].tail = &sample[0][j];
		sample[0][j].cost = 0.0;
	}
	for(i=1;i<num_sets;i++) { //randomly draw the other n-1 sample sets from the data
		sample[i]=(item*)malloc(sample_size*sizeof(item));
		for(j=0;j<sample_size;j++) {
			r = rand() % len;		// len is the number of objects in the data_set
			sample[i][j].data_item = data_set[r];
			sample[i][j].next = NULL;
			sample[i][j].tail = &sample[i][j];
			sample[i][j].cost = 0.0;
		}
	}
}

// free the individual sample sets
void free_samples(item **samples, int num_sets) {
	int i;
	for (i=0;i<num_sets;i++)
		free(samples[i]);
}

// list heads are somewhere in the last sample set array, find and return them
item **get_list_heads(item **sample,int n, int s, int k) {
	int i,j;
	item **list_head = (item**)malloc(k*sizeof(item*));
	int count = 0;
	for (j=0;j<s;j++) {
		if (sample[n-1][j].next) {
			list_head[count++] = &sample[n-1][j];
		}
	}
	return list_head;
}

// find the most center object in the list
item* center_object(item *head,float (*dist)(int, float *, float *),int dim) {
	item *current = head;			// is the current object.
	while(current!=magic_stop) { 		// while there are objects in the list
		item *compare = current;	// point the compare object to the current object
		while(compare!=magic_stop) {
			float d =(*dist)(dim,current->data_item,compare->data_item);
			current->cost+=d;	// find the distance to all other objects ...
			compare->cost+=d;	// in a triangular computation of the distances ...
			compare=compare->next;	// to avoid recomputing the same distance
		}
		current = current->next;	// make current object the next object in the list
	}
	// use the distances to find the center
	item *center = head;
	current = head->next;
	while(current!=magic_stop) {
		if (current->cost<center->cost) {
			center = current;
		}
		current = current->next;
	}
	return center;
}

/*! Inputs:
	data_set	a dim x len matrix
	dim		dimension of data
	len		the number of data items
	k		number of clusters
	n		number of sample sets do draw	
 	s		sample size, if zero the expected value is used
	dist		a distance metric
	randseed	the random seed to use for sampling
 Returns:
	medoids		k pointers to the cluster centers in the data_set
 Note: 
	The membership function is not calculated here.
*/
float **seq_agoras(float **data_set, int dim, int len, int k, int n, int *s, int *nRestart, float (*dist)(int, float *, float *), unsigned int randseed)
{
	srand(randseed);	// seed the random number generator for sampling
	(*nRestart) = 0;	// the number of times the mapping phase restarts

	item **sample;
	sample = (item**)malloc(n*sizeof(item*));	// the sample sets

	float **medoids = (float **)malloc(k*sizeof(float*));	// these are pointers to the data_set

	if ((*s)==0)	// set s to E(W) if left to zero
		(*s)=(int)ceil((double)k*log((double)k)+GAMMA*(double)k); // E(W)

        //////// PHASE 1 -- Sampling and mapping --////////

	while(1) {	// loop until mapping successfully completes
		int i, list_count;
		draw_samples(sample,n,(*s),data_set,len);

		// find the nearest neighbors and map to them
		float d,t;
		for(i=0;i<(n-1);i++) {	// for each sample set, except for the last
			list_count=merge_lists(sample[i],sample[i+1],dist,dim,(*s));
			if (list_count<k) {	// there is no need to continue mapping
				break;	// break out of the for loop over the sample sets
			}
		} // end of for loop over sample sets
		if (list_count<k) { // increase sample size
			//printf("Restarting\n");
			(*s) = (*s)+(int)ceil((float)(*s)*(GC*((float)n-(float)i)/(float)n));
			(*nRestart)++;
			free_samples(sample,n);
		}
		else if(list_count==k) {
			break;	// break out of the the while loop -- continue to phase 2
		}
		else {	// reduce the sample size, it was too large
			(*s) =(int)ceil((float)(*s)*RC);
			(*nRestart)++;
			free_samples(sample,n);
		}
	} // end of mapping while loop

        //////// PHASE 2 -- Finding centers--//////// 

	item **list_head = get_list_heads(sample,n,(*s),k);

	// find the centermost object in each list and set it as a medoid
	int i;
	for (i=0;i<k;i++) {
		item *center = center_object(list_head[i],dist,dim);
		medoids[i] = center->data_item;
	}
	free_samples(sample,n);
	free(sample);
	free(list_head);
	return medoids;
}
