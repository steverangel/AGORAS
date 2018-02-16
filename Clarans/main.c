#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <sys/types.h>
#include <unistd.h>
#include <fcntl.h>
#include <assert.h>
#include <math.h>

#define INLINE_FILE_READ
#include "file_read.h"
#define INLINE_DIST
#include "distances.h"
#define INLINE_GOOD
#include "goodness.h"

#include "clarans.h"

double wtime(void);

static void usage(char *argv0) {
    char *help =
	"Usage: %s [switches] -i filename -n num_clusters\n"
	"       -i filename    	: file containing data to be clustered\n"
	"       -b             	: input file is in binary format (default no)\n"
	"       -n num_clusters	: number of clusters (K must > 1)\n"
	"       -l numlocal  	: numlocal\n"
	"       -x maxneighbor   : maxneighbor\n"
	"       -r rand_seed  	: Random Seed for initialization (default is time)\n"
	"       -c has_class  	: Last column has the class (default no)\n"
	"       -v verbose	: detailed output\n";
    fprintf(stderr, help, argv0);
    exit(-1);
}					

int	main(int argc, char **argv)
{
	int opt;
	extern char *optarg;
	extern int optind;
	int maxneighbor;
	int numlocal;
	char *filename;
	int numClusters, numCoords, numObjs;
	bool hasClass, hasID, calcGoodness, calcMembership, genOutputFiles, verbose;
	int *membership;
	int *true_class;
	float **objects;
	int* clusters;
	double	timing, clusterTiming;
	unsigned int randseed;	
	char medoids_file_name [256];
	char membership_file_name[256];
	
	/* default values for input arguments */
	float percent_maxneighbor = 0.0125;
	numlocal = 2;
	numClusters = 0;
	filename = NULL;
	randseed = time(NULL)+getpid();	// by default use a time seed
	hasClass = false;
	hasID = false;
	calcGoodness = false;
	calcMembership = false;
	genOutputFiles = false;
	verbose = false;

	while ((opt=getopt(argc,argv,"i:b:n:l:x:r:c:d:g:m:f:v:"))!= EOF)
	{
		switch (opt)
		{
			case 'i': 
				filename = optarg;
				break;
			case 'n':
				numClusters = atoi(optarg);
				break;
			case 'x':
				percent_maxneighbor = atof(optarg);
				break;
			case 'l':
				numlocal = atoi(optarg);
				break;
			case 'r':
				randseed = atoi(optarg);
				break;
			case 'c':
				hasClass = atoi(optarg);
				break;
			case 'd':
				hasID = atoi(optarg);
				break;
			case 'g':
				calcGoodness = atoi(optarg);
				break;
			case 'm':
				calcMembership = atoi(optarg);
				break;
			case 'f':
				genOutputFiles = atoi(optarg);
				break;
			case 'v':
				verbose = atoi(optarg);
				break;
			case '?':
				usage(argv[0]);
				break;
			default:
				usage(argv[0]);
				break;
		}
	}
	
	if (filename == 0 || numClusters <= 1) usage(argv[0]);
	
	objects = file_read(filename, &numObjs, &numCoords, hasClass, hasID, &true_class);

	if (objects == NULL) exit(1);

	maxneighbor = (int)ceil(percent_maxneighbor*(float)numObjs);
	
	/* allocate memory for membership array */
	membership = (int *) malloc(numObjs * sizeof(int));
	assert(membership != NULL);

	sprintf(medoids_file_name,"%d_medoids.csv",randseed);
	sprintf(membership_file_name,"%d_membership.csv",randseed);

	float cost = 0.0;

	timing = 0.0;
	clusterTiming = 0.0;
	timing = wtime();
	clusterTiming = timing;
	// ************	START CLARANS	************
	clusters = seq_clarans(objects, numCoords, numObjs, numClusters, membership, L1_norm, maxneighbor, numlocal, &cost, randseed);
	// ************	END CLARANS	************
	timing = wtime();
	clusterTiming = timing - clusterTiming;

	// find cost with the euclid_dist
	float part_cost = 0.0;
	if (calcMembership) {
		part_cost = simple_partition_v2(numClusters, numObjs, numCoords, objects, clusters, membership, L1_norm);
	}
	else {
		calcGoodness = false;
	}
	if (!hasClass)
		calcGoodness = false;

	if (verbose) 
		printf("\n **** Regular CLARANS (sequential version) ****\n");
	int i,j,n;
	float goodness = -1.0;
	if (calcGoodness)
		//goodness = f_measure(true_class,membership,numObjs,1.0);
		goodness = fowlkes_mallows(true_class,membership,numObjs);
	
	if (verbose)
	{
		printf("Random Seed:    %d\n", randseed);
		printf("Input file:     %s\n", filename);
		printf("numObjs:        %d\n", numObjs);
		printf("numCoords:      %d\n", numCoords);
		printf("numClusters:    %d\n", numClusters);
		printf("Computation timing:  %12.8f sec\n", clusterTiming);
		printf("Partition Cost: %lf\n",part_cost);
		if (hasClass==1) printf("Goodness: ,%12.8f\n",goodness);
		printf("Max Neighbor: %d\n",maxneighbor);
		printf("Num Local: %d\n",numlocal);
	}
	else // non-verbose
	{
		printf("%d,%s,%d,%d,%d,%lf,%lf,%lf,%d,%d\n",
			randseed,filename,numObjs,numCoords,numClusters,clusterTiming,part_cost,goodness,maxneighbor,numlocal);
	}
	
	if (genOutputFiles) { // output files for medoids and membership
		FILE *medoids_file = fopen(medoids_file_name,"w");
		int i,j;
		for(i=0;i<numClusters;i++) {
			fprintf(medoids_file,"%d,",i);
			for(j=0;j<numCoords-1;j++)
				fprintf(medoids_file,"%lf,",objects[clusters[i]][j]);
			fprintf(medoids_file,"%lf\n",objects[clusters[i]][numCoords-1]);
		}
		fclose(medoids_file);
		if (calcMembership) {
			FILE *membership_file = fopen(membership_file_name,"w");
			for (i=0;i<numObjs;i++) {
				fprintf(membership_file,"%d,%d\n",i,membership[i]);
			}
			fclose(membership_file);
		}
	}

	free(objects[0]);
        free(objects);	
	free(membership);
        free(clusters);
	if (hasClass == 1) free(true_class);
	
	return(0);
}

	
	
