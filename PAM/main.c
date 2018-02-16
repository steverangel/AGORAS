#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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

#include "pam.h"

double wtime(void);

static void usage(char *argv0) {
    char *help =
	"Usage: %s [switches] -i filename -n num_clusters\n"
	"       -i filename    	: file containing data to be clustered\n"
	"       -n num_clusters	: number of clusters (K must > 1)\n"
	"       -m max_pass  	: Max number of passes (default 100)\n"
	"       -t threshold  	: Threshold for cost percent decrease (default 0.0)\n"
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
	int maxiter;
	char *filename;
	int numClusters, numCoords, numObjs, isBinaryFile, hasClass;
	int *membership;
	int *true_class;
	float **objects;
	int* clusters;
	double timing, clusterTiming;
	unsigned int randseed;
	int verbose;
	
	/* default values for input arguments */	
	numClusters = 0;
	filename = NULL;
	maxiter = 100;
	randseed = time(NULL);	// by default use a time seed
	hasClass = 0;
	verbose = 0;
	
	while ((opt=getopt(argc,argv,"i:b:n:m:r:c:v:"))!= EOF)
	{
		switch (opt)
		{
			case 'i': 
				filename = optarg;
				break;
			case 'n':
				numClusters = atoi(optarg);
				break;
			case 'm':
				maxiter = atoi(optarg);
				break;
			case 'r':
				randseed = atoi(optarg);
				break;
			case 'c':
				hasClass = atoi(optarg);
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

	//objects = file_read(filename, &numObjs, &numCoords, hasClass, &true_class);
	objects = file_read(filename, &numObjs, &numCoords, hasClass, 1, &true_class);
	if (objects == NULL) exit(1);

	/* allocate memory for membership array */
	membership = (int *) malloc(numObjs * sizeof(int));
	assert(membership != NULL);

	float cost;

	timing = 0.0;
	clusterTiming = 0.0;	
	timing = wtime();
	clusterTiming = timing;
	
	// ************	START PAM	************
	clusters = seq_pam(objects, numCoords, numObjs, numClusters, membership, 
				L1_norm, maxiter, randseed, &cost);
	// ************	END PAM		************	
	timing = wtime();
	clusterTiming = timing - clusterTiming;

	int q;

	float part_cost;
	part_cost = 0.0;
	// find cost with the euclid_dist
	part_cost = simple_partition_v2(numClusters, numObjs, numCoords, objects, clusters, 
						membership, L1_norm);

	// ************	OUTPUT RESULTS	************
	if (verbose == 1) printf("\n **** Regular PAM (sequential version) ****\n");
	int i,j,n;
	float goodness = -1.0;
	if (hasClass==1) 
	//	goodness = f_measure(true_class,membership,numObjs,1.0);
		goodness = 0;//fowlkes_mallows(true_class,membership,numObjs);
	
	if (verbose == 1)
	{
		printf("Random Seed:    %d\n", randseed);
		printf("Input file:     %s\n", filename);
		printf("numObjs:        %d\n", numObjs);
		printf("numCoords:      %d\n", numCoords);
		printf("numClusters:    %d\n", numClusters);
		printf("Comp timing:    %12.8f sec\n", clusterTiming);
		printf("Partition Cost: %lf\n",part_cost);
		if (hasClass==1) printf("Goodness: ,%12.8f\n",goodness);
		printf("maxiter:        %d\n",maxiter);
		printf("Medoids:\n");
		printf("ID\t");
		for (i=0;i<numCoords;i++)
			printf("Coord %d  \t",i);
		if (hasClass == 1)
			printf("True Class\n");
		else
			printf("\n");	
		for(i=0;i<numClusters;i++)
		{
			printf("%d\t",i);
			for(j=0;j<numCoords;j++)
				printf("%lf\t",objects[clusters[i]][j]);
			printf("\n");
		}
	}
	else // non-verbose
	{
		printf("%d,%s,%d,%d,%d,%lf,%lf,%lf,%d\n",
			randseed,filename,numObjs,numCoords,numClusters,clusterTiming,part_cost,goodness,maxiter);
	}

	free(objects[0]);
        free(objects);
	free(membership);
        //free(clusters[0]);
        free(clusters);
	if (hasClass == 1) free(true_class);

	return(0);
}	
