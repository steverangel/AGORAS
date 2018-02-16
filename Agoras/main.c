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
#define INLINE_UTIL
#include "cluster_utils.h"

#include "agoras.h"

double wtime(void);

static void usage(char *argv0) {
    char *help =
	"Usage: %s [switches] -i filename -n num_clusters\n"
	"       -i filename    	: file containing data to be clustered (csv file)\n"
	"       -n num_clusters	: number of clusters k (k must > 1)\n"
	"       -l num_samples  : number of samples sets\n"
	"       -r rand_seed  	: random Seed for initialization (default is time)\n"
	"       -c has_class  	: last column has the class (default 0)\n"
	"       -g goodness  	: calculate the goodness F1 measure based on class (default 0)\n"
	"       -m membership  	: calculate the membership function (default 0)\n"
	"       -f output files : generate output files for medoids and membership (default 0)\n"
	"       -v verbose	: detailed output\n";
    fprintf(stderr, help, argv0);
    exit(-1);
}					

int	main(int argc, char **argv) {
	int opt;
	extern char *optarg;
	extern int optind;
	int nSamples;
	int sampleSize;
	int nRestarts;
	char *filename;
	int numClusters, numCoords, numObjs;
	bool hasClass, hasID, calcGoodness, calcMembership, genOutputFiles, verbose;
	int *membership;
	int *true_class;
	float **objects;
	float **clusters;
	double	timing, clusterTiming, partTiming;
	unsigned int randseed;
	float cost, goodness;
	char medoids_file_name [256];
	char membership_file_name[256];
	
	/* default values for input arguments */
	nSamples = -1;		//later set to (int)ceil((double)k*log((double)k)+0.577215665*(double)k);
	numClusters = 0;
	filename = NULL;
	randseed = time(NULL) + getpid();	// by default use a time seed
	hasClass = false;
	hasID = false;
	calcGoodness = false;
	calcMembership = false;
	genOutputFiles = false;
	verbose = false;
	sampleSize = 0; // not an argument
	
	while ((opt=getopt(argc,argv,"i:n:l:r:c:d:g:m:f:v:"))!= EOF)
	{
		switch (opt)
		{
			case 'i': 
				filename = optarg;
				break;
			case 'n':
				numClusters = atoi(optarg);
				break;
			case 'l':
				nSamples = atoi(optarg);
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

	if (nSamples==-1)
		nSamples = (int)ceil((double)numClusters*log((double)numClusters)+0.577215665*(double)numClusters);
	
	if (filename == 0 || numClusters <= 1) 
		usage(argv[0]);
	
	objects = file_read(filename, &numObjs, &numCoords, hasClass, hasID, &true_class);
	if (objects == NULL) exit(1);

	int repeat;
	for (repeat=0;repeat<1;repeat++)
	{
	sampleSize = 0;
	//randseed = rand() % 5000;
	clusterTiming = wtime();
	// ************	START AGORAS	************
	// use the squared_euclid_dist
	clusters = seq_agoras(objects,numCoords,numObjs,numClusters,nSamples,&sampleSize,&nRestarts,squared_euclid_dist,randseed);
	// ************	END AGORAS	************
	clusterTiming = wtime() - clusterTiming;

	/* allocate memory for membership array */
	membership = (int *) malloc(numObjs * sizeof(int));
	assert(membership != NULL);

	sprintf(medoids_file_name,"%d_medoids.csv",randseed);
	sprintf(membership_file_name,"%d_membership.csv",randseed);

	cost = 0.0;
	partTiming = 0.0;
	if (calcMembership) {
		partTiming = wtime();
		// find cost with the euclid_dist
		cost = simple_partition(numClusters,numObjs,numCoords,objects,clusters,membership, squared_euclid_dist);
		partTiming = wtime() - partTiming + clusterTiming;
	}
	else {
		calcGoodness = false;	// goodness needs membership function
	}
	if (!hasClass)	// cannot calc the goodness without ground truth
		calcGoodness = false;
	
	if (verbose) 
		printf("\n **** Regular AGORAS (sequential version) ****\n");

	goodness = -1.0;
	if (calcGoodness)
		goodness = f_measure(true_class,membership,numObjs,1.0);

	if (verbose) {
		printf("     Input file: %s\n", filename);
		printf("    Random seed: %d\n", randseed);
		printf("       Num Objs: %d\n", numObjs);
		printf("     Num Coords: %d\n", numCoords);
		printf("    numClusters: %d\n", numClusters);
		printf("   Compute time: %12.8f sec\n", clusterTiming);
		printf("   Compute time: %12.8f sec (including partitioning)\n", partTiming);
		if (calcMembership) 
		printf("           Cost: %lf\n",cost);
		if (calcGoodness) 
		printf("       Goodness: %12.8f\n",goodness);
		printf("    Sample Size: %d\n", sampleSize);
		printf("    Num Samples: %d\n", nSamples);
		printf("   Num Restarts: %d\n", nRestarts);
		if (genOutputFiles) {
		printf("   Medoids file: %s\n",medoids_file_name);
			if (calcMembership) 
		printf("Membership file: %s\n",membership_file_name);
		}
	}
	else { // non-verbose
		printf("%d,%s,%d,%d,%d,%lf,%lf,%lf,%lf,%d,%d,%d\n",
		randseed,filename,numObjs,numCoords,numClusters,clusterTiming,partTiming,cost,goodness,sampleSize,nSamples,nRestarts);
	}

	if (genOutputFiles) { // output files for medoids and membership
		FILE *medoids_file = fopen(medoids_file_name,"w");
		int i,j;
		for(i=0;i<numClusters;i++) {
			fprintf(medoids_file,"%d,",i);
			for(j=0;j<numCoords-1;j++)
				fprintf(medoids_file,"%lf,",clusters[i][j]);
			fprintf(medoids_file,"%lf\n",clusters[i][numCoords-1]);
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


        free(clusters);
	free(membership);
}
	free(objects[0]);
        free(objects);

	if (hasClass)
		free(true_class);
	
	return(0);
}
