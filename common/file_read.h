#ifndef INLINE_FILE_READ
# define INLINE_FILE_READ extern inline
#endif
#define MAX_CHAR_PER_LINE 128
INLINE_FILE_READ float** file_read(
		char *filename,		/* input file name */
                int  *numObjs,		/* no. data objects (local) */
                int  *numCoords,	/* no. coordinates */
		int hasClass,		/* last column has a class column */
		int hasID,		// first column has ID
		int **objectClass)	// hold the class of the returned objects
	{
	float **objects;		// returned objects that were read
	int i, j, len;
	FILE *infile;
	char *line, *ret;
	int   lineLen;
	ssize_t numBytesRead;

	if ((infile = fopen(filename, "r")) == NULL) 
	{
		fprintf(stderr, "Error: no such file (%s)\n", filename);
		return NULL;
	}
	
	/* first find the number of objects */
	lineLen = MAX_CHAR_PER_LINE;
	line = (char*) malloc(lineLen);
	assert(line != NULL);
	
	(*numObjs) = 0;
	while (fgets(line, lineLen, infile) != NULL) 
	{
		/* check each line to find the max line length */ 
		while (strlen(line) == lineLen-1) 
		{
			/* this line read is not complete */
			len = strlen(line);
			fseek(infile, -len, SEEK_CUR);
			
			/* increase lineLen */
			lineLen += MAX_CHAR_PER_LINE;
			line = (char*) realloc(line, lineLen);
			assert(line != NULL);
			
			ret = fgets(line, lineLen, infile);
			assert(ret != NULL);
		}
		if (strtok(line, " \t\n") != 0)
			(*numObjs)++;
	}
	rewind(infile);
	
	/* find the no. objects of each object */
	(*numCoords) = 0;
	while (fgets(line, lineLen, infile) != NULL) 
	{
		if (hasID!=1)
			(*numCoords) = 1;			
		if (strtok(line, " ,\t\n") != 0) 
		{
			//ignore the id (first coordiinate): numCoords = 1;
			while (strtok(NULL, " ,\t\n") != NULL) (*numCoords)++;
			break; /* this makes read from 1st object */
		}
	}
	rewind(infile);

	(*numCoords)-=hasClass; // Default hasClass is 0
	
	// allocate space for objects[][] and read all objects
	len = (*numObjs) * (*numCoords);
	objects    = (float**)malloc((*numObjs) * sizeof(float*));
	assert(objects != NULL);
	objects[0] = (float*) malloc(len * sizeof(float));
	assert(objects[0] != NULL);
	for (i=1; i<(*numObjs); i++)
		objects[i] = objects[i-1] + (*numCoords);
	i = 0;

	if (hasClass==1) /* allocate memory for class array */
	{
		(*objectClass) = (int *) malloc((*numObjs) * sizeof(int));
		assert(objectClass != NULL);		
	}

	// read all objects 
	while (fgets(line, lineLen, infile) != NULL) 
	{
		if (hasID==1) {
			strtok(line, " ,\t\n");
			j=0;
		}
		else {
			objects[i][0] = atof(strtok(line, " ,\t\n"));
			j=1;
		}
		for (j; j<(*numCoords); j++)
			objects[i][j] = atof(strtok(NULL, " ,\t\n"));
		if (hasClass==1)
			(*objectClass)[i] = atoi(strtok(NULL, " ,\t\n"));		
		i++;
	}
	
	fclose(infile);
	free(line);

    return objects;
}

INLINE_FILE_READ float** mnist_read(int **class_label) {
	FILE *fp1 = fopen("../datasets/MNIST/t10k-images.idx3-ubyte", "rb");
	FILE *fp2 = fopen("../datasets/MNIST/t10k-labels.idx1-ubyte", "rb");

	if(!fp1 || !fp2) {
		printf("Error opening file.\n");
		return 0; 
	}

	int magic_number;
	int num_images;
	int num_rows;
	int num_cols;

	fread(&magic_number,sizeof(magic_number),1,fp1);
	fread(&num_images,sizeof(num_images),1,fp1);
	fread(&num_rows,sizeof(num_rows),1,fp1);
	fread(&num_cols,sizeof(num_cols),1,fp1);

	magic_number = ntohl(magic_number);
	num_images = ntohl(num_images);
	num_rows = ntohl(num_rows);
	num_cols = ntohl(num_cols);

	unsigned char **image;
	image = (unsigned char**)malloc(num_images*sizeof(unsigned char*));
	image[0] = (unsigned char*)malloc(num_images*num_rows*num_cols*sizeof(unsigned char));

	int i,j;
	for (i=1;i<num_images;i++)
		image[i] = image[i-1] + num_rows*num_cols;

	fread(image[0],sizeof(unsigned char),num_images*num_rows*num_cols,fp1);

	fread(&magic_number,sizeof(magic_number),1,fp2);
	fread(&num_images,sizeof(num_images),1,fp2);

	magic_number = ntohl(magic_number);
	num_images = ntohl(num_images);

	unsigned char *label;
	label = (unsigned char*)malloc(num_images*sizeof(unsigned char));
	fread(label,sizeof(unsigned char),num_images,fp2);

	fclose(fp1);
	fclose(fp2);

	float **data = (float**)malloc(num_images*sizeof(float*));
	data[0] = (float*)malloc(num_images*num_rows*num_cols*sizeof(float));
	for (i=1;i<num_images;i++)
		data[i] = data[i-1] + num_rows*num_cols;

	(*class_label) = (int*)malloc(num_images*sizeof(int));

	for (i=0;i<num_images;i++) {
		for (j=0;j<num_rows*num_cols;j++)
			data[i][j] = (float)image[i][j];
		(*class_label)[i] = (int)label[i];
	}

	free(image[0]);
	free(image);
	free(label);

	return data;
}
