#ifndef INLINE_GOOD
# define INLINE_GOOD extern inline
#endif

INLINE_GOOD void stat_errors(int* true_membership, int* predicted_membership, int nobjs, int* TP,int* FP,int* TN,int* FN)
{
	(*TP) = 0;
	(*TN) = 0;
	(*FP) = 0;
	(*FN) = 0;
	int i,j;
	// for each pair in the predicted memberships
	float percent_complete = 0.0;
	//for (i=0;i<100;i++)
		//printf("=");
	//printf("\n=");
	for (i=0;i<nobjs;i++)
	{
		float cur_complete = floorf(((float)i/(float)nobjs)*100.0);
		//printf("%f\n",cur_complete);
		if (cur_complete>percent_complete) {
			percent_complete = cur_complete;
			//printf("=");
			//fflush(stdout);
		}
		for (j=i+1;j<nobjs;j++)
		{	// TP same in pred and same in true
			if      ((predicted_membership[i]==predicted_membership[j])&&(true_membership[i]==true_membership[j]))
				(*TP)++;
			else if ((predicted_membership[i]!=predicted_membership[j])&&(true_membership[i]!=true_membership[j]))
				(*TN)++;
			else if ((predicted_membership[i]==predicted_membership[j])&&(true_membership[i]!=true_membership[j]))
				(*FP)++;
			else if ((predicted_membership[i]!=predicted_membership[j])&&(true_membership[i]==true_membership[j]))
				(*FN)++;
		}
	}
	//printf("\n");
}

INLINE_GOOD float rand_index(int* true_membership, int* predicted_membership, int nobjs)
{
	float RI = 0.0;
	int TP, TN, FP, FN;
	stat_errors(true_membership,predicted_membership,nobjs,&TP,&FP,&TN,&FN);
	RI=(float)(TP+TN)/(float)(TP+FP+FN+TN);
	return RI;
}

INLINE_GOOD float f_measure(int* true_membership, int* predicted_membership, int nobjs, float beta)
{
	float f_b = 0.0;
	float P, R;
	int TP, TN, FP, FN;
	stat_errors(true_membership,predicted_membership,nobjs,&TP,&FP,&TN,&FN);
	P = (float)TP/(float)(TP+FP);
	R = (float)TP/(float)(TP+FN);
	f_b = (((beta*beta)+1)*P*R)/(((beta*beta)*P)+R);
	return f_b;
}

INLINE_GOOD float fowlkes_mallows(int* true_membership, int* predicted_membership, int nobjs) {
	int TP, TN, FP, FN;
	stat_errors(true_membership,predicted_membership,nobjs,&TP,&FP,&TN,&FN);
	return sqrt(((float)TP/(float)(TP+FP))*((float)TP/(float)(TP+FN)));
}
