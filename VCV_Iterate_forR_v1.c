#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <R.h>


//Created January 19, 2012
//Modified for R 14 Feb 2012

// memory allocation for double matrix

double **mem_alloc(int width, int length){
	int i;
	double **result;
	
	result=malloc(width*sizeof(double));
	if(!result){ printf("Memory allocation failure.\n"); exit(1);}
	
	for(i=0;i<width;i++){ 
		result[i]=malloc(length*sizeof(double));
		if(!result[i]){ printf("Memory allocation failure.\n"); exit(1);}
	}
	return result;
}
// free memory allocated to a double matrix

void free_matrix(double **matrix, int n, int m){
	
	int i;
	
	for(i=0;i<n;i++)
		free(matrix[i]);
	
	free(matrix);
	
} 
// memory allocation for double vector

double *mem_1d(int length){
	double *result;
	
	result=malloc(length*sizeof(double));
	if(!result){ printf("Memory allocation failure.\n"); exit(1);}
	
	return result;
}


double VCVP[1000][1000],VCVBar[1000];


void competitionInterval(double *mm, double *VCV, double *S, double *Psi, double *Eps, double *theta, int *Nspecies, int *Gmax)
{
	int i,j,k,G;

	double *VCVBar, *VCVP, *mp;

	VCVP=mem_1d(*Nspecies * *Nspecies);
	VCVBar=mem_1d(*Nspecies);
    mp=mem_1d(*Nspecies);
	
	for(G=0;G<*Gmax;G++)//evolutionary replicates
	{

		//first calculate covbars
		for(i=0;i<*Nspecies;i++)
		{

			VCVBar[i]=0;
			for(j=0;j<*Nspecies;j++)
			{

				VCVBar[i]=VCVBar[i]+VCV[*Nspecies*i+j];

			}

			VCVBar[i]=VCVBar[i]/(1.0 * *Nspecies);

		}
		//Now calculate VCV in next generation
		for(i=0;i < *Nspecies;i++)
		{

			for(j=0;j < *Nspecies;j++)
			{
				if(i==j)
				{
					VCVP[*Nspecies*i+j]=VCV[*Nspecies*i+j]-2.0*(*S + *Psi)*VCV[*Nspecies*i+j]+ *S *(VCVBar[i]+VCVBar[j])+ *Eps;
				}
				else
				{
					VCVP[*Nspecies*i+j]=VCV[*Nspecies*i+j]-2.0*(*S + *Psi)*VCV[*Nspecies*i+j]+ *S *(VCVBar[i]+VCVBar[j]);
				}
			}
		}
        
        //Now calculate VCV in next generation
        for(i=0;i < *Nspecies;i++)
		{
            mp[i]=mm[i]+ *Psi * (*theta - mm[i]);
        }
        
		//Now hand off prime array to starting array
		for(i=0;i < *Nspecies;i++)
		{
			for(j=0;j < *Nspecies;j++)
			{
				VCV[*Nspecies*i+j]=VCVP[*Nspecies*i+j];
			}
		}
        
        for(i=0;i < *Nspecies;i++)
		{
            mm[i]=mp[i];
        }
        
	}
	//end generational recursion

	free(VCVBar);
	free(VCVP);
    free(mp);
	
return;

}

