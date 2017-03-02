#include "stdio.h"
#include <stdlib.h>
#include "grads_names.h"
 
writeb3d(int *nzp, int *nxp, int *nyp, float *array, char *fname)
{
  
         FILE *f1;
	 int i,j,k,point,icheck,size;
     
     
     size=(*nxp)*(*nyp)*(*nzp);
     
/* let's open the binary file */
   f1=fopen(fname,"ab");
   if (f1==NULL) 
   {
   printf("\n error opening file \n");
   exit(-1);
 }
     fwrite(array,sizeof(float),size,f1);
     fclose(f1);
     return;
   }
