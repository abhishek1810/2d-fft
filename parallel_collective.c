#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/times.h>
#include <sys/time.h>
#include <time.h>
#include <mpi.h>

#define N 512

typedef struct {float r; float i;} complex;
static complex ctmp;

#define C_SWAP(a,b) {ctmp=(a);(a)=(b);(b)=ctmp;}

void print_array(complex *array)
{
  int row, col;

  if (N < 17)
  {
    printf("\nA =\n\t");
    for (row = 0; row < N; row++)
    {
      for (col = 0; col < N; col++)
      {
        printf("%5.2f%s", array[row*N+col].r, (col < N - 1) ? ", " : ";\n\t");
      }
    }
    printf("\n");
  }
}

void readFile(char fileName[15], complex *array){
    FILE *fp = fopen(fileName, "r");

    int i, j, result;

    for (i=0;i<N;i++){
        for (j=0;j<N;j++){
            result = fscanf(fp,"%g",&array[i*N+j].r);
            array[i*N+j].i = 0.00;
        }
    }

    fclose(fp);
}

void writeToFile(char fileName[15], complex *array){
    FILE *fp = fopen(fileName, "w");
    int i, j;
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++){
            fprintf(fp,"   %.7e",array[i*N+j].r);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

void transpose(complex *r, int n ) {
    for ( int i = 0; i < (n-1); i++ ) {
        for ( int j =i+1; j < n; j++ ) {
            C_SWAP(r[i*n+j],r[j*n+i]);
        }
    }
}

void c_fft1d(complex *r, int n, int isign)
{
   int     m,i,i1,j,k,i2,l,l1,l2;
   float   c1,c2,z;
   complex t, u;

   if (isign == 0) return;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n/2;i+=2) {
      if (i < j)
         C_SWAP(r[i], r[j]);
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* m = (int) log2((double)n); */
   for (i=n,m=0; i>1; m++,i/=2);

   /* Compute the FFT */
   c1 = -1.0;
   c2 =  0.0;
   l2 =  1;
   for (l=0;l<m;l++) {
      l1   = l2;
      l2 <<= 1;
      u.r = 1.0;
      u.i = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;

            /* t = u * r[i1] */
            t.r = u.r * r[i1].r - u.i * r[i1].i;
            t.i = u.r * r[i1].i + u.i * r[i1].r;

            /* r[i1] = r[i] - t */
            r[i1].r = r[i].r - t.r;
            r[i1].i = r[i].i - t.i;

            /* r[i] = r[i] + t */
            r[i].r += t.r;
            r[i].i += t.i;
         }
         z =  u.r * c1 - u.i * c2;

         u.i = u.r * c2 + u.i * c1;
         u.r = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (isign == -1) /* FWD FFT */
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for inverse transform */
   if (isign == 1) {       /* IFFT*/
      for (i=0;i<n;i++) {
         r[i].r /= n;
         r[i].i /= n;
      }
   }
}

int main(int argc, char **argv){
    /* Timing variables */
    struct timeval etstart, etstop;  /* Elapsed times using gettimeofday() */
    struct timezone tzdummy;
    clock_t etstart2, etstop2;  /* Elapsed times using times() */
    unsigned long long usecstart, usecstop;
    struct tms cputstart, cputstop;  /* CPU times for my processes */

    double starttime = 0.0;
    double endtime = 0.0;

    int i, j, index, div;
    int my_rank, p, source = 0, dest, x;
    int tag = 1;

    MPI_Status status;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);

    complex *A, *B, *C;

    A = (complex *)malloc(sizeof(complex)*N*N);
    B = (complex *)malloc(sizeof(complex)*N*N);
    C = (complex *)malloc(sizeof(complex)*N*N);

    MPI_Datatype complexStruct;
    int blockcounts[2] = { 1, 1 };
    MPI_Aint indices[2] = { 0, sizeof(float) };
    MPI_Datatype types[2] = { MPI_FLOAT, MPI_FLOAT };
    MPI_Type_struct( 2, blockcounts, indices, types, &complexStruct );

    MPI_Type_commit( &complexStruct );

    if ( my_rank == 0 ) {
        char inputA[15] = "sample/1_im1";
        char inputB[15] = "sample/1_im2";

        readFile(inputA, A);
        readFile(inputB, B);

            /* Start Clock */
        printf("\nStarting clock.\n");
        gettimeofday(&etstart, &tzdummy);
        etstart2 = times(&cputstart);
        starttime = MPI_Wtime();
    }

    /* Distribute Data for Step 1 */  
    div = N/p;
    complex *subA, *subB, *subC;
    subA = (complex *)malloc(sizeof(complex)*N*div);
    subB = (complex *)malloc(sizeof(complex)*N*div);
    subC = (complex *)malloc(sizeof(complex)*N*div);
    MPI_Scatter(&A[0], N*div, complexStruct, &subA[0], N*div, complexStruct, 0, MPI_COMM_WORLD);
    MPI_Scatter(&B[0], N*div, complexStruct, &subB[0], N*div, complexStruct, 0, MPI_COMM_WORLD);
    // printf("%d process has A[0] - %f\n", my_rank, subA[0].r);
    /*------------------------------*/


    /* Step 1.a - Perform 1D FFT on A and B */  
    for( i=0; i<div; i++){   
        c_fft1d(&subA[i*N], N, -1);
    }  
    for( i=0; i<div; i++){   
        c_fft1d(&subB[i*N], N, -1);
    }  
    MPI_Barrier(MPI_COMM_WORLD);
    /*------------------------------*/
    

    /* Step 1.b - Collect Data */  
    MPI_Gather(&subA[0], N*div, complexStruct, &A[0], N*div, complexStruct, 0, MPI_COMM_WORLD);
    MPI_Gather(&subB[0], N*div, complexStruct, &B[0], N*div, complexStruct, 0, MPI_COMM_WORLD);
    /*------------------------------*/


    /* Step 1.c - Transpose */  
    if (my_rank == 0) {
        transpose(A, N);
        transpose(B, N);
    }
    /*------------------------------*/


    /* Step 1.d - Distribute Data for Step 1 after transpose*/  
    MPI_Scatter(&A[0], N*div, complexStruct, &subA[0], N*div, complexStruct, 0, MPI_COMM_WORLD);
    MPI_Scatter(&B[0], N*div, complexStruct, &subB[0], N*div, complexStruct, 0, MPI_COMM_WORLD);
    /*------------------------------*/


    /* Step 1.e - Perform 1D FFT on A and B after transpose */  
    for( i=0; i<div; i++){   
        c_fft1d(&subA[i*N], N, -1);
    }  

    for( i=0; i<div; i++){   
        c_fft1d(&subB[i*N], N, -1);
    }
    /*------------------------------*/


    /* Step 2 - Matrix Mutipication */  
    for (i = 0; i < div ; i++) {
        for(j = 0; j < N ; j++) {
            index = i*N+j;
            subC[index].r = (subA[index].r * subB[index].r) - (subA[index].i * subB[index].i);
            subC[index].i = (subA[index].r * subB[index].i) + (subA[index].i * subB[index].r);
        }
    }
    /*------------------------------*/


    /* Step 3.a - Perform 1D FFT on C */  
    for( i=0; i<div; i++){   
        c_fft1d(&subC[i*N], N, 1);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    /*------------------------------*/


    /* Step 3.b - Collect Data */  
    MPI_Gather(&subC[0], N*div, complexStruct, &C[0], N*div, complexStruct, 0, MPI_COMM_WORLD);
    /*------------------------------*/


    /* Step 3.c - Transpose */  
    if (my_rank == 0) {
        transpose(C, N);
    }
    /*------------------------------*/


    /* Step 3.d - Distribute Data for Step 1 after transpose*/  
    MPI_Scatter(&C[0], N*div, complexStruct, &subC[0], N*div, complexStruct, 0, MPI_COMM_WORLD);
    /*------------------------------*/


    /* Step 3.e - Perform 1D FFT on C after transpose */  
    for( i=0; i<div; i++){   
        c_fft1d(&subC[i*N], N, 1);
    }
    /*------------------------------*/


    /* Collect Result Data */  
    MPI_Gather(&subC[0], N*div, complexStruct, &C[0], N*div, complexStruct, 0, MPI_COMM_WORLD);
    /*------------------------------*/


    if (my_rank == 0) {
        /* Stop Clock */
        gettimeofday(&etstop, &tzdummy);
        etstop2 = times(&cputstop);
        endtime = MPI_Wtime();
        printf("Stopped clock.\n");
        usecstart = (unsigned long long)etstart.tv_sec * 1000000 + etstart.tv_usec;
        usecstop = (unsigned long long)etstop.tv_sec * 1000000 + etstop.tv_usec;

        char output[15] = "out_test";
        writeToFile(output, C);

        /* Display timing results */
        printf("\nMPI_Wtime = %f s.\n", (endtime - starttime));

        /* Display timing results */
        printf("\nElapsed time = %g s.\n",
        (float)(usecstop - usecstart)/(float)1000000);

        printf("(CPU times are accurate to the nearest %g ms)\n",
        1.0/(float)CLOCKS_PER_SEC * 1000.0);
        printf("My total CPU time for parent = %g ms.\n",
        (float)( (cputstop.tms_utime + cputstop.tms_stime) -
            (cputstart.tms_utime + cputstart.tms_stime) ) /
        (float)CLOCKS_PER_SEC * 1000);
        printf("My system CPU time for parent = %g ms.\n",
        (float)(cputstop.tms_stime - cputstart.tms_stime) /
        (float)CLOCKS_PER_SEC * 1000);
        printf("My total CPU time for child processes = %g ms.\n",
        (float)( (cputstop.tms_cutime + cputstop.tms_cstime) -
            (cputstart.tms_cutime + cputstart.tms_cstime) ) /
        (float)CLOCKS_PER_SEC * 1000);
            /* Contrary to the man pages, this appears not to include the parent */
        printf("--------------------------------------------\n");
    }

    MPI_Finalize();

    free(A);
    free(B);
    free(C);

    return 0;
}