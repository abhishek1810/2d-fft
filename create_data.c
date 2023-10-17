#include <assert.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>

typedef struct {float r; float i;} complex;
static complex ctmp;

void writeToFile(char fileName[15], int N){
    FILE *fp = fopen(fileName, "w");
    int i, j;
    for (i=0;i<N;i++) {
        for (j=0;j<N;j++){
            fprintf(fp,"   %.7e", 1.0);
        }
        fprintf(fp,"\n");
    }
    fclose(fp);
}

void readFile(char fileName[15], complex *array, int N){
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

int main(int argc, char **argv){
    writeToFile("sample/4_im1", 1024);
    writeToFile("sample/4_im2", 1024);

    writeToFile("sample/5_im1", 2048);
    writeToFile("sample/5_im2", 2048);

    writeToFile("sample/6_im1", 4096);
    writeToFile("sample/6_im2", 4096);

    complex *array1, *array2;

    array1 = (complex *)malloc(sizeof(complex)*1024*1024);
    array2 = (complex *)malloc(sizeof(complex)*2048*2048);

    readFile("sample/4_im1", array1, 1024);
    readFile("sample/5_im1", array2, 2048);

    printf("4_im1 A[1024][1024] - %f\n", array1[1024*1024-1].r);
    printf("5_im1 A[2048][2048] - %f\n", array2[2048*2048-1].r);

    return 0;
}