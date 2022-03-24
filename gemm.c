/***************************************************************************
Copyright (c) 2014, The OpenBLAS Project
All rights reserved.
Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:
1. Redistributions of source code must retain the above copyright
notice, this list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions and the following disclaimer in
the documentation and/or other materials provided with the
distribution.
3. Neither the name of the OpenBLAS project nor the names of
its contributors may be used to endorse or promote products
derived from this software without specific prior written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE OPENBLAS PROJECT OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*****************************************************************************/

//#include "bench.h"
#include <stdio.h>
#include <time.h>
#include <malloc.h>
#include <sched.h>
#include <math.h>
#include "cblas.h"

#define COMPSIZE  1
#define RAND_MAX 1000
struct timespec start = { 0, 0 }, stop = { 0, 0 };

double getsec()
{
    return (double)(stop.tv_sec - start.tv_sec) + (double)((stop.tv_nsec - start.tv_nsec)) * 1.e-9;

}

void begin() {
    clock_gettime(CLOCK_REALTIME, &start);
}

void end() {
    clock_gettime(CLOCK_REALTIME, &stop);
}

//int main(int argc, char *argv[]){
int main(){
  printf("###############start\n");
  printf("###############start\n");
  printf("###############start\n");

  double *a, *b;
  double *c;
  double alpha = 1.0;
  double beta  = 0.0;
  char transa = 'N';
  char transb = 'N';
  int m, n, k, i, j, lda, ldb, ldc;
  int loops = 1;
  char *p;


  double time1, timeg;
/*
  m = atoi(argv[1]);
  n = atoi(argv[2]);
  k = atoi(argv[3]);
  printf("%s %s %s\n", argv[1], argv[2], argv[3]);*/
  m=300;
  n=300;
  k=300;
  printf("%d %d %d\n", m, n, k);
  //fprintf(stderr, "From : %3d  To : %3d Step=%d : Transa=%c : Transb=%c\n", from, to, step, transa, transb);

  
  if (( a = (double *)malloc(sizeof(double) * m * k * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  if (( b = (double *)malloc(sizeof(double) * k * n * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }
  if (( c = (double *)malloc(sizeof(double) * m * n * COMPSIZE)) == NULL) {
    fprintf(stderr,"Out of Memory!!\n");exit(1);
  }



  for (i = 0; i < m * k * COMPSIZE; i++) {
    a[i] = ((double) rand() / (double) RAND_MAX) - 0.5;
  }
  for (i = 0; i < k * n * COMPSIZE; i++) {
    b[i] = ((double) rand() / (double) RAND_MAX) - 0.5;
  }
  for (i = 0; i < m * n * COMPSIZE; i++) {
    c[i] = ((double) rand() / (double) RAND_MAX) - 0.5;
  }

  fprintf(stderr, "          SIZE                   Flops             Time\n");


    
  timeg=0;


  if (transa == 'N') { lda = m; }
  else { lda = k; }
  if (transb == 'N') { ldb = k; }
  else { ldb = n; }
  ldc = m;

  printf(" M=%4d, N=%4d, K=%4d : ", (int)m, (int)n, (int)k);
  begin();

  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m, n, k, alpha, a, lda, b, ldb, beta, c, ldc);

  end();
  time1 = getsec();

  timeg = time1/loops;
  printf(
	   " %10.2f GFlops/s %10.6f sec\n",
	   COMPSIZE * COMPSIZE * 2. * (double)k * (double)m * (double)n / timeg * 1.e-9, time1);
  for(i=0; i<5; i++){
      double sum = 0;
      for(j=0; j<k; j++){
          sum += a[j] * b[j*n+i];
      }
      printf("@@@ test %f %f \n", c[i], sum);
  }
  printf("###############end\n");
  printf("###############end\n");
  printf("###############end\n");
  return 0;
}

