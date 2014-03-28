///////////          Assginment 3            //////////
/////////// Parallel Matrix Multiply on BG/Q  /////////
///////////////////////////////////////////////////////

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>


#if defined(__i386__)

static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int x;
  __asm__ volatile (".byte 0x0f, 0x31" : "=A" (x));
  return x;
}
#elif defined(__x86_64__)


static __inline__ unsigned long long rdtsc(void)
{
  unsigned hi, lo;
  __asm__ __volatile__ ("rdtsc" : "=a"(lo), "=d"(hi));
  return ( (unsigned long long)lo)|( ((unsigned long long)hi)<<32 );
}

#elif defined(__powerpc__)
static __inline__ unsigned long long rdtsc(void)
{
  unsigned long long int result=0;
  unsigned long int upper, lower,tmp;
  __asm__ volatile(
                "0:                  \n"
                "\tmftbu   %0           \n"
                "\tmftb    %1           \n"
                "\tmftbu   %2           \n"
                "\tcmpw    %2,%0        \n"
                "\tbne     0b         \n"
                : "=r"(upper),"=r"(lower),"=r"(tmp)
		   );
  result = upper;
  result = result<<32;
  result = result|lower;

  return(result);
}
#endif

/***********************************************************************/
/* START: MT 19937******************************************************/
/***********************************************************************/

/* Period parameters */  
#define N 624
#define M 397
#define MATRIX_A 0x9908b0dfUL   /* constant vector a */
#define UPPER_MASK 0x80000000UL /* most significant w-r bits */
#define LOWER_MASK 0x7fffffffUL /* least significant r bits */

static unsigned long mt[N]; /* the array for the state vector  */
static int mti=N+1; /* mti==N+1 means mt[N] is not initialized */

/* initializes mt[N] with a seed */
void init_genrand(unsigned long s)
{
    mt[0]= s & 0xffffffffUL;
    for (mti=1; mti<N; mti++) {
        mt[mti] = 
	    (1812433253UL * (mt[mti-1] ^ (mt[mti-1] >> 30)) + mti); 
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
        /* In the previous versions, MSBs of the seed affect   */
        /* only MSBs of the array mt[].                        */
        /* 2002/01/09 modified by Makoto Matsumoto             */
        mt[mti] &= 0xffffffffUL;
        /* for >32 bit machines */
    }
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
/* slight change for C++, 2004/2/26 */
void init_by_array(unsigned long init_key[], int key_length)
{
    int i, j, k;
    init_genrand(19650218UL);
    i=1; j=0;
    k = (N>key_length ? N : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1664525UL))
          + init_key[j] + j; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++; j++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=N-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 30)) * 1566083941UL))
          - i; /* non linear */
        mt[i] &= 0xffffffffUL; /* for WORDSIZE > 32 machines */
        i++;
        if (i>=N) { mt[0] = mt[N-1]; i=1; }
    }

    mt[0] = 0x80000000UL; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0,0xffffffff]-interval */
unsigned long genrand_int32(void)
{
    unsigned long y;
    static unsigned long mag01[2]={0x0UL, MATRIX_A};
    /* mag01[x] = x * MATRIX_A  for x=0,1 */

    if (mti >= N) { /* generate N words at one time */
        int kk;

        if (mti == N+1)   /* if init_genrand() has not been called, */
            init_genrand(5489UL); /* a default initial seed is used */

        for (kk=0;kk<N-M;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+M] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        for (;kk<N-1;kk++) {
            y = (mt[kk]&UPPER_MASK)|(mt[kk+1]&LOWER_MASK);
            mt[kk] = mt[kk+(M-N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
        }
        y = (mt[N-1]&UPPER_MASK)|(mt[0]&LOWER_MASK);
        mt[N-1] = mt[M-1] ^ (y >> 1) ^ mag01[y & 0x1UL];

        mti = 0;
    }
  
    y = mt[mti++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* generates a random number on [0,0x7fffffff]-interval */
long genrand_int31(void)
{
    return (long)(genrand_int32()>>1);
}

/* generates a random number on [0,1]-real-interval */
double genrand_real1(void)
{
    return genrand_int32()*(1.0/4294967295.0); 
    /* divided by 2^32-1 */ 
}

/* generates a random number on [0,1)-real-interval */
double genrand_real2(void)
{
    return genrand_int32()*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on (0,1)-real-interval */
double genrand_real3(void)
{
    return (((double)genrand_int32()) + 0.5)*(1.0/4294967296.0); 
    /* divided by 2^32 */
}

/* generates a random number on [0,1) with 53-bit resolution*/
double genrand_res53(void) 
{ 
    unsigned long a=genrand_int32()>>5, b=genrand_int32()>>6; 
    return(a*67108864.0+b)*(1.0/9007199254740992.0); 
} 
/* These real versions are due to Isaku Wada, 2002/01/09 added */

/***********************************************************************/
/* END: MT 19937 *******************************************************/
/***********************************************************************/

int main (int argc, char** argv)
{
  unsigned long long start_cycles=rdtsc(), 
    finish_cycles=rdtsc(), 
    start_cycles_Irecv=rdtsc(),
    finish_cycles_Irecv=rdtsc(),
    start_cycles_Isend=rdtsc(),
    finish_cycles_Isend=rdtsc(),
    start_cycles_mm=rdtsc(),
    finish_cycles_mm=rdtsc();

  double total_cycles=0,
    total_cycles_Irecv=0,
    total_cycles_Isend=0,
    total_cycles_mm=0;


  double clock_rate =2667000000.0  ;  //change to 2667000000.0/1600000000.0
  double execution_time=0,
    Isend_time=0,
    Isend_time_max=0,
    Isend_time_min=0,
    Isend_time_sum=0,
    Isend_time_avg=0,
    Irecv_time=0,
    Irecv_time_max=0,
    Irecv_time_min=0,
    Irecv_time_sum=0,
    Irecv_time_avg=0,
    mm_time=0,
    mm_time_max=0,
    mm_time_min=0,
    mm_time_sum=0,
    mm_time_avg=0,
    Isend_bandwidth=0,
    Isend_bandwidth_max=0,
    Isend_bandwidth_min=0,
    Isend_bandwidth_sum=0,
    Isend_bandwidth_avg=0,
    Irecv_bandwidth=0,
    Irecv_bandwidth_max=0,
    Irecv_bandwidth_min=0,
    Irecv_bandwidth_sum=0,
    Irecv_bandwidth_avg=0;

  unsigned int matrix_size = 4096;
  unsigned int partitioned_size;

  double **A=NULL;
  double **B=NULL;
  double *B_send=NULL;
  double *B_recv=NULL;
  double **C=NULL;

  int mysize,
    myrank,
    Dest,
    Source,
    count,
    flag_send,
    flag_recv;

  MPI_Request request_send, request_recv;
  MPI_Status status_send, status_recv;
    
  start_cycles = rdtsc();

  //MPI initialization
  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD, &mysize);
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

  //get the partitioned size N/P
  partitioned_size = matrix_size / mysize;

  unsigned long rng_init_seeds[6] = { 0x0, 0x123, 0x234, 0x345, 0x456, 0x789 };
  unsigned long rng_init_length = 6;

  rng_init_seeds[0] = myrank;
  init_by_array(rng_init_seeds, rng_init_length);
  
  int i,j,k;

  //each rank generates N/P * N part of A
  A = (double **)calloc( partitioned_size, sizeof(double*));
  for (i=0; i < partitioned_size; i++)
    A[i] = (double *)calloc( matrix_size, sizeof(double));

  //each rank generates N * N/P part of B
  B = (double **)calloc( matrix_size, sizeof(double*));
  for (i=0; i < matrix_size; i++)
    B[i] = (double *)calloc( partitioned_size, sizeof(double));

  //B_send, B_recv for sending/receiving
  B_send = (double *)calloc( matrix_size*partitioned_size, sizeof(double*));
  B_recv = (double *)calloc( matrix_size*partitioned_size, sizeof(double*));

  //each rank computes N/P * N part of C
  C = (double **)calloc( partitioned_size, sizeof(double*));
  for (i=0; i < partitioned_size; i++)
    C[i] = (double *)calloc( matrix_size, sizeof(double));

  //initialize A & B
  for (i=0; i<partitioned_size; i++) {
    for (j=0; j<matrix_size; j++) {
      A[i][j] = genrand_res53();
      B[j][i] = genrand_res53();
      C[i][j] = 0;
    }
  }

/*
if (myrank == 0) {
  printf("rank 0 Matrix A is \n");
  for (i=0;i<partitioned_size; i++) {
    for (j=0; j<matrix_size; j++) {
      printf("%f ", A[i][j]);
    }
    printf("\n");
  }

  printf("rank 0 Matrix B is \n");
  for (i=0;i<matrix_size; i++) {
    for (j=0; j<partitioned_size; j++) {
      printf("%f ", B[i][j]);
    }
    printf("\n");
  }
}
*/

  //do the matrix multiplication
  start_cycles_mm = rdtsc();
  for (i=0; i<partitioned_size; i++) {
    for (j=0; j<partitioned_size; j++) {
      for (k=0; k<matrix_size; k++) {
        C[i][j+partitioned_size*myrank] += A[i][k] * B[k][j];
      }
    }
  }
  finish_cycles_mm = rdtsc();
  total_cycles_mm += (double)finish_cycles_mm - (double)start_cycles_mm;

 
  /*perform MPI_Isend/Irecv pair between processors in a "ring"*/
  count = matrix_size * partitioned_size;
  Source = (myrank - 1 + mysize) % mysize;
  Dest = (myrank + 1) % mysize;
  int loop=1;
  do {
    for (i=0; i<matrix_size; i++) 
      for (j=0; j<partitioned_size; j++)
        B_send[i*partitioned_size+j] = B[i][j];

    start_cycles_Irecv = rdtsc();
//    printf("start Irecv cycles is %lld \n", start_cycles_Irecv);
    MPI_Irecv (B_recv, count, MPI_DOUBLE, Source, MPI_ANY_TAG, MPI_COMM_WORLD, &request_recv);

    start_cycles_Isend = rdtsc();
//    printf("start Isend cycles is %lld \n", start_cycles_Isend);
    MPI_Isend (B_send, count, MPI_DOUBLE, Dest, 1, MPI_COMM_WORLD, &request_send);

    flag_recv = false;
    flag_send = false;
    while (flag_recv == false || flag_send == false) { 
      if (flag_recv == false) {
        MPI_Test (&request_recv, &flag_recv, &status_recv);
        if (flag_recv == true) { 
          finish_cycles_Irecv = rdtsc();
          total_cycles_Irecv += (double)finish_cycles_Irecv - (double)start_cycles_Irecv; 
        }
//    printf("finish Irecv cycles is %lld \n", finish_cycles_Irecv);
      }
      if (flag_send == false) {
        MPI_Test (&request_send, &flag_send, &status_send);
        if (flag_send == true) {
          finish_cycles_Isend = rdtsc();
          total_cycles_Isend += (double)finish_cycles_Isend - (double)start_cycles_Isend;
        }
//    printf("finish Isend cycles is %lld \n", finish_cycles_Isend);
      }
    } 


//    printf("total Irecv cycles is %lld \n", total_cycles_Irecv); 
  
//    printf("total Isend cycles is %lld \n", total_cycles_Isend);
//printf("test\n");

    for (i=0; i<matrix_size; i++) {
      for (j=0; j<partitioned_size; j++) {
        B[i][j] = B_recv[i*partitioned_size+j];
      }
    }

    start_cycles_mm = rdtsc();
    for (i=0; i<partitioned_size; i++) {
      for (j=0; j<partitioned_size; j++) {
        for (k=0; k<matrix_size; k++) {
          C[i][(j+partitioned_size*(myrank-loop)+matrix_size)%matrix_size] += A[i][k] * B[k][j];
        }
      }
    } 
//printf("test1\n");   
    finish_cycles_mm = rdtsc();
    total_cycles_mm += (double)finish_cycles_mm - (double)start_cycles_mm;

    loop ++;
  } while (loop < mysize);
  /**/

/*
if (myrank == 0) {
  printf("rank 0 Matrix C is \n");
  for (i=0;i<partitioned_size; i++) {
    for (j=0; j<matrix_size; j++) {
      printf("%f ", C[i][j]);
    }
    printf("\n");
  }
}
*/


  Isend_time = total_cycles_Isend / clock_rate;
//printf("Isend time is %f \n", Isend_time);
  Irecv_time = total_cycles_Irecv / clock_rate;
//printf("Irecv time is %f \n", Irecv_time);
  mm_time = total_cycles_mm / clock_rate;

  //compute max, min, sum of time
  MPI_Allreduce (&Isend_time, &Isend_time_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (&Isend_time, &Isend_time_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce (&Isend_time, &Isend_time_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce (&Irecv_time, &Irecv_time_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (&Irecv_time, &Irecv_time_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce (&Irecv_time, &Irecv_time_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce (&mm_time, &mm_time_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (&mm_time, &mm_time_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce (&mm_time, &mm_time_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //compute average of time
  Isend_time_avg = Isend_time_sum / mysize;
  Irecv_time_avg = Irecv_time_sum / mysize;
  mm_time_avg = mm_time_sum / mysize;

  //compute max, min, sum of bandwidth for MPI sending/receving
  Isend_bandwidth = (mysize - 1) * partitioned_size * matrix_size * sizeof(double) / Isend_time;
  Irecv_bandwidth = (mysize - 1) * partitioned_size * matrix_size * sizeof(double) / Irecv_time;

  MPI_Allreduce (&Isend_bandwidth, &Isend_bandwidth_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (&Isend_bandwidth, &Isend_bandwidth_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce (&Isend_bandwidth, &Isend_bandwidth_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  MPI_Allreduce (&Irecv_bandwidth, &Irecv_bandwidth_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  MPI_Allreduce (&Irecv_bandwidth, &Irecv_bandwidth_min, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
  MPI_Allreduce (&Irecv_bandwidth, &Irecv_bandwidth_sum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

  //compute average of bandwidth
  Isend_bandwidth_avg = Isend_bandwidth_sum / mysize;
  Irecv_bandwidth_avg = Irecv_bandwidth_sum / mysize;
 
  MPI_Finalize();
  finish_cycles = rdtsc();
  total_cycles = finish_cycles - start_cycles;
  execution_time = total_cycles / clock_rate;

  //print out the result
  if (myrank == 0) {
    printf ("The total execution time is %f seconds \n", execution_time);
    printf ("\n");
    printf ("The MAX Isend time is %f seconds \n", Isend_time_max);
    printf ("The MIN Isend time is %f seconds \n", Isend_time_min);
    printf ("The AVERAGE Isend time is %f seconds \n", Isend_time_avg);
    printf ("\n");
    printf ("The MAX Irecv time is %f seconds \n", Irecv_time_max);
    printf ("The MIN Irecv time is %f seconds \n", Irecv_time_min);
    printf ("The AVERAGE Irecv time is %f seconds \n", Irecv_time_avg);
    printf ("\n");
    printf ("The MAX matrix multiplication time is %f seconds \n", mm_time_max);
    printf ("The MIN matrix multiplication time is %f seconds \n", mm_time_min);
    printf ("The AVERAGE matrix multiplication time is %f seconds \n", mm_time_avg);
    printf ("\n");
    printf ("The MAX Isend bandwidth is %f bytes \n", Isend_bandwidth_max);
    printf ("The MIN Isend bandwidth is %f bytes \n", Isend_bandwidth_min);
    printf ("The AVERAGE Isend bandwidth is %f bytes \n", Isend_bandwidth_avg);
    printf ("\n");
    printf ("The MAX Irecv bandwidth is %f bytes \n", Irecv_bandwidth_max);    
    printf ("The MIN Irecv bandwidth is %f bytes \n", Irecv_bandwidth_min);
    printf ("The AVERAGE Irecv bandwidth is %f bytes \n", Irecv_bandwidth_avg);
  }
  return 0;
}
