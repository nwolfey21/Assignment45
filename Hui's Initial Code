#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

/***********************************************************************/
/* RTDSC ***************************************************************/
/***********************************************************************/

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

#define MATRIX_SIZE 4096 // change to 16384 for Blue Gene/Q
#define CLOCK_RATE 2666700000 // change to 1600000000.0 for Blue Gene/Q

int num_of_processors, myrank;


// print out matrix
void printMatrix(double** matrix, int rows, int columns) {
	printf("Rank %d:\n", myrank);	

	int i, j;
	for (i = 0; i < rows; i++) {
		for (j = 0; j < columns; j++) {
			printf("%lf ", matrix[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}


// create the matrix with calloc
double** createMatrix(int rows, int columns) {
	double* data = (double*) calloc(rows * columns, sizeof(double));
	double** matrix = (double**) calloc(rows, sizeof(double*));

	int i;
	for (i = 0; i < rows; i++) {
		matrix[i] = &(data[columns*i]);
	}

	return matrix;
}

// assign matrix with random generator
void populateMatrix(double** matrix, int rows, int columns) {
	int i, j;
	for(i = 0; i < rows; i++) {
		for(j = 0; j < columns; j++) {
			matrix[i][j] = genrand_res53();
		}
	}
}


void printActivityTimes(long long mpi_isend_cycles, long long mpi_irecv_cycles, long long matrix_multiply_cycles) {
	double mpi_isend_time = mpi_isend_cycles * 1.0 / CLOCK_RATE;
	double mpi_irecv_time = mpi_irecv_cycles * 1.0 / CLOCK_RATE;
	double matrix_multiply_time = matrix_multiply_cycles * 1.0 / CLOCK_RATE;
	double min_send_time;
	double max_send_time;
	double total_send_time;
	double send_average;
	double min_recv_time;
	double max_recv_time;
	double total_recv_time;
	double recv_average;
	double min_mm_time;
	double max_mm_time;
	double total_mm_time;
	double mm_average;

	// MPI_Allreduce
	MPI_Allreduce(&mpi_isend_time, &min_send_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&mpi_isend_time, &max_send_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&mpi_isend_time, &total_send_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	send_average = total_send_time / num_of_processors;

	MPI_Allreduce(&mpi_irecv_time, &min_recv_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&mpi_irecv_time, &max_recv_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&mpi_irecv_time, &total_recv_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	recv_average = total_recv_time / num_of_processors;

	MPI_Allreduce(&matrix_multiply_time, &min_mm_time, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);
	MPI_Allreduce(&matrix_multiply_time, &max_mm_time, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
	MPI_Allreduce(&matrix_multiply_time, &total_mm_time, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	mm_average = total_mm_time / num_of_processors;

	// calculate bandwidth
	double min_send_bandwidth = MATRIX_SIZE * MATRIX_SIZE / min_send_time;
	double max_send_bandwidth = MATRIX_SIZE * MATRIX_SIZE / max_send_time;
	double send_bandwidth_average = MATRIX_SIZE * MATRIX_SIZE / send_average;

	double min_recv_bandwidth = MATRIX_SIZE * MATRIX_SIZE / min_recv_time;
	double max_recv_bandwidth = MATRIX_SIZE * MATRIX_SIZE / max_recv_time;
	double recv_bandwidth_average = MATRIX_SIZE * MATRIX_SIZE / recv_average;

    // print out the MPI and multiplication statistics
	if (myrank == 0) {
	   printf ("===================The send time ================================\n");
	   printf ("The MAX send time is %f seconds\n", max_send_time);
	   printf ("The MIN send time is %f seconds\n", min_send_time);
	   printf ("The average send time is %f seconds\n", send_average);

	   printf ("===================The recv time ================================\n");
	   printf ("The MAX recv time is %f seconds\n", max_recv_time);
	   printf ("The MIN recv time is %f seconds\n", min_recv_time);
	   printf ("The average recv time is %f seconds\n", recv_average);

	   printf("===================The Matrix Multiply time ======================\n");
       printf ("The MAX matrix multiply time is %f seconds\n", max_mm_time);
	   printf ("The MIN matrix multiply time is %f seconds\n", min_mm_time);
	   printf ("The average matrix multiply time is %f seconds\n", mm_average);

	   printf ("=================Send bandwidth ============\n");
	   printf ("The MAX send bandwidth is %f Mbps\n", max_send_bandwidth);
	   printf ("The MIN send bandwidth is %f Mbps\n", min_send_bandwidth);
	   printf ("The average send bandwidth is %f Mbps\n", send_bandwidth_average);
	   
	   printf ("=================Recv bandwidth ============\n");
	   printf ("The MAX recv bandwidth is %f Mbps\n", max_recv_bandwidth);
	   printf ("The MIN recv bandwidth is %f Mbps\n", min_recv_bandwidth);
	   printf ("The average recv bandwidth is %f Mbps\n", recv_bandwidth_average);
	}
}

void calculate(double** A, double** B, double** C, int start_row, int start_col, int size)
{
	int i, j, k, m;
	unsigned long long last = rdtsc();
	unsigned long long current = rdtsc();

	for(i = 0; i < size; i++) {
		for(j = 0; j < size; j++) {
			for(k = 0; k < MATRIX_SIZE; k++) {
				C[i][j+start_col] += A[i][k] * B[k][j];
			}
		}
	}
}

void matrixMultiply(double** A, double** B, double** C, int start_row, int start_col, int size) {
	// determine the data source and destination
	int dest = myrank + 1;
	if (dest == num_of_processors) {
		dest = 0;
	}
	int source = myrank - 1;
	if (source == -1) {
		source = num_of_processors - 1;
	}

	long long mpi_isend_cycles = 0;
	long long mpi_irecv_cycles = 0;
	long long matrix_multiply_cycles = 0;
	
	MPI_Request send_request, recv_request;
	int i;

	for (i = 0; i < num_of_processors; i++) {

	// start the cycle counter
		long long matrix_multiply_start_cycles = rdtsc();
		calculate(A, B, C, start_row, start_col, size);
		matrix_multiply_cycles += rdtsc() - matrix_multiply_start_cycles;

		double** tempB = createMatrix(MATRIX_SIZE, MATRIX_SIZE / num_of_processors);
		int send_flag = 0;
		int recv_flag = 0;
		long long mpi_isend_start_cycles = rdtsc();
		long long mpi_irecv_start_cycles = rdtsc();
		MPI_Isend(&(B[0][0]), MATRIX_SIZE * MATRIX_SIZE / num_of_processors, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD, &send_request);
		MPI_Irecv(&(tempB[0][0]), MATRIX_SIZE * MATRIX_SIZE / num_of_processors, MPI_DOUBLE, source, MPI_ANY_TAG, MPI_COMM_WORLD, &recv_request);

		while (1) {
			if (send_flag != 0 && recv_flag != 0) {
				B = tempB;
				break;
			}
			if (send_flag == 0) {	
				MPI_Status send_status;
				MPI_Test(&send_request, &send_flag, &send_status);
				if (send_flag != 0) {
					mpi_isend_cycles = rdtsc() - mpi_isend_start_cycles;
				}
			}
			if (recv_flag == 0) {
				MPI_Status recv_status;
				MPI_Test(&recv_request, &recv_flag, &recv_status);
				if (recv_flag != 0) {
					mpi_irecv_cycles = rdtsc() - mpi_irecv_start_cycles;
				}
			}
		}

		// adjust the start row and start column
		start_row += size;
		start_col += size;
		if (start_row == MATRIX_SIZE) {
			start_row = 0;
		}
		if (start_col == MATRIX_SIZE) {
			start_col = 0;
		}
	}

	printActivityTimes(mpi_isend_cycles, mpi_irecv_cycles, matrix_multiply_cycles);
}


int main(int argc, char *argv[]) {
	long long start_cycles = rdtsc();

	// MPI Initialize
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &num_of_processors);
	MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
	
	unsigned long rng_init_seeds[6] = {myrank, 0x123, 0x234, 0x345, 0x456, 0x789};
	unsigned long rng_init_length = 6;
	init_by_array(rng_init_seeds, rng_init_length);

	double** A = createMatrix(MATRIX_SIZE / num_of_processors, MATRIX_SIZE);
	double** B = createMatrix(MATRIX_SIZE, MATRIX_SIZE / num_of_processors);
	double** C = createMatrix(MATRIX_SIZE / num_of_processors, MATRIX_SIZE);
	populateMatrix(A, MATRIX_SIZE / num_of_processors, MATRIX_SIZE);
	populateMatrix(B, MATRIX_SIZE, MATRIX_SIZE / num_of_processors);

	// obtain the matrix
	matrixMultiply(A, B, C, MATRIX_SIZE / num_of_processors * myrank, MATRIX_SIZE / num_of_processors * myrank, MATRIX_SIZE / num_of_processors);
	if (myrank == 0) {
	}
	MPI_Finalize();

	long long total_cycles = rdtsc() - start_cycles;
	double execution_time = total_cycles * 1.0 / CLOCK_RATE;

	// get the total execution time
	if (myrank == 0) {	
       printf ("====================The execution time============================\n");
	   printf ("The total execution time is %f seconds\n", execution_time);
	   printf ("The no. of processors is %d\n", num_of_processors);

	   	}

	return 0;
}
