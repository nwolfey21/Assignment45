#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <mpi.h>
#include <math.h>
#include <pthread.h>

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

typedef struct 
{
  double** A;
  double** B;
  double** C;
  int numthrds;
  int datapernode;
} MULTIPDATA;

/***********************Start Global Variables**************************/
#define MAXTHRDS 4
MULTIPDATA mdata;
pthread_t threads[MAXTHRDS];

unsigned long rng_init_seeds[6]={0x0, 0x123, 0x234, 0x345, 0x456, 0x789};
unsigned long rng_init_length=6;
double clock_rate=2000000000.0;//3200000000.0;//2666700000.0; // change to 1600000000.0 for Blue Gene/Q

int startcolumn = 0;
/************************End Global Variables***************************/

void *multip(void *arg)
{ 
    /*define and use local variables for convenience*/
    long mythrd;  
    int numthrds, dataperthread, numpernode, start, end;
    double** A, ** B, ** C;

    /*the number of threads and nodes defines the size of the matrix slice*/
    mythrd=(long)arg;

    numthrds = mdata.numthrds;
    datapernode = mdata.datapernode;
    dataperthread = datapernode / numthrds;
    A = mdata.A;
    B = mdata.B;
    C = mdata.C;
    start = dataperthread*mythrd;
    end = start + dataperthread;

    /*perform matrix multiply*/
    for (int i=start; i<end; i++) {
        for (int j=0; j<datapernode; j++) {
            for (int k=0; k<matrix_size; k++) {
                C[i][j+startcolumn] += A[i][k] * B[k][j];
            }
        }
    }
    
    pthread_exit((void*)0);
}


/***********************************************************************/
/* Start															   */
/* Standard matrix multiplication 									   */
/***********************************************************************/
/*
void matrix_multiply( int rank, int id, double **A, double **B, double **C, int size, int slice, int count )
{
    int i=0, j=0, k=0, jj=0;
    unsigned long long last=rdtsc();
    unsigned long long current=rdtsc();

	for( i = 0; i < slice; i++ )
	{
		if( i != 0 && i % 16 == 0 && rank == 0)
		{
			last = current;
			current = rdtsc();
			printf("Completed %d in time %lf secs \n", i, ((double)current - (double)last)/clock_rate);
		}
		for( j = 0; j < slice; j++ )
		{
			jj = j+count*slice;
			for( k = 0; k < size; k++ )
			{
				C[i][jj] += A[i][k] * B[k][j];
			}
//			if(rank == 0)
//			printf("rank:%d j:%d jj:%d C[%d][%d]=%f\n",rank,j,jj,i,jj,C[i][jj]);
		}
	}
}
*/
/***********************************************************************/
/* End  															   */
/* Standard matrix multiplication 									   */
/***********************************************************************/

/*void parallelIO(MPI_File fh, double **C, int matrix_size, int datapernode, int myRank, char *filename)
{
	MPI_Offset offset = 0;
	MPI_Status status;
	double *tempC = (double *)calloc(datapernode,sizeof(double*));
	int count=0,i,j,ii,err=0;
	for(i=0;i<datapernode;i++)
	{	
		MPI_Barrier(MPI_COMM_WORLD);
		for(j=0;j<matrix_size;j++)
		{
			tempC[j] = C[i][j];
		}
		ii = i+myRank*datapernode;		//compute i index into global C matrix
		offset = ii*matrix_size;	//Set the offset into the file for first data location
		MPI_File_write_at(fh, offset*8, tempC, matrix_size, MPI_DOUBLE, &status);	//Write vector C[i]
*/		/****Start Verify Data was written to file****/
/*		err = MPI_Get_elements(&status, MPI_BYTE, &count);
		if(err != MPI_SUCCESS)
			printf("Write Failed\n");
		if(count != matrix_size*sizeof(double))
			printf("Did not write the same number of bytes as requested\n");
		else
			printf("Rank:%d Wrote %d bytes(%d doubles) in %s at offset %d\n", myRank,count,count/8,filename,(int)offset);
*/		/****End Verify Data was written to file****/
//	}
//}

/***********************************************************************/
/* Start															   */
/* Main 									   						   */
/***********************************************************************/
int main(int argc, char** argv)
{
	/***Matrix Variables***/
	unsigned int matrix_size=4;
	double **A=NULL;
	double **B=NULL;
//	double **outB=NULL;
	double *outB=NULL;
	double *inB=NULL;
	double **C=NULL;
	/***Other Variables***/
	int i, ii, j, n, datapernode, flag, tasks, nodes, ranksPerFile=0;
	int commSize, myRank, send, recv, count, id, blocking=0;
	long long start_time, time_cycles, sum_cycles, transfer_cycles, start_transfer;
	long long start_compute, compute_cycles, min_cycles, max_cycles;
	long long sum_copy, sum_compute, sum_transfer;
	long long max_copy, max_compute, max_transfer;
	long long min_copy, min_compute, min_transfer;
	long long start_copy, copy_cycles;
    unsigned long long last=rdtsc();
    unsigned long long current=rdtsc();
	compute_cycles = 0;
	copy_cycles = 0;
	MPI_Request srequest, rrequest;
	MPI_Status status;
	char filePath[75];
	FILE *dataFile;
	#define BLOCK_SIZE 8388608

	//---Setup output file---//
	for(i = 1; i < argc; ++i) 
	{
		if(argv[i][0] == '-') 
		{
			switch(argv[i][1]) 
			{
				case 't': 	tasks = atoi(argv[i+1]);
							break;
				case 'n': 	nodes = atoi(argv[i+1]);
							break;
				case 'r':	ranksPerFile = atoi(argv[i+1]);	//1=all ranks use one file, 4=4ranks per file, etc
							break;
				case 'b':	blocking = 1;					//1=write files to 8MB blocks, 0=write data compact
			}
		}
	}

	//***Initialize MPI***/
        int provided = 0;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	MPI_Comm_size(MPI_COMM_WORLD, &commSize);
	MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
	//---Compute Partial Matrix Initialization Parameters---//
	datapernode = matrix_size/commSize;
	id = myRank;
	outB = (double*) calloc(datapernode*matrix_size,sizeof(double));
	inB = (double*) calloc(datapernode*matrix_size,sizeof(double));
	rng_init_seeds[0] = myRank;
	init_by_array(rng_init_seeds, rng_init_length);
	sprintf(filePath, "data/data.dat");

        //    Pthreads    //
        int numthrds = MAXTHRDS;

        pthread_attr_t attr;
        pthread_attr_init(&attr);
        pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

	//-----Initialize A Matrix------//
	A = (double **)calloc( datapernode, sizeof(double*));
	for( i = 0; i < datapernode; i++ ) 
	{
		A[i] = (double *)calloc( matrix_size, sizeof(double));
	}
	//-----Initialize B Matrix------//
	B = (double **)calloc( matrix_size, sizeof(double*));
//	outB = (double **)calloc( matrix_size, sizeof(double*));
	for( i = 0; i < matrix_size; i++ ) 
	{
		B[i] = (double *)calloc( datapernode, sizeof(double));
//		outB[i] = (double *)calloc( datapernode, sizeof(double));
	}
	//-----Initialize C Matrix------//
	C = (double **)calloc( datapernode, sizeof(double*));
	for( i = 0; i < datapernode; i++ ) 
	{
		C[i] = (double *)calloc( matrix_size, sizeof(double));
	}

	if(myRank == 0)
	{
		printf("\n\nMatrix Size:%d datapernode:%d commSize:%d\n",matrix_size, datapernode, commSize);
	}

	//---Initialize Random Values---//
	for( i = 0; i < datapernode; i++ )  
	{
		for( j = 0; j < matrix_size; j++) 
		{
			A[i][j] = genrand_res53(); // (double)i*j;
			B[j][i] = genrand_res53(); //A[i][j] * A[i][j];
//			printf("rank:%d A[%d][%d]=%f\n",myRank,i,j,A[i][j]);
//			printf("rank:%d B[%d][%d]=%f\n",myRank,j,i,B[j][i]);
		}
	}

        mdata.A = A;
        mdata.B = B;
        mdata.C = C;
        mdata.numthrds = numthrds;
        mdata.datapernode = datapernode;

	if(myRank == 0)
		printf("\nComputing C Matrix\n");
	//--Start Timer--//
	start_time = rdtsc();
	//---Compute and Transfer Loop---//
	for(count=0; count<commSize; count++)
	{
		start_compute = rdtsc();
		//---Compute Partial Matrix Multiplication---//
		//matrix_multiply(myRank, id, A, B, C, matrix_size, datapernode, count);
           
                startcolumn = (datapernode * (myRank - count) + matrix_size) % matrix_size;
                for (i=0; i<numthrds; i++) {
			pthread_create(&threads[i], &attr, multip, (void*)i);
                }         
		for (i=0; i<numthrds; i++) {
			pthread_join(threads[i], NULL)
		}
 
		compute_cycles += rdtsc() - start_compute;
		if(count<commSize-1)
		{
			//---Calculate sender and reciever---//
			send = myRank+1;
			recv = myRank-1;
			if(myRank == commSize-1)
			{
				send = 0;
			}
			if(myRank == 0)
			{
				recv = commSize-1;
			}
			id = recv;
			//---Copy B matrix data to output buffer---//
			start_copy = rdtsc();
			for( i = 0; i < matrix_size; i++ )
			{
				for( j = 0; j < datapernode; j++ )
				{
					n = j + i*datapernode;
					outB[n] = B[i][j];
				}
			}
			copy_cycles += rdtsc() - start_copy;
			//---Send Output Data Buffer---//
			start_transfer = rdtsc();
			MPI_Isend(outB, datapernode*matrix_size, MPI_DOUBLE, send, 1, MPI_COMM_WORLD, &srequest);
			//---Recieve Data in B Matrix---//
			MPI_Irecv(inB, datapernode*matrix_size, MPI_DOUBLE, recv, MPI_ANY_TAG, MPI_COMM_WORLD, &rrequest);
			flag = 0;
			//---Verify Data was Received---//
			while(flag != 1)
			{
		  		MPI_Test(&rrequest,&flag,&status);
			}
			transfer_cycles += rdtsc() - start_transfer;
			start_copy = rdtsc();
			//---Copy input buffer data into B matrix---//
			for( i = 0; i < matrix_size; i++ )
			{
				for( j = 0; j < datapernode; j++ )
				{
					n = j + i*datapernode;
					B[i][j] = inB[n];
				}
			}
			copy_cycles += rdtsc() - start_copy;
		}
		if(myRank == 0 /*&& count%64 == 0*/)
		{
			last = current;
			current = rdtsc();
			printf("Completed %d out of %d in time %lf secs \n", count+1, commSize, ((double)current - (double)last)/clock_rate);
			for(i=0;i<datapernode;i++)
			{
				for(j=0;j<matrix_size;j++)
				{
//					printf("rank:%d C[%d][%d]=%f\n",myRank,i,j,C[i][j]);
				}
			}
		}
	}
	MPI_Barrier(MPI_COMM_WORLD);
	if(myRank ==0)
		printf("C Matrix Computed\n\n");

	/********Start Parallel I/O (Writing computed C matrix to file)********/
	if(myRank == 0)
		printf("Writing Data to File\n");
	MPI_File fh;
	MPI_Offset offset = 0;
	char filename[30];
	double *tempC = (double *)calloc(matrix_size,sizeof(double*));
	double *temp = (double *)calloc(matrix_size,sizeof(double*));
	int err;

	if(ranksPerFile !=0)	//All cases other than 1 file for all ranks
	{
		int fileID = floor(myRank/ranksPerFile); //ranksPerFIle many sequential MPI ranks will have the same fileID
		sprintf(filename, "data/cMatrix%d.dat",fileID);	//Brand each file with it's id
		err = MPI_File_open( MPI_COMM_SELF, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
		if (err) 								//Verify File was succesfully opened
		{
			printf( "Unable to open file %s\n",filename );
		}

		for(i=0;i<datapernode;i++)
		{	
			for(j=0;j<matrix_size;j++)
			{
				tempC[j] = C[i][j];
			}
			ii = i+myRank*datapernode;		//compute i index into global C matrix
			if(blocking)
				offset = ii*matrix_size*BLOCK_SIZE;		//Move offset over 8MB (1block)
			else
  				offset = ii*matrix_size*8;	//Set the offset into the file for first data location without blocking
			MPI_File_write_at(fh, offset, tempC, matrix_size, MPI_DOUBLE, &status);	//Write vector C[i]
			/****Start Verify Data was written to file****/
			err = MPI_Get_elements(&status, MPI_BYTE, &count);
			if(err != MPI_SUCCESS)
				printf("Write Failed\n");
			if(count != matrix_size*sizeof(double))
				printf("Did not write the same number of bytes as requested\n");
			else
				printf("Rank:%d Wrote %d bytes or %d double values at offset %d\n", myRank,count,count/8,(int)offset);
			/****End Verify Data was written to file****/
		}
	}
	else 	//Case with all ranks writing to one file
	{
		sprintf(filename, "data/cMatrix.dat");
		err = MPI_File_open( MPI_COMM_WORLD, filename, MPI_MODE_RDWR | MPI_MODE_CREATE, MPI_INFO_NULL, &fh );
		if (err) 								//Verify File was succesfully opened
		{
			printf( "Unable to open file \n" );
		}
//		parallelIO(fh, C, matrix_size, datapernode, myRank, filename);
		for(i=0;i<datapernode;i++)
		{	
			MPI_Barrier(MPI_COMM_WORLD);
			for(j=0;j<matrix_size;j++)
			{
				tempC[j] = C[i][j];
			}
			ii = i+myRank*datapernode;		//compute i index into global C matrix
			if(blocking)
				offset = ii*matrix_size*BLOCK_SIZE;		//Move offset over 8MB (1block)
			else
  				offset = ii*matrix_size*8;	//Set the offset into the file for first data location without blocking
			MPI_File_write_at(fh, offset, tempC, matrix_size, MPI_DOUBLE, &status);	//Write vector C[i]
			/****Start Verify Data was written to file****/
			err = MPI_Get_elements(&status, MPI_BYTE, &count);
			if(err != MPI_SUCCESS)
				printf("Write Failed\n");
			if(count != matrix_size*sizeof(double))
				printf("Did not write the same number of bytes as requested\n");
			else
				printf("Rank:%d Wrote %d bytes(%d doubles) in %s at offset %d\n", myRank,count,count/8,filename,(int)offset);
			/****End Verify Data was written to file****/
		}
	}
	/*********End Parallel I/O (Writing computed C matrix to file)*********/

	temp = (double *)calloc(matrix_size*matrix_size,sizeof(double*));
	MPI_Barrier(MPI_COMM_WORLD);
	if(myRank == 0)
	{
			printf("\nVerification\n");
			offset = 0;
			/***Start Verify correct data was written to file***/
			MPI_File_read_at(fh,offset,temp,matrix_size*matrix_size,MPI_DOUBLE,&status);
			/****Start Verify Data was written to file****/
			err = MPI_Get_elements(&status, MPI_BYTE, &count);
			if(err != MPI_SUCCESS)
				printf("----Read Failed\n");
			if(count != matrix_size*matrix_size*sizeof(double))
				printf("----Did not Read the same number of bytes as requested (read %d bytes)\n",count);
			else
				printf("----Rank:%d Read %d bytes or %d double values at offset %d\n", myRank,count,count/8,(int)offset);
			/****End Verify Data was written to file****/
			for(j=0;j<matrix_size*matrix_size;j++)
			{
				printf("----offset:%d temp[%d]:%f \n",(int)offset+8*j,j,temp[j]);
			}
			/***End Verify correct data was written to file***/
	}
//	printf("myrank:%d\n",myRank);
	if(myRank==0)
	MPI_File_close( &fh );					//Close the MPI File
//	printf("myrank:%d\n",myRank);

	time_cycles = rdtsc() - start_time;
//	printf("myrank0:%d\n",myRank);
	//---Total program execution time reduction---//
	MPI_Allreduce(&time_cycles,&sum_cycles,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
//	printf("myrank0.5:%d\n",myRank);
	MPI_Allreduce(&time_cycles,&max_cycles,1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&time_cycles,&min_cycles,1,MPI_LONG_LONG,MPI_MIN,MPI_COMM_WORLD);
//	printf("myrank1:%d\n",myRank);
	//---MPI transfer time reduction---//
	MPI_Allreduce(&transfer_cycles,&sum_transfer,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&transfer_cycles,&max_transfer,1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&transfer_cycles,&min_transfer,1,MPI_LONG_LONG,MPI_MIN,MPI_COMM_WORLD);
	//---Matrix Copy time Reduction---//
	MPI_Allreduce(&copy_cycles,&sum_copy,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&copy_cycles,&max_copy,1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&copy_cycles,&min_copy,1,MPI_LONG_LONG,MPI_MIN,MPI_COMM_WORLD);
	//---Matrix Compute Time Reduction---//
	MPI_Allreduce(&compute_cycles,&sum_compute,1,MPI_LONG_LONG,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(&compute_cycles,&max_compute,1,MPI_LONG_LONG,MPI_MAX,MPI_COMM_WORLD);
	MPI_Allreduce(&compute_cycles,&min_compute,1,MPI_LONG_LONG,MPI_MIN,MPI_COMM_WORLD);

	if(myRank == 0)
	{
		printf("\nWriting Timing File\n");
		//---Open Output File for writing---//
		if((dataFile = fopen(filePath,"a")) == NULL) 
		{
		    fprintf(stdout,"Couldn't open data output file\n");
		    exit(1);
		}
/*		printf("Cycles:%lld in time %lf secs\n",time_cycles,(double)time_cycles/clock_rate);
		printf("Max_Cycles:%lld in time %lf secs\n",max_cycles,(double)max_cycles/clock_rate);
		printf("Min_Cycles:%lld in time %lf secs\n",min_cycles,(double)min_cycles/clock_rate);
		printf("Avg_Cycles:%lld in time %lf secs\n",sum_cycles/commSize,(double)(sum_cycles/commSize)/clock_rate);
		printf("Total_Cycles:%lld in time %lf secs\n",sum_cycles,(double)sum_cycles/clock_rate);
		printf("Max_Transfer_Cycles:%lld in time %lf secs\n",max_transfer,(double)max_transfer/clock_rate);
		printf("Min_Transfer_Cycles:%lld in time %lf secs\n",min_transfer,(double)min_transfer/clock_rate);
		printf("Avg_Transfer_Cycles:%lld in time %lf secs\n",sum_transfer/commSize,(double)(sum_transfer/commSize)/clock_rate);
		printf("Total_Transfer_Cycles:%lld in time %lf secs\n",sum_transfer,(double)sum_transfer/clock_rate);
		printf("Max_Copy_Cycles:%lld in time %lf secs\n",max_copy,(double)max_copy/clock_rate);
		printf("Min_Copy_Cycles:%lld in time %lf secs\n",min_copy,(double)min_copy/clock_rate);
		printf("Avg_Copy_Cycles:%lld in time %lf secs\n",sum_copy/commSize,(double)(sum_copy/commSize)/clock_rate);
		printf("Total_Copy_Cycles:%lld in time %lf secs\n",sum_copy,(double)sum_copy/clock_rate);
		printf("Max_Compute_Cycles:%lld in time %lf secs\n",max_compute,(double)max_compute/clock_rate);
		printf("Min_Compute_Cycles:%lld in time %lf secs\n",min_compute,(double)min_compute/clock_rate);
		printf("Avg_Compute_Cycles:%lld in time %lf secs\n",sum_compute/commSize,(double)(sum_compute/commSize)/clock_rate);
		printf("Total_Compute_Cycles:%lld in time %lf secs\n",sum_compute,(double)sum_compute/clock_rate);
*/
		//---Write timing data to the output file---//
		fprintf(dataFile,"%d %d ",nodes, tasks);
		fprintf(dataFile,"%lld %lf ",time_cycles,(double)time_cycles/clock_rate);
		fprintf(dataFile,"%lld %lf ",max_cycles,(double)max_cycles/clock_rate);
		fprintf(dataFile,"%lld %lf ",min_cycles,(double)min_cycles/clock_rate);
		fprintf(dataFile,"%lld %lf ",sum_cycles/commSize,(double)(sum_cycles/commSize)/clock_rate);
		fprintf(dataFile,"%lld %lf ",sum_cycles,(double)sum_cycles/clock_rate);
		//Transfer times
		fprintf(dataFile,"%lld %lf ",max_transfer,(double)max_transfer/clock_rate);
		fprintf(dataFile,"%lld %lf ",min_transfer,(double)min_transfer/clock_rate);
		fprintf(dataFile,"%lld %lf ",sum_transfer/commSize,(double)(sum_transfer/commSize)/clock_rate);
		fprintf(dataFile,"%lld %lf ",sum_transfer,(double)sum_transfer/clock_rate);
		//Copy Times
		fprintf(dataFile,"%lld %lf ",max_copy,(double)max_copy/clock_rate);
		fprintf(dataFile,"%lld %lf ",min_copy,(double)min_copy/clock_rate);
		fprintf(dataFile,"%lld %lf ",sum_copy/commSize,(double)(sum_copy/commSize)/clock_rate);
		fprintf(dataFile,"%lld %lf ",sum_copy,(double)sum_copy/clock_rate);
		//Compute times
		fprintf(dataFile,"%lld %lf ",max_compute,(double)max_compute/clock_rate);
		fprintf(dataFile,"%lld %lf ",min_compute,(double)min_compute/clock_rate);
		fprintf(dataFile,"%lld %lf ",sum_compute/commSize,(double)(sum_compute/commSize)/clock_rate);
		fprintf(dataFile,"%lld %lf\n",sum_compute,(double)sum_compute/clock_rate);
	}
	MPI_Finalize();

	return 0;
}
/***********************************************************************/
/* End																   */
/* Main 									   						   */
/***********************************************************************/




















