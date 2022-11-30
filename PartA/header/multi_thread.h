#include <pthread.h>

// Create other necessary functions here

#define nt 8

struct Args{
	int i;
	int N;
	int *matA; 
	int *matB; 
	int *output;
	Args(int i1,int N1,int *matA1,int *matB1,int *output1){   //Constructor for Struct
		i=i1;
		N=N1;
		matA=matA1;
		matB=matB1;
		output=output1;
	}
}*ag[nt]; 

int rowA=0;
void* singleThread1(void *ag){
  	struct Args *t1=(struct Args *)ag;  //temp struct
  	int i=t1->i;
  	int N=t1->N;
  	int *matA,*matB,*output;
  	matA=t1->matA;
  	matB=t1->matB;
  	output=t1->output;
  
  	assert( N>=4 and N == ( N &~ (N-1)));

  	union{
     		__m256i sums;
     		int sum[8];
  	};
  	int k = N/nt;  //load on each thread
  	for(int rowA=i*k;rowA<(i+1)*k;rowA+=2){
    		for(int iter=0;iter<N;iter+=2){
      			__m256i a1=_mm256_set1_epi32(matA[rowA*N + iter]+matA[(rowA+1)*N + iter]);
      			__m256i a2=_mm256_set1_epi32(matA[rowA*N+(iter+1)]+matA[(rowA+1)*N+(iter+1)]);
      			int rowC=rowA>>1;
      			int index=rowC*(N>>1);
      			
      			for(int colB=0;colB<N;colB+=8){
        			__m256i b1= _mm256_loadu_si256((__m256i*)&matB[iter*N+colB]);
        			__m256i b2= _mm256_loadu_si256((__m256i*)&matB[(iter+1)*N+colB]);
        			sums=_mm256_add_epi32(_mm256_mullo_epi32(a1,b1),_mm256_mullo_epi32(a2,b2));
        			int colC=colB>>1;
        			int C =index + colC;
        
        			output[C]+=sum[0]+sum[1];        
        			output[C+1]+=sum[2]+sum[3];        
        			output[C+2]+=sum[4]+sum[5];        
        			output[C+3]+=sum[6]+sum[7];
      			}
    		}
  	}
  /*unoptimised
   int k = N/nt;  //load on each thread
  for(int rowA=i*k;rowA<(i+1)*k;rowA+=2){
    for(int colB = 0; colB < N; colB += 2){
      int sum = 0;
      for(int iter = 0; iter < N; iter++) 
      {
      	//Initial Code
        sum += matA[rowA * N + iter] * matB[iter * N + colB];
        sum += matA[(rowA+1) * N + iter] * matB[iter * N + colB];
	sum += matA[rowA * N + iter] * matB[iter * N + (colB+1)];
        sum += matA[(rowA+1) * N + iter] * matB[iter * N + (colB+1)];
        
        //1342
        sum += matA[rowA * N + iter] * matB[iter * N + colB];
        sum += matA[rowA * N + iter] * matB[iter * N + (colB+1)];
        sum += matA[(rowA+1) * N + iter] * matB[iter * N + (colB+1)];
        sum += matA[(rowA+1) * N + iter] * matB[iter * N + colB];
        
      }

      // compute output indices
      int rowC = rowA>>1;
      int colC = colB>>1;
      int indexC = rowC * (N>>1) + colC;
      output[indexC] = sum;
    }
  }
  */
  //unopt ends
  
  /*loop interchange
  int k = N/nt;  //load on each thread
  for(int rowA=i*k;rowA<(i+1)*k;rowA+=2){
        int rowC = rowA>>1;
        for(int iter = 0; iter < N; iter++){
            	for(int colB = 0; colB < N; colB+=2) {
                	int sum = 0;
                	sum += matA[rowA * N + iter] * matB[iter * N + colB];
                	sum += matA[(rowA+1) * N + iter] * matB[iter * N + colB];
                	sum += matA[rowA * N + iter] * matB[iter * N + (colB+1)];
                	sum += matA[(rowA+1) * N + iter] * matB[iter * N + (colB+1)];
        
                	// compute output indices
                	int colC = colB>>1;
                	int indexC = rowC * (N>>1) + colC;
                	output[indexC] += sum;
            	}
    	}
  }*/
  	return 0;	  
}

// Fill in this function
void multiThread(int N, int *matA, int *matB, int *output){
	//int nt=8;
	pthread_t thread[nt];
	
	for(int i=0;i<nt;i++){
		ag[i] = new Args(i,N,matA,matB,output);  //args for thread i  //already defined above
		/*ag[i]->i=i;
		ag[i]->N=N;
		ag[i]->matA=matA;
		ag[i]->matB=matB;
		ag[i]->output=output;
		*/
		int iret=pthread_create(&thread[i],NULL,&singleThread1,(void *)ag[i]);
	}   
	for(int i=0;i<nt;i++){
		pthread_join(thread[i],NULL);
	}
}
