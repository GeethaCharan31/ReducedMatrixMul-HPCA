// Optimize this function
#include<immintrin.h>

void singleThread(int N, int *matA, int *matB, int *output){
	
  	assert( N>=4 and N == ( N &~ (N-1)));
  	union{
     		__m256i sums;
     		int sum[8];
  	};
  	for(int rowA=0;rowA<N;rowA+=2){
    		for(int i=0;i<N;i+=2){
      			__m256i a1=_mm256_set1_epi32(matA[rowA*N + i]+matA[(rowA+1)*N + i]);
      			__m256i a2=_mm256_set1_epi32(matA[rowA*N+(i+1)]+matA[(rowA+1)*N+(i+1)]);
      			int rowC=rowA>>1;
      			int index=rowC*(N>>1);
      			
      			for(int colB=0;colB<N;colB+=8){
        			__m256i b1= _mm256_loadu_si256((__m256i*)&matB[i*N+colB]);
        			__m256i b2= _mm256_loadu_si256((__m256i*)&matB[(i+1)*N+colB]);
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

  
  /* Method 1
  for(int rowA = 0; rowA < N; rowA +=2) {
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
  
  
  /* Method 2
  for(int rowA = 0; rowA < N ; rowA +=2) {
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
  }
  */
}
