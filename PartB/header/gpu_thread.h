__global__ void ReducedMatrixMul(const int *a, const int *b, int *c, int N) {
  	int n = N>>1;
  	int row = blockIdx.y * blockDim.y + threadIdx.y;
  	int col = blockIdx.x * blockDim.x + threadIdx.x;

  	int sum = 0;

  	for (int iter = 0; iter < N; iter++) {
    		sum += a[(row<<1) * N + iter] * b[iter * N + (col<<1)];
    		sum += a[(row<<1) * N + iter] * b[iter * N + (col<<1)+1];
    		sum += a[((row<<1)+1) * N + iter] * b[iter * N + (col<<1)];
    		sum += a[((row<<1)+1) * N + iter] * b[iter * N + (col<<1)+1];
  	}
  	c[row * n + col] = sum;
}


void gpuThread(int N, int *matA, int *matB, int *output){
	int n = N>>1;  //red_mat
    	size_t matrix_bytes = N * N * sizeof(int);
    	size_t result_bytes = n * n * sizeof(int);
    	
    	int *a, *b, *c;
  	cudaMalloc(&a, matrix_bytes);
  	cudaMalloc(&b, matrix_bytes);
  	cudaMalloc(&c, result_bytes);

  	cudaMemcpy(a, &matA, matrix_bytes, cudaMemcpyHostToDevice);
  	cudaMemcpy(b, &matB, matrix_bytes, cudaMemcpyHostToDevice);

  	int nt = 16;     //no.of threads
  	int nb = n / nt; //no.of blocks

  	dim3 threads(nt, nt);
  	dim3 blocks(nb, nb);

  	ReducedMatrixMul<<<blocks, threads>>>(a, b, c, N);

  	cudaMemcpy(&output, c, result_bytes, cudaMemcpyDeviceToHost);

  	cudaFree(a);
  	cudaFree(b);
  	cudaFree(c);
}

