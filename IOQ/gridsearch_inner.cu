#include <math.h>

#define SZ 32

__global__ 
void lattice_inner_loop(
    float * out,
    float * Lp, 
    const float * dlp, 
    const float * b, 
    const size_t rows, 
    const size_t cols, 
    const size_t numel)
{
    size_t const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= rows)
        return;
    size_t const j = blockDim.y * blockIdx.y + threadIdx.y;
    if (j >= cols)
        return;
    size_t const idx = i + j * rows;
    if (idx >= numel)
        return;
    out[idx] = dlp[i] + dlp[j] + b[i] - b[j] - 2*Lp[idx];
}

__device__ float infinity(void)
{
    return 0x7f800000;
}

__device__ float minus_infinity(void)
{
    return 0xff800000;
}

/*
Find the minimum of
	b[r] - b[c] + R_rc
	r = 1...rows
	c = 1...cols
where R_rc is the resistance defined by (Lp)_rr + (Lp)_cc - 2 (Lp)_rc
*/
__global__ 
void reduce_cols(
    float * out_min,
    unsigned int * out_r,
    const float * Lp,
    const float * dlp,
    const float * b,
    unsigned int rows,
    unsigned int cols,
    unsigned int numel)
{
    unsigned int c = blockDim.y * blockIdx.y + threadIdx.y;
    if (c >= cols)
        return;
    unsigned int idx = c * rows;

    float min = infinity();
    unsigned int min_r = 0;

    for (unsigned int r = 0; r < rows; r++, idx++)
    {
        float z = dlp[r] + dlp[c] + b[r] - b[c] - 2*Lp[idx];
		//float z = dlp[r] + b[r] - 2*Lp[idx];
        if (z < min)
        {
            min = z;
            min_r = r;
        }
    }

    out_min[c] = min;
    out_r[c] = min_r;
}

/*
Find the minimum of
	b[r] - b[c] + R_rc
	r = 1...rows
	c = 1...cols
where R_rc is the resistance defined by (Lp)_rr + (Lp)_cc - 2 (Lp)_rc. Lp is the pseudo-inverse of the Laplacian.

R is a symmetric matrix given as an array of size nv*(nv-1)/2. This array was created by taking the columns of lower triangle part (not including the diagonal!)
For example:
	* * * *
	1 * * *
	2 5 * *
	3 6 8 *
will be given as 
    [1, 2, ..., 6]
*/
__global__ 
void reduce_cols_R_symmetric(
    float * out_min,
    unsigned int * out_r,
    const float * R,
    const float * b,
    unsigned int rows,
    unsigned int cols,
    unsigned int numel)
{
	float resistance;
    unsigned int c = blockDim.y * blockIdx.y + threadIdx.y;
    if (c >= cols)
        return;
    unsigned int idx;

    float min = infinity();
    unsigned int min_r = 0;

    for (unsigned int r = 0; r < rows; r++)
    {
		if (r < c) {
			idx = ((2*rows - r - 1)*r) / 2.0 + c - r - 1;
			resistance = R[idx];
		}
		else if (r > c) {
			idx = ((2*rows - c - 1)*c) / 2.0 + r - c - 1;
			resistance = R[idx];
		}
		else {
			resistance = 0;
		}
		
		float z = b[r] - b[c] + resistance;
        if (z < min)
        {
            min = z;
            min_r = r;
        }
    }

    out_min[c] = min;
    out_r[c] = min_r;
}

__global__ 
void reduce_cols_R(
    float * out_min,
    unsigned int * out_r,
    const float * R,
    const float * b,
    unsigned int rows,
    unsigned int cols,
    unsigned int numel)
{
    unsigned int c = blockDim.y * blockIdx.y + threadIdx.y;
    if (c >= cols)
        return;
    unsigned int idx = c * rows;;

    float min = infinity();
    unsigned int min_r = 0;

    for (unsigned int r = 0; r < rows; r++, idx++)
    {
        //float z = dlp[r] + dlp[c] + b[r] - b[c] - 2*Lp[idx];
		//float z = dlp[r] + b[r] - 2*Lp[idx];
		float z = b[r] - b[c] + R[idx];
        if (z < min)
        {
            min = z;
            min_r = r;
        }
    }

    out_min[c] = min;
    out_r[c] = min_r;
}

/* 
Find the minimum of 
	b_r - b_c + R_rc
	r = 1...rows
	c = 1...cols
where R_rc is the resistance matrix, approximated by 
	| Ztilde(:, r) - Ztilde(:, c) |^2
*/
__global__ 
void reduce_cols_Ztilde(
	float * out_min,
	unsigned int * out_r,
	const float * Ztilde,
	const unsigned int k,
	const float * b,
	const unsigned int rows,
	const unsigned int cols,
	const unsigned int numel)
{
	unsigned int c = blockDim.y * blockIdx.y + threadIdx.y;
    if (c >= cols)
        return;
    unsigned int idx = c * rows;
	float min = infinity();
    unsigned int min_r = 0;
	
	const unsigned int ck = c * k;
		
	for (unsigned int r = 0; r < rows; r++, idx++)
	{
		float resistance = 0;
		unsigned int zidx1 = ck;
		unsigned int zidx2 = r * k;
		//for (unsigned int i = 0; i < k; i++, zidx1++, zidx2++)
		for(; zidx1 < ck + k; zidx1++, zidx2++)
		{
			//resistance += ( Ztilde[zidx1] - Ztilde[zidx2] ) * ( Ztilde[zidx1] - Ztilde[zidx2] );
			resistance += pow( Ztilde[zidx1] - Ztilde[zidx2], 2 );
		}
			
		float z = b[r] - b[c] + resistance;
		//float z = b[r] + resistance;
		if (z < min)
		{
			min = z;
			min_r = r;
		}
	}
	
	out_min[c] = min;
    out_r[c] = min_r;
}

/* 
*Row major* pairwise Euclidean distances (squared)
input matrix in of size n x m
output matrix out of size n x n
*/
__global__ void matrix_euclidean_distance_kernel_fast(float* out, float* in, int n, int m){
	__shared__ float Ys[SZ][SZ];
	__shared__ float Xs[SZ][SZ];

	int bx = blockIdx.x, by = blockIdx.y;
	int tx = threadIdx.x, ty = threadIdx.y;

	int yBegin = by * SZ * m;
	int xBegin = bx * SZ * m;

	int yEnd = yBegin + m - 1, y, x, k, o;

	float tmp, s = 0;

	for (y = yBegin, x = xBegin;
		y <= yEnd;
		y += SZ, x += SZ){
		Ys[ty][tx] = in[y + ty * m + tx];
		Xs[tx][ty] = in[x + ty * m + tx];
		__syncthreads();

		for (k = 0; k<SZ; k++){
			tmp = Ys[ty][k] - Xs[k][tx];
			s += tmp * tmp;
		}
		__syncthreads();
	}
	o = by * SZ * n + ty * n + bx * SZ + tx;
	out[o] = s;
}

__global__ void
gpuPdist1(float *out, float *in, int n, int m)
{
	extern __shared__ float Rs[];
	float tmp, s;
	int myRow = blockIdx.x*256 + threadIdx.x;
	
	for(int r=0; r<n; r++) {
		s = 0;
		for(int i=0; i<=m/256; i++) {
			if (i*256+threadIdx.x < m)
				Rs[i*256+threadIdx.x] = in[r*m+i*256+threadIdx.x];
		}
		__syncthreads();
		
		for(int i=0; i<m && myRow<n; i++) {
			tmp = Rs[i] - in[myRow*m + i];
			s += tmp*tmp;
		}
		if (myRow < n)
			out[myRow*n+r] = s; // not sqrtf(s)
		__syncthreads();
	}
}

__global__ void
gpuPdist2(float *out, float *in, int n, int m) {
	__shared__ float Ys[SZ][SZ];
	__shared__ float Xs[SZ][SZ];
	int bx = blockIdx.x, by = blockIdx.y;
	int tx = threadIdx.x, ty = threadIdx.y;
	int yBegin = by*SZ*m;
	int xBegin = bx*SZ*m;
	int yEnd = yBegin + m - 1, y, x, k, o;
	float tmp, s = 0;
	
	for(y=yBegin, x=xBegin; y<=yEnd; y+=SZ, x+=SZ) {
		Ys[ty][tx] = in[y + ty*m + tx];
		Xs[tx][ty] = in[x + ty*m + tx];
		__syncthreads();
		
		for(k=0; k<SZ; k++) {
			tmp = Ys[ty][k] - Xs[k][tx];
			s += tmp*tmp;
		}
		__syncthreads();
	}
	o = by*SZ*n + ty*n + bx*SZ + tx;
	out[o] = s; // not sqrtf(s)
}

__global__ void 
gpuPdist_col( float* out, float* in, int rows, int cols )
{
	int i, squareeucldist = 0;
    //int c = blockDim.x * blockIdx.x + threadIdx.x; // cols
    int r = blockDim.y * blockIdx.y + threadIdx.y; // rows

	if( r < rows  ){
		for ( i = 0; i < cols; i++ ) //row-major order
            squareeucldist  += ( in[i + cols*r] - in[i + cols*r] ) * ( in[i + cols*r] - in[i + cols*r] );
		out[r] = squareeucldist;
		squareeucldist = 0;
    }
}