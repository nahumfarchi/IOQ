__global__ void add2( double * v1, const double * v2 ) 
{
    int idx = threadIdx.x;
    v1[idx] += v2[idx];
}

__global__ void add_mat(double * Lp, const double * dlp, const double * b, int rows, int cols, int numel)
{
    // Which block are we?
    //size_t const globalBlockIndex = blockIdx.x + blockIdx.y * gridDim.x;
    // Which thread are we within the block?
    //size_t const localThreadIdx = threadIdx.x + blockDim.x * threadIdx.y;
    // How big is each block?
    //size_t const threadsPerBlock = blockDim.x*blockDim.y;
    // Which thread are we overall?
    //size_t const globalThreadIdx = localThreadIdx + globalBlockIndex*threadsPerBlock;

    //if (globalThreadIdx >= numel) {
    //    return;
    //}

    //A[globalThreadIdx] = A[globalThreadIdx] + B[globalThreadIdx];
    size_t const i = blockDim.x * blockIdx.x + threadIdx.x;
    if (i >= rows)
        return;
    size_t const j = blockDim.y * blockIdx.y + threadIdx.y;
    if (j >= cols)
        return;
    size_t const idx = i + j * rows;
    if (idx >= numel)
        return;
    Lp[idx] = dlp[i] + dlp[j] + b[i] - b[j] - 2*Lp[idx];
}

__global__ void lattice_inner_loop(
    float * out,
    float * Lp, 
    const float * dlp, 
    const float * b, 
    int rows, 
    int cols, 
    int numel)
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