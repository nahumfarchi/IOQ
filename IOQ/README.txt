Compiling the CUDA kernel on windows:

1. Install cuda drivers
2. Add cl.exe to your path (C:\Program Files\Microsoft Visual Studio 10.0\VC\bin)
3. To compile:
	nvcc -ptx gridsearch_inner.cu