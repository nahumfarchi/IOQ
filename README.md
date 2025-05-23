# Integer-only Cross Field Computation

Matlab code for the SIGGRAPH paper [Integer-only cross field computation (2018)](https://dl.acm.org/doi/abs/10.1145/3197517.3201375).

IOQ is an iterative algorithm for computing smooth cross fields on triangle meshes that is simple, easily parallelizable on the GPU, and finds solutions with lower energy and fewer cone singularities than state-of-the-art methods. The approach is based on a formal equivalence between two formulations of the optimization problem. This equivalence allows us to eliminate the real variables and design an efficient grid search algorithm for the cone singularities. It leverages a recent graph-theoretical approximation of the resistance distance matrix of the triangle mesh to speed up the computation and enable a trade-off between the computation time and the smoothness of the output.

Includes:

- Implementation with and without GPU acceleration
- Mesh data structure with support for various differential operators (divergence, gradient, laplacian, etc)
- Experiment scripts
- ... and more

![Cross field](https://github.com/nahumfarchi/IOQ/blob/master/img/sample.PNG)

![Quad mesh](https://github.com/nahumfarchi/IOQ/blob/master/img/quad.png)

# Example usage:

```
FACE0 = 1; THETA0 = 0; DEGREE = 4; GVEC = [1,0,0]; SEED = 112;
USE_GPU = false;
m = Mesh('mesh.off');
V = m.V; F = m.F; nv = m.nV; ne = m.nE; nf = m.nF;

% run IOQ
[alpha, beta] = IOQ_highgenus_gpu(...
    V, F, ...
    'UseGPU', USE_GPU, ...
    'Iterations', 2000);

% run TCODS
k = [alpha; beta];
res_tc = TCODS(m, ...
    'k', k, ...
    'f0', FACE0, ...
    'theta0', THETA0, ...
    'degree', DEGREE, ...
    'CreateFField', true, ...
    'Duplicate', true, ...
    'gConstraintVec', GVEC);
[E_edges, Emiq] = per_edge_energy(res_tc);

% Plot result
figure
res_tc.draw('FaceAlpha', 1, 'EdgeAlpha', 0)
```
