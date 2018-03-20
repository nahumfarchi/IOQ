function [A, K, d0, d1, H] = tcods_gsystem(V, F)
    %function [A, K, d0, d1, H] = tcods_gsystem(V, F)
	% Return the geometry part of the tcods system A x = -K + 2 pi k 
	% (i.e., not including the singularities).
	%
	% Input:
	%	V - vertex positions
	%	F - faces
	%
	% Output:
	%	A - [d0'; H']
	%		A is (nv+2g) x ne
	%	K - [Kg; z]
	%		Kg - gaussian curvature 
	%		z - generators curvature
	%	d0, d1 - exterior derivatives on (primal) 0-forms and 1-forms.
	%		d0 is ne x nv
	%		d1 is ne x nf
	%	H - generators (column wise, ne x (2 ng)).
	
	mesh = Mesh(V, F);
	if nnz(mesh.nbE > 0)
		error('Mesh has boundary')
	end

	H = mesh.H;
	[d0, d1] = get_exterior_derivatives(mesh);
	A = [d0'; H'];

	Kg = get_gaussian_curvature(mesh);
	%z = -flatten_generators(mesh);
    z = wrapToPi(generator_angle_defects(mesh));
    %z = generator_angle_defects(mesh);
	K = [Kg; z];
end

