#include "mex.h"
#include<iostream>
#include <math.h>
#include <time.h>
#include <string.h>
#include <map>
#include <list>
using namespace std;

/*
 Given a vertex list V, and a face list F, compute a cotangent Laplacian 
 L, possibly based on Inrinsic Delaunay Triangulation, and return it. 
 Usage: L = compute_cot_laplacian(V, F). 

 At present, no checks are performed for the vertices and faces! It is assumed
 that V is a 3 x nverts array and F is a 3 x nfaces array, where nverts and
 nfaces are the number of vertices and faces respectively. V(:, i) contains
 the 3D coordinates of vertex i; F(:, i) contains the indices of vertices
 that make up face i.
*/
double myangle(double *u, double *v);
double tanalpha(double a, double b, double c);

#define V prhs[0] // vertices
#define F prhs[1] // faces
#define idt_flag prhs[2] // flag whether to do Intrinsic Delaunay Triangulation
#define L plhs[0] // output Laplacian
#define idtL plhs[1] // output IDT Laplacian (optional)

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    double *V_pr = (double *) mxGetData(V), *F_pr = (double *) mxGetData(F);
    mwIndex *F_ir = mxGetIr(F), *F_jc = mxGetJc(F);
    mwIndex *V_ir = mxGetIr(V), *V_jc = mxGetJc(V);
    int do_idt = (int) *mxGetPr(idt_flag);

    const mwSize *F_size = mxGetDimensions(F), *V_size = mxGetDimensions(V);
    int nverts = V_size[1], nfaces = F_size[1];

    // dimensionality must be 3 x nverts and 3 x nfaces
    if (V_size[0] != 3 & V_size[1] < V_size[0])
        mexErrMsgTxt("V array must be 3 x number of vertices");

    if (F_size[0] != 3 & F_size[1] < F_size[0])
        mexErrMsgTxt("F array must be 3 x number of faces");

    fprintf(stderr, "dims are %d %d; nverts %d; nfaces %d\n", F_size[0], F_size[1], nverts, nfaces);

    // Allocate space for cot Laplacian and optionally the
    // IDT-Laplacians Currently, maximum number of non-zeros hardocoded to 10*n - must be changed later
    L = mxCreateSparse(nverts, nverts, nverts*10, mxREAL);
    double *L_pr = (double *) mxGetData(L), *idtL_pr;
    mwIndex *L_ir = mxGetIr(L), *L_jc = mxGetJc(L), *idtL_ir, *idtL_jc;

    if (do_idt) {
        idtL = mxCreateSparse(nverts, nverts, nverts*10, mxREAL);
        double *idtL_pr = (double *) mxGetData(idtL);
        mwIndex *idtL_ir = mxGetIr(L), *idtL_jc = mxGetJc(idtL);
    } else {
        idtL = NULL;
    }

    // IDT code - currently disabled
    /*// from a list of faces, compute the list of unique edges
    // e2v is a 2 x E array that gives the 2 vertices forming
    // each unique edge given it's ID;
    // v2e is the converse: given two vertices, it returns the ID of the edge
    int tmp[F_size[1]*2][2];
    for (int e1 = 0; e1 < F_size[0]; e1++) {
        int e2 = (e1 % F_size[0]) + 1;
        for (int i = A_jcol[e1]; k < A_jcol[e1 + 1]; k++) {
        }
    }
    
	
	if (do_idt) {
		// compute a map from vertices to edge index
	}*/

	map<int, list<int> > ring;
	// compute the faces adjacent to each vertex
	for (int i = 0; i < nfaces; i++) {
		for (int k = 0; k < 3; k++) {
			int vertex = F_pr[i*3 + k] - 1; 
			ring[vertex].push_back(i); 
		}
	}

	int nz = 0; 
	for (int i = 0; i < nverts; i++) {
		// sorted entries of row i
		// first entry is neighbor and second entry is value
		map<int, double> row;
		row.clear();
		list<int>::iterator it;

		ring[i].sort();
		
		it = ring[i].begin();
      while (it != ring[i].end()) {
         // find the vertices on this face
			int vertices[3], start_ptr = (*it)*3;
			for (int k = 0; k< 3; k++) {
				vertices[k] = F_pr[start_ptr + k] - 1; 
			}
			it++;
			int vertex[2];

			if (vertices[0] == i) {
				vertex[0] = vertices[1]; vertex[1] = vertices[2];
			} else if (vertices[1] == i) {
				vertex[0] = vertices[0]; vertex[1] = vertices[2];
			} else if (vertices[2] == i) {
				vertex[0] = vertices[0]; vertex[1] = vertices[1];
			} else {
            fprintf(stderr, "i %d vertices = %d %d %d\n", i, vertices[0], vertices[1], vertices[2]);
				mexErrMsgTxt("Problem in face ring");
			}

			double vi[3], vj[3], vk[3];
         int start_ptr1[3]; 
         start_ptr1[0] = i*3;
         start_ptr1[1] = vertex[0]*3;
         start_ptr1[2] = vertex[1]*3;
         // coordinates of the 3 vertices
			for (int k = 0; k < 3; k++) {
				vi[k] = V_pr[start_ptr1[0] + k];		
				vj[k] = V_pr[start_ptr1[1] + k];		
				vk[k] = V_pr[start_ptr1[2] + k];		
			}

         // Naiive implementation
			/*double vki[3], vkj[3], vji[3], vjk[3];
			for (int k = 0; k < 3; k++) {
				vki[k] = vk[k] - vi[k];
				vkj[k] = vk[k] - vj[k];
				vji[k] = vj[k] - vi[k];
				vjk[k] = vj[k] - vk[k];
			}
			double alpha = myangle(vki, vkj);
			double beta = myangle(vji, vjk);

			row[vertex[0]] += 1/tan(alpha);
			row[vertex[1]] += 1/tan(beta);*/

         // More efficient - avoids tan calls
         double vik = 0, vkj = 0, vji = 0;
         for (int k = 0; k < 3; k++) {
            double tmp = (vi[k] - vk[k]);
            vik += tmp*tmp;
            tmp = (vk[k] - vj[k]);
            vkj += tmp*tmp;
            tmp = (vj[k] - vi[k]);
            vji += tmp*tmp;
         }
         vik = sqrt(vik);
         vkj = sqrt(vkj);
         vji = sqrt(vji);

         // alpha is angle opposite i-j; beta is angle opposite i-k
         double tan_alpha2 = tanalpha(vji, vik, vkj);
         double tan_beta2 = tanalpha(vik, vkj, vji);
         double cot_alpha = (1 - tan_alpha2*tan_alpha2)/(2*tan_alpha2);
         double cot_beta = (1 - tan_beta2*tan_beta2)/(2*tan_beta2);

			row[vertex[0]] += cot_alpha;
			row[vertex[1]] += cot_beta;
		}

		// now add entries to matrix
      // make negative and add diagonal also
		L_jc[i] = nz;
		map<int, double>::iterator it1 = row.begin(); 
      double diag = 0;
      int diag_pos = -1;
		while (it1 != row.end()) {
         // check if we have gone past the diagonal
         if ((it1->first > i) && (diag_pos < 0)) {
            // save a space for the diagonal
            diag_pos = nz;
            nz++;
         }
			L_ir[nz] = it1->first;
			L_pr[nz] = -it1->second;
			nz++;
         diag += it1->second;
         it1++;
		}
      // check if all neighbors are less than i
      if (diag_pos < 0) {
        diag_pos = nz++;
      }
      L_ir[diag_pos] = i;
      L_pr[diag_pos] = diag;
	}
	L_jc[nverts] = nz;
}

/*
 Compute angle between two 3x1 vectors
*/
double myangle(double *u, double *v)
{
	double du = 0, dv = 0, uv = 0;
	for (int k = 0; k < 3; k++) {
		du += u[k]*u[k];
		dv += v[k]*v[k];
		uv += u[k]*v[k];
	}
	du = sqrt(du);
	dv = sqrt(dv);
	if (du < 1e-16) {
		du = 1e-16;
	}
	if (dv < 1e-16) {
		dv = 1e-16;
	}
	return acos(uv/(du * dv));
}

/*
 Given triangle with sides a, b, and c, return
 tan (alpha/2) where alpha is the angle opposite
 side with length a

 Taken from Section 2.2 of Fisher et. al.'s Intrinsic Delaunay
 Triangulation paper
*/

double tanalpha(double a, double b, double c)
{
	return sqrt(((a - b + c)*(a + b - c))/((a + b + c)*(-a + b + c)));
}

/*
 Given a quadrilateral with edges of length a, b, c, d; and edge e 
 such that (a, b, e) and (c, d, e) are triangles, compute the
 length of the flipped edge f, such that (a, d, f) and (b, c, f) are
 triangles

 Taken from Section 2.2 of Fisher et. al.'s Intrinsic Delaunay
 Triangulation paper
*/
double length_flip_edge(double a, double b, double c, double d, double e)
{
	double tanalpha2 = tanalpha(a, b, e);
	double tandelta2 = tanalpha(c, d, e);
	double tanad2 = (tanalpha2 + tandelta2)/(1 - tanalpha2*tandelta2);
	tanad2 *= tanad2;
	double cosad = (1 - tanad2)/(1 + tanad2);
	return sqrt(b*b + c*c - 2*cosad);
}
