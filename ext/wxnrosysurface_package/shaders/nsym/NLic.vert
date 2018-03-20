uniform int N;

varying vec3 posEyeC;
varying vec3 posTexC;

// This shader does the actual LIC tracing and convolution

void main(void){
	
	//gl_Position = gl_Vertex;//gl_ModelViewProjectionMatrix * gl_Vertex;
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	posEyeC = vec3(gl_ModelViewProjectionMatrix * gl_Vertex);
	posTexC = posEyeC + 1.0;
	posTexC *= .5;
}
