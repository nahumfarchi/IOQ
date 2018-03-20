uniform int N;

//varying vec3 posWorldC;
varying vec3 posEyeC;
varying vec3 posTexC;
//varying float blend;


void main(void){
	
	gl_Position = gl_Vertex;//gl_ModelViewProjectionMatrix * gl_Vertex;
	posEyeC = vec3(gl_ModelViewProjectionMatrix * gl_Vertex);
	posTexC = posEyeC + 1.0;
	posTexC *= .5;
}
