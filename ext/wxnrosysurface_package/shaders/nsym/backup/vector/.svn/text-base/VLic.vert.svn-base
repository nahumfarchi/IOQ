uniform int N;

varying vec3 posWorldC;
varying vec3 posEyeC;
varying vec3 posTexC;
varying vec3 normal;
//varying float blend;


void main(void){
	posWorldC = vec3(gl_ModelViewMatrix * gl_Vertex);
	normal = normalize(vec3(gl_NormalMatrix * gl_Normal));
	
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	posEyeC = vec3(gl_Position);
	posTexC = posEyeC + 1.0;
	posTexC /= 2.0;
}
