varying vec2 posTexC;
varying vec3 posWorldC;
varying vec3 normal;

void main(void){
	posWorldC = vec3(gl_ModelViewMatrix * gl_Vertex);
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	normal = gl_NormalMatrix * gl_Normal;
	posTexC = vec2(gl_Position);
	posTexC += 1.0;
	posTexC *= .5;
}
