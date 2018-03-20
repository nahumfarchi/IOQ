varying vec3 normal;

// This shader generates a normal image for a given model

void main(void){
	normal = normalize(gl_NormalMatrix * gl_Normal);
	
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
