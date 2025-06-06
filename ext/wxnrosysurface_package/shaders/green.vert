varying vec3 world;
varying vec3 normal;

void main(void){
	world = vec3(gl_ModelViewMatrix * gl_Vertex);
	normal = normalize(vec3(gl_NormalMatrix * gl_Normal));
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}

