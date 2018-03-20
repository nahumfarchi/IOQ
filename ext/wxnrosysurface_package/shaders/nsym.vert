varying vec3 worldCoords;
varying vec3 eyeCoords;
varying vec3 normal;


void main(void){
	worldCoords = vec3(gl_ModelViewMatrix * gl_Vertex);
	normal = normalize(vec3(gl_NormalMatrix * gl_Normal));
	
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	eyeCoords = vec3(gl_Position);
}

