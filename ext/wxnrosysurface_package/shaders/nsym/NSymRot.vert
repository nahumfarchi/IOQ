varying vec2 texC;

// rotates the given field by 2\pi/N about the surface normal

void main(void){
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	texC = vec2(gl_Position);
	texC += 1.0;
	texC *= 0.5;
}
