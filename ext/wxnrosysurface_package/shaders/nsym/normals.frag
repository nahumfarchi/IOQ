varying vec3 normal;

// This shader generates a normal image for a given model

void main (void){
	normal += 1.0;
	normal *= 0.5;
	gl_FragColor = vec4(normal, 1.0);
}


