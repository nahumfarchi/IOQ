varying vec3 K3Screen;
varying float blendWeight;

// picks the member vector for the RoSy at each point on the surface that best
// matches the given global direction

void main (void){
	K3Screen += 1.0;
	K3Screen *= 0.5;

	gl_FragColor = vec4(K3Screen, blendWeight);
}


