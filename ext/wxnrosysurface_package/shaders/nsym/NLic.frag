uniform sampler2D field;
uniform sampler2D noise;
uniform sampler2D normals;

varying vec3 posWorldC;
varying vec3 posEyeC;
varying vec3 posTexC;
varying vec3 normal;

// This shader does the actual LIC tracing and convolution

void main (void){
	int res = 512;
	float texScale = 1.5;
	int len = 40;
	
	vec2 st0 = vec2(posTexC);
	vec2 st;
	vec3 color;
	
	float oneOverNum = 1./(len+len+1.);	
	vec2 v;
	vec3 K3Screen;
	float dmax = 1./res;
	
	color = vec3(texture2D(noise,  texScale*st0));

	st = st0;
	
	for(int i = 0; i < len; i++){
		K3Screen = vec3(texture2D(field, st));
		K3Screen *= 2.0;
		K3Screen -= 1.0;

		v = vec2(K3Screen);

		v *= dmax;
		st += v;
		color += vec3(texture2D(noise,  texScale*st));
	}

	st = st0;

	for(int i = 0; i < len; i++){
		K3Screen = vec3(texture2D(field, st));
		K3Screen *= 2.0;
		K3Screen -= 1.0;

		v = vec2(K3Screen);

		v *= dmax;
		st -= vec2(v);
		color += vec3(texture2D(noise,  texScale*st));
	}

	color *= oneOverNum;

	gl_FragColor = vec4(color, 1.0);
}

