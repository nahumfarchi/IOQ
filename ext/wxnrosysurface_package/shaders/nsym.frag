varying vec3 worldCoords;
varying vec3 eyeCoords;
varying vec3 normal;
uniform sampler2D noise;

void main (void){
	//vec3 light0 = normalize(gl_LightSource[0].position - worldCoords);
	vec3 light1 = normalize(gl_LightSource[1].position - worldCoords);
	float ambient = .4;
	float intensity = ambient + dot( gl_LightSource[1].diffuse * light1, normal);
	int len = 20;
	float oneOverNum = 1./(len+len+1.);

	vec2 st = vec2(eyeCoords);
	
	vec2 tmp = normalize(st - vec2(.5, .5));
	vec2 v = vec2(tmp.y, -tmp.x);
	v *= 2./512.;
	vec3 rgb = vec3(texture2D(noise, st));
	
	for(int i = 0; i < len; i++){
		st += v;
		rgb += vec3(texture2D(noise, st));
	}

	st = vec2(eyeCoords);

	for(int i = 0; i < len; i++){
		st -= v;
		rgb += vec3(texture2D(noise, st));
	}

	rgb *= oneOverNum;
	
	gl_FragColor = vec4(intensity * rgb, 1.0);
}

