uniform float greenChannel;
varying vec3 world;
varying vec3 normal;


void main (void){
	vec3 light = normalize(gl_LightSource[0].position - world);
	float intensity = dot(light, normal);
	vec3 rgb = vec3(0.0, greenChannel, 0.0);
	gl_FragColor = vec4(intensity * rgb, 1.0);
}

