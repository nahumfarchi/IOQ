varying vec3 world;
varying vec3 normal;

vec3 getLight();

void main (void){
	vec3 rgb = vec3(1.0, 1.0, 1.0);
	//gl_FragColor = vec4(getLight() * rgb, 1.0);
	float light = length(getLight());
	
	
	gl_FragColor = light*gl_Color;
}

vec3 getLight(){
	vec3 lightVec0 = normalize(vec3(gl_LightSource[0].position) - world);
	vec3 lightVec1 = normalize(vec3(gl_LightSource[1].position) - world);
	float intensity0 = dot(lightVec0, normal);
	float intensity1 = dot(lightVec1, normal);

	vec3 light0 = vec3(gl_LightSource[0].ambient) + vec3(gl_LightSource[0].diffuse)*intensity0;
	vec3 light1 = vec3(gl_LightSource[1].ambient) + vec3(gl_LightSource[1].diffuse)*intensity1;

	return (light0 + light1);
}
