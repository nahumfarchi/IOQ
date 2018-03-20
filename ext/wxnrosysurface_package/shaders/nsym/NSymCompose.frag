uniform sampler2D img0;
uniform sampler2D img1;
uniform sampler2D blend0;
uniform sampler2D blend1;
uniform int N;

varying vec2 posTexC;
varying vec3 posWorldC;
varying vec3 normal;
varying vec3 color;

// blends the final two compositions of LIC images to remove artifacts

vec3 getLight();

void main (void){
	vec3 light = getLight();

	if(N%2==0)
		N = N/2;

	float stretchScaleUni = sqrt(float(N));

	vec3 lum0 = vec3(texture2D(img0, posTexC)).rgb;
	vec3 lum1 = vec3(texture2D(img1, posTexC)).rgb;

	lum0 = (stretchScaleUni*(lum0 - .5)+.5);
	lum1 = (stretchScaleUni*(lum1 - .5)+.5);
	
	float weight = texture2D(blend0, posTexC).a;

	float stretchScale = 1.4;

	vec3 lum = mix(lum1, lum0, weight);
	lum = (stretchScale*(lum - .5) + .5);
	
	gl_FragColor = vec4(light*lum, 1.0);
}

vec3 getLight(){
	vec3 lightVec0 = normalize(vec3(gl_LightSource[0].position) - posWorldC);
	vec3 lightVec1 = normalize(vec3(gl_LightSource[1].position) - posWorldC);
	float intensity0 = dot(lightVec0, normal);
	float intensity1 = dot(lightVec1, normal);

//	if(intensity0 < 0)
//		intensity0 = 0;
//	if(intensity1 < 0)
//		intensity1 = 0;

	vec3 light0 = vec3(gl_LightSource[0].ambient) + vec3(gl_LightSource[0].diffuse)*intensity0;
	vec3 light1 = vec3(gl_LightSource[1].ambient) + vec3(gl_LightSource[1].diffuse)*intensity1;

	return (light0 + light1);
}


