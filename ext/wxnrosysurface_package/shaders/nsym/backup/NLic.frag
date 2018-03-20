uniform sampler2D field;
uniform sampler2D noise;
//uniform sampler2D drawBuffer;
//uniform int numImgs;

varying vec3 posWorldC;
varying vec3 posEyeC;
varying vec3 posTexC;
varying vec3 normal;

const float PI = 3.14159265;
const float twoPI = 6.2831853;

float getAngle(vec3 rgb);

void main (void){
	
	vec3 light1 = normalize(gl_LightSource[1].position - posWorldC);
	float intensity = gl_LightSource[1].ambient + dot(gl_LightSource[1].diffuse * light1, normal);

	vec2 st0 = vec2(posTexC);
	vec2 st;
	vec3 color;
	int len = 20;
	float oneOverNum = 1./(len+len+1.);	
	vec2 v;
	vec3 texRGB;
	float angle;
	float dmax = 2./512;
	float texScale = 3.0;
	//float weight = 1./numImgs;
	
	color = vec3(texture2D(noise,  texScale*st0));

	st = st0;
	
	for(int i = 0; i < len; i++){
		texRGB = vec3(texture2D(field, st));
		angle = getAngle(texRGB);
		v.x = cos(angle);
		v.y = sin(angle);
		v *= dmax;
		st += v;
		color += vec3(texture2D(noise,  texScale*st));
	}

	st = st0;

	for(int i = 0; i < len; i++){
		texRGB = vec3(texture2D(field, st));
		angle = getAngle(texRGB);
		v.x = cos(angle);
		v.y = sin(angle);
		v *= dmax;
		st -= vec2(v);
		color += vec3(texture2D(noise,  texScale*st));
	}

	color *= oneOverNum;
	
	//vec3 oldColor = texture2D(drawBuffer, st0);
	//gl_FragColor = vec4(oldColor + weight * intensity * color, 1.0);
	gl_FragColor = vec4(color, 1.0);
}

float getAngle(vec3 rgb){
	vec3 convert = vec3(256.*256., 256., 1.);
	return dot(convert, rgb*255)*twoPI/16777215;
}
