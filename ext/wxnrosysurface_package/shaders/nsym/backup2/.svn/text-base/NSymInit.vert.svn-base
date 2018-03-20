uniform vec2 globalDir;
uniform int N;

attribute vec2 K;
attribute vec3 LX;
attribute vec3 LY;

varying vec3 K3Screen;
varying float blendWeight;

const float twoPI = 6.2831853;

void main(void){
	float rotAngle = twoPI/N;
	float ca = cos(rotAngle);
	float sa = sin(rotAngle);
	mat2 rotMat = mat2(ca, -sa, sa, ca);
	
	vec4 K4Model = vec4(K.x*LX + K.y*LY, 0.0);
	vec3 v = normalize(gl_NormalMatrix * vec3(K.x*LX + K.y*LY));
	vec2 v2 = normalize(vec2(v));

	vec3 maxV = v;
	float maxTest = dot(v2, globalDir);
	for(int i = 1; i < N; i++){
		K = rotMat * K;
		vec4 K4Model = vec4(K.x*LX + K.y*LY, 0.0);
		v = normalize(gl_NormalMatrix*  vec3(K.x*LX + K.y*LY));
		v2 = normalize(vec2(v));

		float test = dot(v2, globalDir);
		if(test > maxTest){
			maxTest = test;
			maxV = v;
		}
	}
	
	K3Screen = maxV;
	float theta = atan(K3Screen.y, K3Screen.x);
	float value = cos(N * theta / 2);
	blendWeight = value*value;

	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
