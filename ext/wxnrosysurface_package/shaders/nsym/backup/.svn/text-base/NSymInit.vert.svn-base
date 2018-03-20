uniform vec2 globalDir;
uniform int N;

attribute vec2 K;
//attribute vec3 K3D;
attribute vec3 LX;
attribute vec3 LY;


varying vec2 K2Screen;
varying vec3 color;
varying float blendWeight;

const float twoPI = 6.2831853;

void main(void){
	float rotAngle = twoPI/N;
	float ca = cos(rotAngle);
	float sa = sin(rotAngle);
	mat2 rotMat = mat2(ca, -sa, sa, ca);
	
	vec4 K4Model = vec4(K.x*LX + K.y*LY, 0.0);
	vec2 v = normalize(vec2(gl_ModelViewProjectionMatrix*K4Model));

	//v = normalize(vec2(gl_NormalMatrix*K4Model));
	vec2 maxV = v;
	float maxTest = dot(v, globalDir);
	for(int i = 1; i < N; i++){
		K = rotMat * K;
		vec4 K4Model = vec4(K.x*LX + K.y*LY, 0.0);
		v = normalize(vec2(gl_ModelViewProjectionMatrix*K4Model));
		//v = normalize(vec2(gl_NormalMatrix*K4Model));
		float test = dot(v, globalDir);
		if(test > maxTest){
			maxTest = test;
			maxV = v;
		}
	}
	
	K2Screen = maxV;
	float theta = atan(K2Screen.y, K2Screen.x);
	float value = cos(N * theta / 2);
	blendWeight = value*value;
//	color = LX;
	//K2Screen = vec2(dot(LX, LY), dot(LX, LY));
	
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
