uniform int N;

varying vec2 K2Screen;
varying float blendWeight;

const float PI = 3.14159265;

void main (void){
	float twoPI = 2.0*PI;
	float angle = atan(K2Screen.y, K2Screen.x);
	angle = angle >= 0.0 ? angle : angle + twoPI;
	angle /= twoPI;
	
	int sup256_3 = 16777216;
	int sup256_2 = 65536;
	int iAngle = int(angle*(sup256_3-1));

	vec3 angleVec;
	angleVec.r = float(iAngle/sup256_2);
	angleVec.g = float(mod(iAngle,sup256_2)/256);
	angleVec.b = float(mod(mod(iAngle,sup256_2),256));
	angleVec /= 255.;

	//gl_FragColor = vec4(angle1, 0.0, 0.0, 1.0);
	gl_FragColor = vec4(angleVec, blendWeight);
	//gl_FragColor = vec4(angle, 0.0, 0.0, 1.0);
	//gl_FragColor = vec4(color, 1.0);
}


