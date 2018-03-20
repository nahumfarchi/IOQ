uniform sampler2D field;
uniform int N;

attribute vec3 LX;
attribute vec3 LY;

varying vec3 K3Screen;
varying float blendWeight;

const float twoPI = 6.2831853;

void main(void){
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
	vec2 texC = vec2(gl_Position);
	texC += 1.0;
	texC *= .5;
	
	float rotAngle = twoPI/N;
	float ca = cos(rotAngle);
	float sa = sin(rotAngle);
	
	mat3 zRot = mat3(ca, sa, 0.,
		        -sa, ca, 0.,
		         0., 0., 1.);  // this is actually the transpose
			   	      // since mat3 is assigned column major

	vec3 LZ = normalize(gl_NormalMatrix * gl_Normal);
	vec3 LX = gl_NormalMatrix * LX;
	vec3 LY = gl_NormalMatrix * LY;
	mat3 orth = mat3(LX, LY, LZ);
	mat3 orthT = mat3(LX.x, LY.x, LZ.x,
			  LX.y, LY.y, LZ.y,
			  LX.z, LY.z, LZ.z); // this is actually the transpose
			   	             // since mat3 is assigned column major

	mat3 rotMat = orth * zRot * orthT;
	
	vec4 texRGBA = vec4(texture2D(field, texC));
	
	K3Screen = vec3(texRGBA);
	K3Screen *= 2.0;
	K3Screen -= 1.0;
	K3Screen = rotMat*K3Screen;
	
	blendWeight = texRGBA.a;
}
