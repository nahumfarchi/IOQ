uniform sampler2D normals;
uniform sampler2D field;
uniform int N;

const float twoPI = 6.2831853;

varying vec2 texC;

// rotates the given field by 2\pi/N about the surface normal

void main (void){
	vec3 normal = vec3(texture2D(normals, texC));
	normal *= 2.0;
	normal -= 1.0;
	vec4 texFieldRGBA = vec4(texture2D(field, texC));
	vec3 K3Screen  = vec3(texFieldRGBA);
	K3Screen *= 2.0;
	K3Screen -= 1.0;

	// set up local frame
	vec3 LX = normal;
	if(LX.x == 1.0) LX.y = 1.0;
	else		LX.x = 1.0;
	LX = normalize(cross(normal, LX));
	vec3 LY = cross(normal, LX);
	vec3 LZ = normal;

	// construct rotation matrix
	float rotAngle = twoPI/N;
	float ca = cos(rotAngle);
	float sa = sin(rotAngle);
	
	mat3 zRot = mat3(ca, sa, 0.,
		        -sa, ca, 0.,
		         0., 0., 1.);  // this is actually the transpose
			   	      // since mat3 is assigned column major

	mat3 orth = mat3(LX, LY, LZ);
	mat3 orthT = mat3(LX.x, LY.x, LZ.x,
			  LX.y, LY.y, LZ.y,
			  LX.z, LY.z, LZ.z); // this is actually the transpose
			   	             // since mat3 is assigned column major

	mat3 rotMat = orth * zRot * orthT;


	// rotate K in tangent space and repack
	K3Screen = rotMat * K3Screen;

	K3Screen += 1.0;
	K3Screen *= .5;

	gl_FragColor = vec4(K3Screen, texFieldRGBA.a);
}


