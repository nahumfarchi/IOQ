uniform vec2 globalDir;
uniform int N;

attribute vec2 K;
attribute vec3 LX;
attribute vec3 LY;

varying vec3 K3Screen;
varying float blendWeight;
varying vec3 normal;

const float twoPI = 6.2831853;

// picks the member vector for the RoSy at each point on the surface that best
// matches the given global direction

void main(void){
	if(length(K) > 1e-14){
		float rotAngle = twoPI/N;
		float ca = cos(rotAngle);
		float sa = sin(rotAngle);
		mat2 rotMat = mat2(ca, -sa, sa, ca);
	
		K3Screen = normalize(gl_NormalMatrix *  (K.x*LX + K.y*LY));
		vec2 K2Screen = normalize(vec2(K3Screen));
		vec2 K2 = K;
		float len = length(K);
		
		
		normal = normalize(gl_NormalMatrix * gl_Normal);//normalize(vec3(gl_Normal));
		vec3 Z = vec3(0,0,1);

		mat3 viewRotMat;
		if(dot(normal, Z) < .9){

		vec3 rotAxis = cross(normal, Z);
		rotAxis = normalize(rotAxis);
		vec3 viewNormal = cross(Z, rotAxis);
		viewNormal = normalize(viewNormal);
		float viewTransAngle = acos(dot(normal, Z));
		float cva = cos(viewTransAngle);
		float sva = sin(viewTransAngle);
		

		mat3 zRot = mat3(cva, sva, 0.,
		        	-sva, cva, 0.,
		         	 0.,  0.,  1.); // this is actually the transpose
			   		        // since mat3 is assigned column major
		//mat3 zRot = mat3(1,0,0,0,1,0,0,0,1);

		mat3 orth =  mat3(viewNormal,   Z,   rotAxis);
		mat3 orthT = mat3(viewNormal.x, Z.x, rotAxis.x,
				  viewNormal.y, Z.y, rotAxis.y,
				  viewNormal.z, Z.z, rotAxis.z); // this is actually the transpose
		   			             		 // since mat3 is assigned column major

		viewRotMat = orth * zRot * orthT;
		}		
		else
		viewRotMat = mat3(1,0,0,0,1,0,0,0,1);

		vec3 maxK3Screen = K3Screen;
		vec2 maxK3ScreenView = normalize(vec2(viewRotMat*K3Screen));
		float maxTest = dot(maxK3ScreenView, globalDir);
		for(int i = 1; i < N; i++){
			K2 = rotMat * K2;
		
			K3Screen = normalize(gl_NormalMatrix * (K2.x*LX + K2.y*LY));
			K2Screen = normalize(vec2(K3Screen));

			vec2 K3ScreenView = normalize(vec2(viewRotMat*K3Screen));

			float test = dot(K3ScreenView, globalDir);
			if(test > maxTest){
				maxTest = test;
				maxK3Screen = K3Screen;
				maxK3ScreenView = K3ScreenView;
			}
		}
	
		K3Screen = maxK3Screen;
		float theta = atan(maxK3ScreenView.y, maxK3ScreenView.x);

		float value = cos(N * theta / 2);

		blendWeight = value*value;	
	}
	else{
		K3Screen = vec3(0.0, 0.0, 0.0);
		blendWeight = 0.0;
	}

	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
