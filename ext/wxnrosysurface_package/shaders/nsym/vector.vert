uniform vec2 globalDir;
uniform mat4 NRot;
uniform int N;

attribute vec3 K3Model;

varying vec2 K2Screen;

void main(void){
	//K4Screen = 
	/*float maxTest = -1.0;
	for(int i = 1; i < N; i++){
		Vector2 v(texVect[tvij0],texVect[tvij1]);
		float test = dot(v, vTest);
		if(test > maxTest){
			maxTest = test;
			initRot[i] = j;
		}	
	}
	*/

	
	K2Screen = normalize(vec2(gl_ModelViewProjectionMatrix * vec4(K3Model, 0.0)));
	
	
	gl_Position = gl_ModelViewProjectionMatrix * gl_Vertex;
}
