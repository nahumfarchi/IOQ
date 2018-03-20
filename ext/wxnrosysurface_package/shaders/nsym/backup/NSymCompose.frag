uniform sampler2D img0;
uniform sampler2D img1;
uniform sampler2D blend;

varying vec2 posTexC;
varying vec3 posWorldC;
varying vec3 normal;

void main (void){
	vec3 light1 = normalize(gl_LightSource[1].position - posWorldC);
	float intensity = gl_LightSource[1].ambient + dot(gl_LightSource[1].diffuse * light1, normal);

	vec3 color0 = vec3(texture2D(img0, posTexC));
	vec3 color1 = vec3(texture2D(img1, posTexC));

	
	
	float weight = texture2D(blend, posTexC).a;
	gl_FragColor = vec4(intensity*mix(color1, color0, weight), 1.0);
//	gl_FragColor = vec4(weight, 0.0, 0.0, 1.0);
}


