R"(
#version 330 core

struct Light {
	int type; // 0 = point
	vec3 position;
	vec3 color;
};
 
struct Material {
	vec3 amb; // Ambient
	vec3 diff; // Diffuse
	vec3 spec; // Specular
	float shini; // Spec shininess
};

layout(location=0) out vec4 out_fragcolor;

in vec3 vposition;
in vec3 vnormal;
in mat4 mv_mat;

uniform Light lights[8];
uniform Material material;
uniform int num_lights;


//
//	Calculate light contribution from a point light
//
vec3 point_light( int lightIdx, vec3 normal ){
	vec3 result = material.amb * lights[lightIdx].color; // start with ambient
	vec4 lp = mv_mat * vec4(lights[lightIdx].position,1);
	vec3 l = normalize(vec3(lp) - vposition);
	float ndotl = dot(normal, l);
	if( ndotl > 0.0 ){

		// Diffuse component:
		vec3 diffuse = material.diff * lights[lightIdx].color;
		result += diffuse*ndotl;

		// Specular component:
		if( material.shini > 0 ){
			vec3 e = normalize(-vposition);
			vec3 r = normalize((2.0 * dot(l, normal) * normal) - l);
			result += pow(max(dot(r, e), 0.0), material.shini) * material.spec * lights[lightIdx].color;
		}
	}
	return result;
}


//
//	Main fragment shader
//
void main() {

	vec3 normal = normalize(vnormal);
	vec3 result = vec3(0);

	// Loop through lights and add contribution
	for( int i = 0; i < num_lights; ++i ){
		if( lights[i].type == 0 ){
			result += point_light(i, normal);
		}
	}

	out_fragcolor = clamp( vec4(result,1), 0, 1 );
}


)"
