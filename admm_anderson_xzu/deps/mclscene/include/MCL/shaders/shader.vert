R"(
#version 330 core

layout(location=0) in vec3 in_position;
layout(location=1) in vec3 in_color;
layout(location=2) in vec3 in_normal;

out vec3 vposition;
out vec3 vnormal;
out mat4 mv_mat;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform mat4 invmv;

void main() {
	mv_mat = view * model;
	vec4 pos = vec4(in_position,1.0);
	vposition = ( mv_mat * pos ).xyz;
	vec4 mv_normal = invmv * vec4(in_normal,0.0);
	vnormal = normalize( mv_normal.xyz );
	gl_Position = projection * mv_mat * pos;
}

)"
