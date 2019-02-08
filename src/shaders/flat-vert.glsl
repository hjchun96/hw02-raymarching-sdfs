#version 300 es
precision highp float;

// The vertex shader used to render the background of the scene
uniform float u_Time;
// uniform mat4 u_Model;

in vec4 vs_Pos;
out vec2 fs_Pos;
out vec4 fs_LightVec; 

const vec4 lightPos = vec4(5, 0, 3, 1);

void main() {
  
	// vec4 modelposition = u_Model * vs_Pos;
  fs_Pos = vs_Pos.xy;
  fs_LightVec = lightPos;// - modelposition;  // Compute the direction in which the light source lies
  //float zOffset = sin(u_Time/10.0) * 0.1;
  gl_Position = vs_Pos;
}
