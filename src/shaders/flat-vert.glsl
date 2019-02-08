#version 300 es
precision highp float;

// The vertex shader used to render the background of the scene
uniform float u_Time;
uniform float u_Bikespeed;
uniform float u_Background;

in vec4 vs_Pos;
out vec2 fs_Pos;
out vec4 fs_LightVec; 

const vec4 lightPos = vec4(5, 2, 6, 1);

void main() {
  
  fs_Pos = vs_Pos.xy;
  fs_LightVec = lightPos;  // Compute the direction in which the light source lies
  gl_Position = vs_Pos;
}
