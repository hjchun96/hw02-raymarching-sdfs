#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;
uniform float u_Bikespeed;
uniform float u_Lighting;

in vec2 fs_Pos;
in vec4 fs_LightVec;
out vec4 out_Col;

// Cube Struct for BBox
struct Cube {
	vec3 min;
	vec3 max;
};

// ------- Constants ------- //
const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;
int IS_PLANE;
float plane_dist;

// int IS_BLOCK;
// float block_dist[30];

// ------- Primatives ------- //
// All taken from IQ: http://www.iquilezles.org/www/articles/distfunctions/distfunctions.htm

float sdTorus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float sdSphere( vec3 p, float s )
{
  return length(p)-s;
}

float sdRoundedCylinder( vec3 p, float ra, float rb, float h )
{
    vec2 d = vec2( length(p.xz)-2.0*ra+rb, abs(p.y) - h );
    return min(max(d.x,d.y),0.0) + length(max(d,0.0)) - rb;
}

float sdRoundBox( vec3 p, vec3 b, float r )
{
  vec3 d = abs(p) - b;
  return length(max(d,0.0)) - r
         + min(max(d.x,max(d.y,d.z)),0.0); // remove this line for an only partially signed sdf 
}

float sdBox( vec3 p, vec3 b )
{
  vec3 d = abs(p) - b;
  return length(max(d,0.0))
         + min(max(d.x,max(d.y,d.z)),0.0); // remove this line for an only partially signed sdf 
}

float sdCapsule( vec3 p, vec3 a, vec3 b, float r )
{
    vec3 pa = p - a, ba = b - a;
    float h = clamp( dot(pa,ba)/dot(ba,ba), 0.0, 1.0 );
    return length( pa - ba*h ) - r;
}

float sdPlane( vec3 p, vec4 n )
{
  // n must be normalized
  return dot(p,n.xyz) + n.w;
}

// ------- Rotate Operations ------- //
mat3 rotateX(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(1, 0, 0),
        vec3(0, c, -s),
        vec3(0, s, c)
    );
}

mat3 rotateY(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, 0, s),
        vec3(0, 1, 0),
        vec3(-s, 0, c)
    );
}

mat3 rotateZ(float theta) {
    float c = cos(theta);
    float s = sin(theta);
    return mat3(
        vec3(c, -s, 0),
        vec3(s, c, 0),
        vec3(0, 0, 1)
    );
}

// ------- Operations ------- //
// All taken from IQ, constants modified appropriately

float opU( float d1, float d2 ) { return min(d1,d2); }

float opI( float d1, float d2 ) { return max(d1,d2); }

float opSmoothUnion( float d1, float d2, float k ) {
  float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
  return mix( d2, d1, h ) - k*h*(1.0-h); 
}

float opSmoothSubtraction( float d1, float d2, float k ) {
    float h = clamp( 0.5 - 0.5*(d2+d1)/k, 0.0, 1.0 );
    return mix( d2, -d1, h ) + k*h*(1.0-h); 
}

vec3 opCheapBend(vec3 p)
{
    float c = cos(0.5*p.y);
    float s = sin(0.5*p.y);
    mat2  m = mat2(c,-s,s,c);
    return  vec3(m*p.xy,p.z);
}

float opDisplace(vec3 p)
{
  float total = floor(p.x*float(u_Dimensions.x)) +
        floor(p.y*float(u_Dimensions.y));
  bool isEven = mod(total, 2.0)==0.0;
  return (isEven)? 0.4:0.1;
}

// ------- Toolbox Functions ------- //
float impulse( float k, float x )
{
    float h = k*x;
    return h*exp(1.0-h);
}


float cubicPulse( float c, float w, float x )
{
    x = abs(x - c);
    if( x > w ) return 0.0;
    x /= w;
    return 1.0 - x*x*(3.0-2.0*x);
}

float square_wave(float x, float freq, float amplitude) {
	return abs(mod(floor(x * freq), 2.0)* amplitude);
}


// ------- BVH Optimization Functions ------- //

bool rayCubeIntersect(vec3 eye, vec3 dir, Cube cube) {

	float tnear = -500.0;
	float tfar = 500.0;

	for (int i = 0; i < 3; i++) {
		if (eye[i] == 0.0) {
			if (eye[i] < cube.min[i] || eye[i] > cube.max[i]) {
				return false;
			}
		}

		float t0 = (cube.min[i] - eye[i]) / dir[i];
		float t1 = (cube.max[i] - eye[i]) / dir[i];
		if (t0 > t1) {
			float tmp = t0;
			t0 = t1;
			t1 = tmp;
		}
		tnear = max(t0, tnear);
		tfar = min(t1, tfar);
	}

	if (tnear > tfar) {
		return false;
	}
	return true;
}

Cube createBikeCube() {
	Cube cube;
	cube.min = vec3(-4.0, -5.0, -3.0);
	cube.max = vec3(3.0, 5.0, 3.0);
	return cube;
}

Cube createBlockCube() {
	Cube cube;
	cube.min = vec3(-30.0, -1.0, -30.0);
	cube.max =  vec3(30.0, 1.0, 30.0);
	return cube;
}


// ------- SDF & Ray Marching Core Functions ------- //

vec3 preProcessPnt(vec3 pnt, mat3 rotation) {
	vec3 new_pnt = rotation * pnt;
  vec3 centerOffset = vec3(-2.0, 0.0, 0.0);
  return new_pnt - centerOffset;
}

float sceneSDF(vec3 og_pnt) {

	vec3 pnt = preProcessPnt(og_pnt, rotateX(1.57));
	float timeVar =sin(u_Bikespeed * 0.3 * u_Time); 
	vec3 timeOffset = vec3(0.0, 0.0, timeVar * 0.1);
	float wheelSizeOffset = u_Bikespeed * timeVar * -0.2;

  // Define Components and Position
  float neck = sdRoundBox(rotateY(-0.3) * opCheapBend(pnt) - timeOffset - vec3(0.0, 0.0, 0.5), vec3(0.8, 0.2, 0.0), 0.3);
  float wheel1 = sdTorus(pnt + timeOffset - vec3(4.0, 0.0, 3.2), vec2(1.0 + wheelSizeOffset, 0.2));
  float pipe1 = sdCapsule(pnt + timeOffset - vec3(2.0, 0.0, 1.0), vec3(2.0, 0.0, 2.2), vec3(1.0, 0.0, 0.0), 0.2);
  float seatBase = sdRoundBox(pnt - timeOffset - vec3(2.0, 0.0, 1.0), vec3(1.0, 0.3, 0.0), 0.7);
  float seatSurr = sdRoundBox((pnt- timeOffset - vec3(2.1, 0.0, -0.3)) * rotateX(1.57), vec3(0.2, 0.3, 0.4), 0.8);									
  float wheel2 = sdTorus(pnt  + timeOffset - vec3(0.0, 0.0, 3.2), vec2(1.0 +  wheelSizeOffset, 0.2));
  float pipe2 = sdCapsule(pnt + timeOffset - vec3(2.0, 0.0, 1.0), vec3(-1.0, 0.0, 0.0), vec3(-2.0, 0.0, 2.2), 0.2);

  // Combine components using Operations
  float headlight = opI(sdBox(pnt - timeOffset + vec3(1.7, 0.0, 0.0), vec3(0.6, 0.6, 0.6)),
  											sdSphere(pnt - timeOffset + vec3(2.2, 0.0, 0.0), 0.4));

  float handle = opSmoothUnion(sdSphere(pnt - timeOffset + vec3(1.3, 0.0, 0.0), 0.9), 
  														 sdRoundedCylinder(pnt - timeOffset + vec3(1.5, 0.0, 0.0), 0.1, 0.1, 1.8), 0.7);
  float seat = opSmoothSubtraction(seatSurr, seatBase, 0.2);

  float res = opSmoothUnion(seat, neck, 0.4);
  res = opU(res, handle);
  res = opU(res, headlight);
  res = opU(res, wheel1); 
  res = opU(res, wheel2);
  res = opU(res, pipe1);
  res = opU(res, pipe2);

  plane_dist = sdPlane(pnt - vec3(0.0, 0.0, 10.0), vec4(0.0, 0.0, -1.0, -1.0));

  float plane = sdPlane(pnt - vec3(0.0, 0.0, 10.0), vec4(0.0, 0.0, -1.0, -1.0));
  res = opU(res, plane);
	for (int i = 0; i < 30; i++)  {
		float block_bool = square_wave(1.0, 1.0, 5.0);
		float x_displacement = -30.0 + block_bool * 1.2 * float(i) + mod(u_Time, 10.0) * u_Bikespeed/2.0;
  	float block = sdBox(pnt - vec3(x_displacement, 0.0, 7.0) , vec3(1.0, 1.0, 0.2));
  	res = opU(res, block);
	} 

  return res;
}

float raymarch(vec3 eye, vec3 raymarchDir, float start, float end) {// eye = ray orientation

	//BVH Optimziation
	Cube objectCube = createBikeCube();
	Cube floorCube = createBlockCube();
	bool hitCube = rayCubeIntersect(eye, raymarchDir, objectCube);
	bool hitPlane = rayCubeIntersect(eye, raymarchDir, floorCube);
	if (!hitPlane && ! hitCube) return end;

	float depth = start;
	for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
		vec3 pnt = eye + depth * raymarchDir;
		float dist = sceneSDF(pnt);

		if (dist < EPSILON) { 
			if (dist  == plane_dist) {
				IS_PLANE = 1;
			} 
			return depth;
		}
		depth += dist; // Move along the view ray, spheremarch optimization
		if (depth >= end) { // abort
			return end;
		}
	}
	return end;
}

// ------- Normal ------- //
vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
        sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
        sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
        sceneSDF(vec3(p.x, p.y, p.z  + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}


// ------- Ray March Direction Calc ------- //
vec3 calculateRayMarchPoint() {
  vec3 forward = u_Ref - u_Eye;
	vec3 right = normalize(cross(u_Up, forward));
	vec3 up = normalize(cross(forward, right));
	float len = length(u_Ref - u_Eye);
	forward = normalize(forward);
	float aspect = u_Dimensions.x / u_Dimensions.y;

  float fovy = 90.0;
	float alpha = radians(fovy/2.0);
	vec3 V = up * len * tan(alpha);
	vec3 H = right * len * aspect * tan(alpha);

	float sx = 1.0 - (2.0 * gl_FragCoord.x/u_Dimensions.x);
	float sy = (2.0 * gl_FragCoord.y/u_Dimensions.y) - 1.0;
	vec3 p = u_Ref + sx * H + sy * V;
	return p;
}

// Creat striped texture on the ground
vec4 createFloorTexture(vec3 p) {

		bool use_sqfctn = false;
		float total = floor(p.x*float(u_Dimensions.x));
		bool isEven = mod(total, 10.0)==0.0;
		if (use_sqfctn) { 
			float block_bool = square_wave(p.x, 300.0, 1.0);
			isEven =  block_bool==1.0;
		}
		vec4 col1 = vec4(225.0/255.0, 225.0/255.0, 10.0/255.0, 1.0);
		vec4 col2 = vec4(136.0/255.0, 206.0/255.0, 235.0/255.0, 1.0);
		return (isEven)? col1:col2;
}

void main() {

	IS_PLANE = 0;
	float ambiance = 0.1;

	//***  Set up Ray Direction  ***//
	vec3 p = calculateRayMarchPoint();
	vec3 dir = normalize(p - u_Eye);

  //***  Ray Marching  ***//
	float dist = raymarch(u_Eye, dir, MIN_DIST, MAX_DIST);
 	
  vec3 lightDirection = normalize(vec3(fs_LightVec.xyz));

  if (dist <= MAX_DIST - EPSILON) {

  	vec4 bike_col = vec4(0.4, 0.4, 0.4, 1.0);
		vec4 floor_col = createFloorTexture(p);
		
		// Calcualte lighting
    vec3 n = estimateNormal(u_Eye + dist * dir);
    float light = dot(-fs_LightVec.xyz, n); // Lambertian Reflective Lighting

    if (u_Lighting > 0.0) {
    	out_Col = bike_col * light * ambiance;
    	if (IS_PLANE == 1) {	
				  out_Col = floor_col;
  		}
    } else {
    	out_Col = bike_col * vec4(n, 1);
    	if (IS_PLANE == 1) {
  			out_Col = floor_col * vec4(n, 1);
  		}
    }
		return;
  }

	vec3 base_col = 0.5 * (dir + vec3(1.0, 1.0, 1.0));
  out_Col = vec4(base_col, 1.0); 

}
