#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;
uniform float u_Bikespeed;
uniform float u_Background;

in vec2 fs_Pos;
in vec4 fs_LightVec;
out vec4 out_Col;


// ------- Constants ------- //
const int MAX_MARCHING_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.0001;
int IS_PLANE;


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

float opS( float d1, float d2 ) { return max(-d1,d2); }

float opI( float d1, float d2 ) { return max(d1,d2); }

float opSmoothUnion( float d1, float d2, float k ) {
  float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
  return mix( d2, d1, h ) - k*h*(1.0-h); 
}

vec3 opCheapBend(vec3 p)
{
    float c = cos(0.5*p.y);
    float s = sin(0.5*p.y);
    mat2  m = mat2(c,-s,s,c);
    return  vec3(m*p.xy,p.z);
}

vec3 opTwist( vec3 p )
{
    float  c = cos(10.0*p.y+10.0);
    float  s = sin(10.0*p.y+10.0);
    mat2   m = mat2(c,-s,s,c);
    return vec3(m*p.xz,p.y);
}

float opOnion( in float sdf, in float thickness )
{
    return abs(sdf)-thickness;
}

// float opDisplace( in sdf3d primitive, in vec3 p )
// {
//     float d1 = primitive(p);
//     float d2 = displacement(p);
//     return d1+d2;
// }


// 	float opRep( in vec3 p, in vec3 c, in sdf3d primitive )
// {
//     vec3 q = mod(p,c)-0.5*c;
//     return primitve( q );
// }


float impulse( float k, float x )
{
    float h = k*x;
    return h*exp(1.0-h);
}

float square_wave(float x, float freq, float amplitude) {
	return abs(mod(floor(x * freq), 2.0)* amplitude);
}

// ------- SDF & Ray Marching Core Functions ------- //

vec3 preProcessPnt(vec3 pnt) {
	vec3 new_pnt = rotateX(1.57) *pnt;
  vec3 centerOffset = vec3(-2.0, 0.0, 0.0);
  return new_pnt - centerOffset;
}

// This is the distance field function.  The distance field represents the closest distance to the surface
// of any object we put in the scene.  If the given point (point p) is inside of an object, we return a
// negative answer.
float sceneSDF(vec3 pnt) { // map fctn

	pnt = preProcessPnt(pnt);
	vec3 timeOffset = vec3(0.0, 0.0, sin(u_Bikespeed * 0.3 * u_Time) * 0.1);

  // Define Components and Position
  float handle = opSmoothUnion(sdSphere(pnt - timeOffset + vec3(1.5, 0.0, 0.0), 1.0), 
  														 sdRoundedCylinder(pnt - timeOffset + vec3(1.5, 0.0, 0.0), 0.1, 0.1, 2.0), 0.7);
  float neck = sdRoundBox(rotateY(-0.3) * opCheapBend(pnt) - timeOffset - vec3(0.0, 0.0, 0.5), vec3(1.0, 0.2, 0.0), 0.3);
  float wheel1 = sdTorus(pnt + timeOffset - vec3(4.0, 0.0, 3.2), vec2(1.2,0.2));
  float pipe1 = sdCapsule(pnt + timeOffset - vec3(2.0, 0.0, 1.0), vec3(2.0, 0.0, 2.2), vec3(1.0, 0.0, 0.0), 0.2);
  float seat = sdRoundBox(pnt - timeOffset - vec3(2.0, 0.0, 1.0), vec3(1.0, 0.3, 0.0), 0.7);
  float wheel2 = sdTorus(pnt  + timeOffset - vec3(0.0, 0.0, 3.2), vec2(1.2,0.2));
  float pipe2 = sdCapsule(pnt + timeOffset - vec3(2.0, 0.0, 1.0), vec3(-1.0, 0.0, 0.0), vec3(-2.0, 0.0, 2.2), 0.2);


  float res = opU(handle, wheel1); 
  res = opU(res, wheel2);
  res = opU(res, seat);
  res = opU(res, pipe1);
  res = opU(res, pipe2);
  res = opU(res, neck);

  float plane = sdPlane(pnt - vec3(0.0, 0.0, 10.0), vec4(0.0, 0.0, -1.0, -1.0));
  res = opU(res, plane);

	for (int i = 0; i < 30; i++)  {
		float block_bool = square_wave(1.0, 1.0, 5.0);
		float x_displacement = -30.0 + block_bool * 1.2 * float(i) + mod(u_Time, 10.0) * u_Bikespeed/2.0;
  	float block = sdBox(pnt - vec3(x_displacement, 0.0, 7.0) , vec3(1.0, 1.0, 0.4));
  	res = opU(res, block);
	} 

  return res;
}

float raymarch(vec3 eye, vec3 marchingDirection, float start, float end) { // aymarch fntn ray orient, ray dir

	//end-=sin(sqrt(fs_Pos.x * fs_Pos.x + fs_Pos.y * fs_Pos.y) - u_Time * 6.);/// * 0.000096;
	float depth = start;

	for (int i = 0; i < MAX_MARCHING_STEPS; i++) {
		vec3 pnt = eye + depth * marchingDirection;
		float dist = sceneSDF(pnt);

		if (dist < EPSILON) { // inside scene surface

			if (dist  == sdPlane(preProcessPnt(pnt) - vec3(0.0, 0.0, 10.0), vec4(0.0, 0.0, -1.0, -1.0))) {
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


// vec3 background(vec3 d )
// {
//     // cheap cubemap
//     vec3 n = abs(d);
//     vec2 uv = (n.x>n.y && n.x>n.z) ? d.yz/d.x: 
//               (n.y>n.x && n.y>n.z) ? d.zx/d.y:
//                                      d.xy/d.z;
    
//     // fancy blur
//     vec3  col = vec3( 0.0 );
//     for( int i=1; i<200; i++ )
//     {
//         float h = float(i)/200.0;
//         float an = 31.0*6.2831*h;
//         vec2  of = vec2( cos(an), sin(an) ) * h;

//         vec3 tmp = vec3( uv*0.25 + 0.0075*of, 4.0 ).yxz;
//         col = max(max( col, tmp), 0.5 );
//     }
    
//     return pow(col,vec3(3.5,3.0,6.0))*0.2;
// }

// float gain(float x, float k) 
// {
//     float a = 0.5*pow(2.0*((x<0.5)?x:1.0-x), k);
//     return (x<0.5)?a:1.0-a;
// }

// vec4 applyTexture(vec4 col) {
// 	return vec4(gain(col.x,0.5), gain(col.y, 0.5), gain(col.z, 0.5), 1.0);
// }

vec3 ProceduralTexture(vec2 uv){
	vec2 pt = 2.0 * uv - 1.0;
	return vec3(length(pt));
}


void main() {

	IS_PLANE = 0;

	// vec2 uv = gl_FragCoord.xy / u_Dimensions;
	// if (uv.y < 0.2) {
	// 	vec2 pt = 2.0 * uv - 1.0;
	// 	//pt.x *= u_Dimensions.x / u_Dimensions.y;
	// 	pt.x = mod(u_Dimensions.x, u_Time * 2.0) ;
	// 	pt.y += 0.75;
	// 	float inside = float(length(pt) < 0.2);

	// 	out_Col = vec4(vec3(inside), 1);
	// 	return;
	// }
	//***  Set up Ray Direction  ***//
	vec3 p = calculateRayMarchPoint();
	vec3 dir = normalize(p - u_Eye);
	// out_Col = vec4(0.5 * (dir + vec3(1.0, 1.0, 1.0)), 1);
	// return;
	// vec3 col = 0.5 * (dir + vec3(1.0, 1.0, 1.0));

  //***  Ray Marching  ***//
	float dist = raymarch(u_Eye, dir, MIN_DIST, MAX_DIST);
  if (IS_PLANE == 1) {
  	out_Col = vec4(0.2, 0.2, 0.2, 1.0);
  	return;
  }

  if (dist <= MAX_DIST - EPSILON) {
    // Hit Something!
    vec3 n = estimateNormal(u_Eye + dist * dir);
    float light = dot(-fs_LightVec.xyz, n);
    // vec4 col = vec4(128.0/255.0, 128.0/255.0, 128.0/255.0, 2.0);
    vec4 col = vec4(0.1, 0.1, 0.1, 1.0);

    //vec2 objUv = GetUVs();
    if (dist < 0.2) {

    }
    // vec4 ret = col * light;
    out_Col = vec4(n, 1);
		return;
  }


	vec3 base_col = 0.5 * (dir + vec3(1.0, 1.0, 1.0));
  // vec3 base_col = vec3(0.7, 0.6, 0.0); 
  out_Col = vec4(base_col, 1.0); 

}
