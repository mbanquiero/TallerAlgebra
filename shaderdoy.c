

const vec3 lb = vec3(-10,-10,-10);
const vec3 rt = vec3(40,40,40);

float ray_aabb(vec3 org,vec3 dir)
{
    float t = -1.0;
    
    vec3 dirfrac;
	dirfrac.x = 1.0 / dir.x;
	dirfrac.y = 1.0 / dir.y;
	dirfrac.z = 1.0 / dir.z;
	// lb is the corner of AABB with minimal coordinates - left bottom, rt is maximal corner
	// r.org is origin of ray
	float t1 = (lb.x - org.x)*dirfrac.x;
    float t2 = (rt.x - org.x)*dirfrac.x;
    float t3 = (lb.y - org.y)*dirfrac.y;
    float t4 = (rt.y - org.y)*dirfrac.y;
    float t5 = (lb.z - org.z)*dirfrac.z;
    float t6 = (rt.z - org.z)*dirfrac.z;
	
    float tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
    float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));
	
    // if tmax < 0, ray (line) is intersecting AABB, but whole AABB is behing us
    if (tmax < 0.0)
    {
        t = -1.0;
    }
    else
		// if tmin > tmax, ray doesn't intersect AABB
		if (tmin > tmax)
		{
			t = -1.0;
		}
		else
		{
			t = tmin;
		}
		
		return t;
}

float ray_sphere(vec3 ray, vec3 dir, vec3 center, float radius)
{
	vec3 rc = ray-center;
	float c = dot(rc, rc) - (radius*radius);
	float b = dot(dir, rc);
	float d = b*b - c;
	float t = -b - sqrt(abs(d));
	float st = step(0.0, min(t,d));
	return mix(-1.0, t, st);
}

void main(void)
{
	
	vec2 uv = (-1.0 + 2.0*gl_FragCoord.xy / iResolution.xy) * 
		vec2(iResolution.x/iResolution.y, 1.0);
	vec3 D = normalize(vec3(uv, 1.0));
    vec3 vLookFrom = vec3(0.0,0.0,-3.0);
    
	float t = ray_sphere(vLookFrom,D , vec3(0.0,0.0,0.0),1.0);
    
    
    if(t>-1.0)
		gl_FragColor = vec4(0.0,0.0,0.0,1.0);
	else
		gl_FragColor = vec4(1.0,1.0,1.0,1.0);
	
}

const vec3 lb = vec3(-10,-10,-10);
const vec3 rt = vec3(40,40,40);

float ray_aabb(vec3 org,vec3 dir)
{
    float t = -1.0;
    
    vec3 dirfrac;
	dirfrac.x = 1.0 / dir.x;
	dirfrac.y = 1.0 / dir.y;
	dirfrac.z = 1.0 / dir.z;
	// lb is the corner of AABB with minimal coordinates - left bottom, rt is maximal corner
	// r.org is origin of ray
	float t1 = (lb.x - org.x)*dirfrac.x;
    float t2 = (rt.x - org.x)*dirfrac.x;
    float t3 = (lb.y - org.y)*dirfrac.y;
    float t4 = (rt.y - org.y)*dirfrac.y;
    float t5 = (lb.z - org.z)*dirfrac.z;
    float t6 = (rt.z - org.z)*dirfrac.z;

    float tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
    float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

    // if tmax < 0, ray (line) is intersecting AABB, but whole 

AABB is behing us
    if (tmax < 0.0)
    {
        t = -1.0;
    }
    else
	// if tmin > tmax, ray doesn't intersect AABB
    if (tmin > tmax)
    {
        t = -1.0;
    }
    else
    {
	    t = tmin;
    }
    
    return t;
}


float ray_sphere(vec3 p0,vec3 D,vec3 centro,float radio)
{
	float rta = -1.0;
	// verifico si interseta con la esfera
	vec3 aux = p0 - centro;
	float c = dot(aux,aux) - radio*radio;
	float B = 2.0*dot(D,aux);
	float disc = B*B - 4.0*c;
	if(disc>=0.0)
	{
		float t0 = (-B-sqrt(disc))/2.0;
		float t1 = (-B+sqrt(disc))/2.0;

		if(t0<0.0)
			rta = t1;
		else
		if(t0>t1)
		{
			float aux_t = t1;
			t1 = t0;
			t0 = aux_t;
		}

		// entra en la esfera por t0 y sale por t1. 
		// se supone que el punto p0 NO esta dentro de la 
		// esfera, asi que t0, t1 son ambos positivos o ambos negativos
		rta = t0;
	}
	return rta;
}



void main(void)
{
	vec2 uv = gl_FragCoord.xy / iResolution.xy;
	float x =2.0*(uv.x-1.0);
    float y = 1.0-2.0*uv.y;

    vec3 vLookFrom = vec3(0.0,0.0,2000);
    vec3 vLookAt = vec3(0,0,0);
    
    vec3 N = normalize(vLookAt-vLookFrom);
    vec3 V = normalize(N * vec3(0,1,0));
    vec3 U = V * N;
    const float fov = 3.1415 / 4.0;
    float k = 2.0*tan(fov/2.0);
	vec3 Dy = U*(k*iResolution.y/iResolution.x);
	vec3 Dx = V*k;
    
	// direccion de cada rayo

	vec3 D = normalize(N + Dx*x + Dy*y);
    
	float t = ray_sphere(vLookFrom,D , vec3(0.0,0.0,0.0),50.0);
    
    
    if(t>0.0)
		gl_FragColor = vec4(0.0,0.0,0.0,1.0);
	else
		gl_FragColor = vec4(1.0,1.0,1.0,1.0);

}