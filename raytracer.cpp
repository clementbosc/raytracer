
#include "raytracer.h"
#include "scene_types.h"
#include "ray.h"
#include "image.h"
#include "kdtree.h"
#include <stdio.h>
#include <math.h>


/// acne_eps is a small constant used to prevent acne when computing intersection
//  or boucing (add this amount to the position before casting a new ray !
const float acne_eps = 1e-4;

bool intersectPlane(Ray *ray, Intersection *intersection, Object *plane) {

  float t;


  //On vérifie si il y a une intersection
  if(dot(ray->dir, plane->geom.plane.normal) == 0){
    return false;
  }


  t = -((dot(ray->orig, plane->geom.plane.normal) + plane->geom.plane.dist)/dot(ray->dir, plane->geom.plane.normal));

  if(t < ray->tmin || t > ray->tmax){
      return false;
  }

  ray->tmax = t;
  intersection->normal = normalize(plane->geom.plane.normal);
  intersection->position = ray->orig + t * ray->dir;
  intersection->mat = &plane->mat;

  return true;
}

bool intersectSphere(Ray *ray, Intersection *intersection, Object *sphere) {

  float t, a, b, c, delta;
  b = 2.f * (dot(ray->dir, (ray->orig - sphere->geom.sphere.center)));
  a = 1.f;
  c = (dot((ray->orig - sphere->geom.sphere.center), (ray->orig - sphere->geom.sphere.center)) - pow(sphere->geom.sphere.radius, 2));
  delta = pow(b, 2.f) - 4.f*a*c;

  if(delta > 0.f){
    //cas 2 solutions

    float x1 = (-b - sqrt(delta))/(2.f*a);
    float x2 = (-b + sqrt(delta))/(2.f*a);


    if((x1 >= ray->tmin && x1 <= ray->tmax) && (x2 >= ray->tmin && x2 <= ray->tmax)){
      if(x1 < x2){
        t = x1;
      }else{
        t = x2;
      }
    }else if (x1 >= ray->tmin && x1 <= ray->tmax){
      t = x1;
    }else if (x2 >= ray->tmin && x2 <= ray->tmax){
      t = x2;
    }else{
      return false;
    }

    ray->tmax = t;
    intersection->position = ray->orig + t * ray->dir;
    intersection->normal = normalize(intersection->position - sphere->geom.sphere.center); 
    intersection->mat = &(sphere->mat);

    return true;

  }else if(delta == 0.f){
    //cas 1 solution

    float t = -b/(2.f*a);

    if(t < ray->tmin || t > ray->tmax){
      return false;
    }

    ray->tmax = t;
    intersection->position = ray->orig + t * ray->dir;
    intersection->normal = normalize(intersection->position - sphere->geom.sphere.center); 
    intersection->mat = &(sphere->mat);

  }else if(delta < 0.f){
    return false;
  }

  return false;
}

bool intersectScene(const Scene *scene, Ray *ray, Intersection *intersection) {
  bool hasIntersection = false;
  size_t objectCount = scene->objects.size();

  for (int i = 0; i < objectCount; i++)
  {
    if(scene->objects[i]->geom.type == SPHERE){
      if(intersectSphere(ray, intersection, scene->objects[i])){
        hasIntersection = true;
      }
    }else if(scene->objects[i]->geom.type == PLANE){
      if(intersectPlane(ray, intersection, scene->objects[i])){
        hasIntersection = true;
      }
    }
  }

  return hasIntersection;
}

float RDM_chiplus(float c) {
  return (c > 0.f) ? 1.f : 0.f;
}


/* --------------------------------------------------------------------------- */
/*
 *	The following functions are coded from Cook-Torrance bsdf model description and are suitable only
 *  for rough dielectrics material (RDM. Code has been validated with Mitsuba renderer)
 */

/** Normal Distribution Function : Beckmann
 * NdotH : Norm . Half
 */
float RDM_Beckmann(float NdotH, float alpha) {
  float tanCarreO = (1.f - pow(NdotH, 2))/(pow(NdotH, 2));

  float d = RDM_chiplus(NdotH) * (exp(-tanCarreO/pow(alpha, 2))/((float)M_PI * pow(alpha, 2) * pow(NdotH, 4)));

  return d;
  //return 0.5f;

}

// Fresnel term computation. Implantation of the exact computation. we can use the Schlick approximation
// LdotH : Light . Half
float RDM_Fresnel(float LdotH, float extIOR, float intIOR) {

  float Rs, Rp, F, cosI, cosT, sinCar_t;
  cosI = LdotH;
  sinCar_t = pow(extIOR/intIOR, 2) * (1.f - pow(cosI, 2));
  if(sinCar_t > 1.f)
    return 1.f;

  cosT = sqrt(1 - sinCar_t);

  Rs = pow((extIOR * cosI - intIOR * cosT), 2)/pow((extIOR * cosI + intIOR * cosT), 2);
  Rp = pow((extIOR * cosT - intIOR * cosI), 2)/pow((extIOR * cosT + intIOR * cosI), 2);

  F = (1.f/2.f) * (Rs + Rp);

  return F;

}


// Shadowing and masking function. Linked with the NDF. Here, Smith function, suitable for Beckmann NDF
// float RDM_chiplus(float c) {
//   return (c > 0.f) ? 1.f : 0.f;
// }

// DdotH : Dir . Half
// HdotN : Half . Norm
float RDM_G1(float DdotH, float DdotN, float alpha) {
  float tan_Ox = (sqrt(1.f - pow(DdotN, 2)) / DdotN);
  float b = 1.f / (alpha * tan_Ox);
  float k = DdotH / DdotN;
  float G1;

  if(k > 0.f && b < 1.6){
    G1 = RDM_chiplus(k) * ((3.535 * b + 2.181 * pow(b, 2)) / (1.f + 2.276 * b + 2.577 * pow(b, 2)));
  }else{
    G1 = RDM_chiplus(k);
  }

  return G1;

}

// LdotH : Light . Half
// LdotN : Light . Norm
// VdotH : View . Half
// VdotN : View . Norm
float RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN, float alpha) {

  return RDM_G1(LdotH, LdotN, alpha) * RDM_G1(VdotH, VdotN, alpha);

}

// Specular term of the Cook-torrance bsdf
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdotN : View . Norm
color3 RDM_bsdf_s(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m) {

  //RDM_Beckmann(float NdotH, float alpha)
  //RDM_Fresnel(float LdotH, float extIOR, float intIOR)
  //RDM_Smith(float LdotH, float LdotN, float VdotH, float VdotN, float alpha)

  float D, F, G;
  D = RDM_Beckmann(NdotH, m->roughness);
  F = RDM_Fresnel(LdotH, m->IOR, 1);
  G = RDM_Smith(LdotH, LdotN, VdotH, VdotN, m->roughness);

  return m->specularColor * ( (D*F*G) / (4.f * LdotN * VdotN) ) ;

  //!\todo specular term of the bsdf, using D = RDB_Beckmann, F = RDM_Fresnel, G = RDM_Smith

  
}
// diffuse term of the cook torrance bsdf
color3 RDM_bsdf_d(Material *m) {
  return m->diffuseColor / (float)M_PI;

  //!\todo compute diffuse component of the bsdf

}

// The full evaluation of bsdf(wi, wo) * cos (thetai)
// LdotH : Light . Half
// NdotH : Norm . Half
// VdotH : View . Half
// LdotN : Light . Norm
// VdtoN : View . Norm
// compute bsdf * cos(Oi)
color3 RDM_bsdf(float LdotH, float NdotH, float VdotH, float LdotN, float VdotN, Material *m) {

  return RDM_bsdf_s(LdotH, NdotH, VdotH, LdotN, VdotN, m) + RDM_bsdf_d(m);

  //! \todo compute bsdf diffuse and specular term
}




/* --------------------------------------------------------------------------- */

color3 shade(vec3 n, vec3 v, vec3 l, color3 lc, Material *mat ){

  //! \todo compute bsdf, return the shaded color taking into account the
  //! lightcolor

  vec3 h = (v + l) / length(v + l);

  float LdotH, NdotH, VdotH, LdotN, VdotN;
  LdotH = dot(l, h);
  NdotH = dot(n, h);
  VdotH = dot(v, h);
  LdotN = dot(l, n);
  VdotN = dot(v, n);

  return lc * RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat) * LdotN;

  if(dot(l,n) < 0){
    return color3(0, 0, 0);
  }else{
    return lc * RDM_bsdf(LdotH, NdotH, VdotH, LdotN, VdotN, mat) * LdotN;
    //return ((mat->diffuseColor / (float)M_PI) * dot(l,n) * lc);
  }
  	    
}

//! if tree is not null, use intersectKdTree to compute the intersection instead of intersect scene
color3 trace_ray(Scene * scene, Ray *ray, KdTree *tree) {  
  color3 color_d = color3(0,0,0);
  color3 ret, color_r;
  Intersection intersection;
  size_t lightCount = scene->lights.size();

  if(intersectScene(scene, ray, &intersection)){

    //Pour chaque source de lumière
    for (int i = 0; i < lightCount; i++){
      vec3 l = normalize(scene->lights[i]->position - intersection.position);

      Ray rayOmbre;
      Intersection interShadow;

      rayInit(&rayOmbre, (intersection.position + l * acne_eps), l, 0, distance(scene->lights[i]->position, intersection.position));

      if(!intersectScene(scene, &rayOmbre, &interShadow)){
        color_d+= shade(intersection.normal, -1.f*ray->dir, l, scene->lights[i]->color, intersection.mat);
      }else{
        color_d+= color3(0,0,0);
      }
    }

    ret = color_d;

    // Ray rayReflect;
    // vec3 dirRayRelect = reflect(, normalize(intersection.normal));

    // rayInit(&rayReflect, intersection.position, dirRayRelect);
    // rayReflect.depth = 0;

    // if(rayReflect.depth < 10){
    //   color_r = trace_ray(scene, &rayReflect, tree);
    //   rayReflect.depth++;
    // }

    // ret = color_d + color_r;

    // if(rayonReflechi.depth < 10){
    //   ret += trace_ray(scene, &rayonReflechi, tree);
    //   rayonReflechi.depth++;
    // }

  }else{
    ret = scene->skyColor;
  }

  return ret;
}

void renderImage(Image *img, Scene *scene) {

  //! This function is already operational, you might modify it for antialiasing and kdtree initializaion
  float aspect = 1.f/scene->cam.aspect;
    
  KdTree *tree =  NULL;


  //! \todo initialize KdTree

  float delta_y = 1.f / (img->height * 0.5f); //! one pixel size
  vec3 dy = delta_y * aspect * scene->cam.ydir; //! one pixel step 
  vec3 ray_delta_y = (0.5f - img->height * 0.5f) / (img->height * 0.5f) * aspect * scene->cam.ydir;

  float delta_x = 1.f / (img->width * 0.5f);
  vec3 dx = delta_x * scene->cam.xdir;
  vec3 ray_delta_x = (0.5f - img->width * 0.5f) / (img->width * 0.5f) *scene->cam.xdir;
  
    
  for(size_t j=0; j<img->height; j++) {
    if(j!=0) printf("\033[A\r");
    float progress = (float)j/img->height*100.f;
    printf("progress\t[");
    int cpt = 0;
    for(cpt = 0; cpt<progress; cpt+=5) printf(".");
    for(       ; cpt<100; cpt+=5) printf(" ");
    printf("]\n");
#pragma omp parallel for
    for(size_t i=0; i<img->width; i++) {
      color3 *ptr = getPixelPtr(img, i,j);
      vec3 ray_dir = scene->cam.center + ray_delta_x + ray_delta_y + float(i)*dx + float(j)*dy;

      Ray rx;
      rayInit(&rx, scene->cam.position, normalize(ray_dir));
      *ptr = trace_ray(scene, &rx, tree);

    }
  }
}
