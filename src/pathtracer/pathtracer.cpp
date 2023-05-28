
#include "pathtracer.h"

#include "scene/light.h"
#include "scene/sphere.h"
#include "scene/triangle.h"


using namespace CGL::SceneObjects;

namespace CGL {

PathTracer::PathTracer() {
  gridSampler = new UniformGridSampler2D();
  hemisphereSampler = new UniformHemisphereSampler3D();

  tm_gamma = 2.2f;
  tm_level = 1.0f;
  tm_key = 0.18;
  tm_wht = 5.0f;
}

PathTracer::~PathTracer() {
  delete gridSampler;
  delete hemisphereSampler;
}

void PathTracer::set_frame_size(size_t width, size_t height) {
  sampleBuffer.resize(width, height);
  sampleCountBuffer.resize(width * height);
}

void PathTracer::clear() {
  bvh = NULL;
  scene = NULL;
  camera = NULL;
  sampleBuffer.clear();
  sampleCountBuffer.clear();
  sampleBuffer.resize(0, 0);
  sampleCountBuffer.resize(0, 0);
}

void PathTracer::write_to_framebuffer(ImageBuffer &framebuffer, size_t x0,
                                      size_t y0, size_t x1, size_t y1) {
  sampleBuffer.toColor(framebuffer, x0, y0, x1, y1);
}

Vector3D
PathTracer::estimate_direct_lighting_hemisphere(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // For this function, sample uniformly in a hemisphere.

  // Note: When comparing Cornel Box (CBxxx.dae) results to importance sampling, you may find the "glow" around the light source is gone.
  // This is totally fine: the area lights in importance sampling has directionality, however in hemisphere sampling we don't model this behaviour.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);

  // This is the same number of total samples as
  // estimate_direct_lighting_importance (outside of delta lights). We keep the
  // same number of samples for clarity of comparison.
  int num_samples = scene->lights.size() * ns_area_light;
  Vector3D L_out;

  // TODO (Part 3): Write your sampling loop here
  // TODO BEFORE YOU BEGIN
  // UPDATE `est_radiance_global_illumination` to return direct lighting instead of normal shading 
  for (int j = 0; j < num_samples; j++) {
    Vector3D w_j = hemisphereSampler->get_sample();
    
    double cos_theta_j = cos_theta(w_j);
    Ray newRay = Ray(hit_p, o2w * w_j);
    newRay.min_t = EPS_F;
    Intersection iSource;

    if (bvh->intersect(newRay, &iSource)) {
      Vector3D L_in = iSource.bsdf->get_emission(),
               f_r = isect.bsdf->f(w_out, w_j);
      L_out += f_r * L_in * cos_theta_j * (2 * PI);
    }
  }

  return L_out / num_samples;
}

Vector3D
PathTracer::estimate_direct_lighting_importance(const Ray &r,
                                                const Intersection &isect) {
  // Estimate the lighting from this intersection coming directly from a light.
  // To implement importance sampling, sample only from lights, not uniformly in
  // a hemisphere.

  // make a coordinate system for a hit point
  // with N aligned with the Z direction.
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  // w_out points towards the source of the ray (e.g.,
  // toward the camera if this is a primary ray)
  const Vector3D hit_p = r.o + r.d * isect.t;
  const Vector3D w_out = w2o * (-r.d);
  Vector3D L_out;

  for (auto light : scene->lights) {
    int num_samples = light->is_delta_light() ? 1 : ns_area_light;
    Vector3D L_this_light = Vector3D(0, 0, 0);

    for (int j = 0; j < num_samples; j++) {
      Vector3D w_j;
      double distToLight;
      double pdf;
      Vector3D L_in = light->sample_L(hit_p, &w_j, &distToLight, &pdf); // w_j is world space
      
      double cos_theta_j = dot(w_j, isect.n);
      Ray shadowRay = Ray(hit_p, w_j);
      shadowRay.min_t = EPS_F;
      shadowRay.max_t = distToLight - EPS_F;
      Intersection iSource;

      if (!bvh->intersect(shadowRay, &iSource)) {
        Vector3D f_r = isect.bsdf->f(w_out, w2o * w_j);
        L_this_light += f_r * L_in * cos_theta_j / pdf;
      }
    }

    L_out += L_this_light / num_samples;
  }

  return L_out;
}

Vector3D PathTracer::zero_bounce_radiance(const Ray &r,
                                          const Intersection &isect) {
  // TODO: Part 3, Task 2
  // Returns the light that results from no bounces of light
  //cout << "fucklgu " << isect.bsdf->get_emission() << endl;
  return isect.bsdf->get_emission();
}

Vector3D PathTracer::one_bounce_radiance(const Ray &r,
                                         const Intersection &isect) {
  // TODO: Part 3, Task 3
  // Returns either the direct illumination by hemisphere or importance sampling
  // depending on `direct_hemisphere_sample`
  
  if (direct_hemisphere_sample)
    return estimate_direct_lighting_hemisphere(r, isect);
  else
    return estimate_direct_lighting_importance(r, isect);
}

Vector3D PathTracer::at_least_one_bounce_radiance(const Ray &r,
                                                  const Intersection &isect) {
  Matrix3x3 o2w;
  make_coord_space(o2w, isect.n);
  Matrix3x3 w2o = o2w.T();

  Vector3D hit_p = r.o + r.d * isect.t;
  Vector3D w_out = w2o * (-r.d);

  Vector3D L_out(0, 0, 0);

  L_out += one_bounce_radiance(r, isect);

  // TODO: Part 4, Task 2
  // Returns the one bounce radiance + radiance from extra bounces at this point.
  // Should be called recursively to simulate extra bounces.
  double roulette_threshold = 1.0/3.0; // probability of termination

  if (coin_flip(roulette_threshold) || r.depth == 1 || PART < 4)
    return L_out;

  Vector3D w_j;
  double pdf;
  Vector3D f = isect.bsdf->sample_f(w_out, &w_j, &pdf);

  Ray newRay = Ray(hit_p, o2w * w_j);
  newRay.depth = r.depth - 1;
  newRay.min_t = EPS_F;

  double cos_theta_j = cos_theta(w_j);// w_j is object space

  Intersection iSource;
  if (bvh->intersect(newRay, &iSource)) {
    Vector3D L_j = at_least_one_bounce_radiance(newRay, iSource);
    L_out += f * L_j * cos_theta_j / pdf / (1-roulette_threshold);
  }

  return L_out;
}

Vector3D PathTracer::est_radiance_global_illumination(const Ray &r) {
  Intersection isect;
  Vector3D L_out;

  // You will extend this in assignment 3-2.
  // If no intersection occurs, we simply return black.
  // This changes if you implement hemispherical lighting for extra credit.

  // The following line of code returns a debug color depending
  // on whether ray intersection with triangles or spheres has
  // been implemented.
  //

  if (!bvh->intersect(r, &isect))
    return envLight ? envLight->sample_dir(r) : L_out;

  L_out = (isect.t == INF_D) ? debug_shading(r.d) : normal_shading(isect.n);
  
  // TODO (Part 3): Return the direct illumination.
  if (PART >= 3)
    L_out = zero_bounce_radiance(r, isect);
  // TODO (Part 4): Accumulate the "direct" and "indirect"
  if (max_ray_depth >= 1 && PART >= 3)
    L_out += at_least_one_bounce_radiance(r, isect);
  // parts of global illumination into L_out rather than just direct

  return L_out;
}

void PathTracer::raytrace_pixel(size_t x, size_t y) {
  // TODO (Part 1.2):
  // Make a loop that generates num_samples camera rays and traces them
  // through the scene. Return the average Vector3D.
  // You should call est_radiance_global_illumination in this function.

  // TODO (Part 5):
  // Modify your implementation to include adaptive sampling.
  // Use the command line parameters "samplesPerBatch" and "maxTolerance"

  int num_samples = ns_aa;          // total samples to evaluate
  int batch_accu = 0;
  Vector2D origin = Vector2D(x, y); // bottom left corner of the pixel
  Vector3D resColor = Vector3D(0, 0, 0); // The result color of the given pixel
  double s1 = 0, s2 = 0;
  
  int n = 1;
  for (; n <= num_samples; n++) {
    Vector2D randomPos = gridSampler->get_sample() + origin;
    Ray sampleRay = camera->generate_ray(1.0 * randomPos.x / sampleBuffer.w,
                                         1.0 * randomPos.y / sampleBuffer.h);
    sampleRay.depth = max_ray_depth;
    Vector3D radiance = est_radiance_global_illumination(sampleRay);
    resColor += radiance;
    
    if (PART == 5) {
      double illum = radiance.illum();
      s1 += illum;
      s2 += illum * illum;

      batch_accu++;
      if (batch_accu == samplesPerBatch) {
        batch_accu = 0;

        double illum_mean = s1 / n,
              illum_variance = (1.0 / (n-1)) * (s2 - s1 * s1 / n),
              I = 1.96 * sqrt(illum_variance / n);

        if (I <= maxTolerance * illum_mean)
          break;
      }
    }
  }

  sampleBuffer.update_pixel(resColor / n, x, y);
  sampleCountBuffer[x + y * sampleBuffer.w] = n;
}

void PathTracer::autofocus(Vector2D loc) {
  Ray r = camera->generate_ray(loc.x / sampleBuffer.w, loc.y / sampleBuffer.h);
  Intersection isect;

  bvh->intersect(r, &isect);

  camera->focalDistance = isect.t;
}

} // namespace CGL

