
#include "bvh.h"

#include "CGL/CGL.h"
#include "triangle.h"

#include <iostream>
#include <stack>

using namespace std;

namespace CGL {
namespace SceneObjects {

BVHAccel::BVHAccel(const std::vector<Primitive *> &_primitives,
                   size_t max_leaf_size) {

  primitives = std::vector<Primitive *>(_primitives);
  root = construct_bvh(primitives.begin(), primitives.end(), max_leaf_size);
}

BVHAccel::~BVHAccel() {
  if (root)
    delete root;
  primitives.clear();
}

BBox BVHAccel::get_bbox() const { return root->bb; }

void BVHAccel::draw(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->draw(c, alpha);
    }
  } else {
    draw(node->l, c, alpha);
    draw(node->r, c, alpha);
  }
}

void BVHAccel::drawOutline(BVHNode *node, const Color &c, float alpha) const {
  if (node->isLeaf()) {
    for (auto p = node->start; p != node->end; p++) {
      (*p)->drawOutline(c, alpha);
    }
  } else {
    drawOutline(node->l, c, alpha);
    drawOutline(node->r, c, alpha);
  }
}

std::vector<Primitive *>::iterator BVHAccel::quickSelect(std::vector<Primitive *>::iterator start,
                                               std::vector<Primitive *>::iterator end,
                                               int axis = 0, int left = 0, int right = 0) {
  if (end - start == 1)
    return start;

  auto t = start, pi = start+1, pj = end-1;
  while (true) {
    BBox bi = (*pi)->get_bbox(),
          bj = (*pj)->get_bbox();
    double ct = (*t)->get_bbox().centroid()[axis];
    while (pi < pj && bj.centroid()[axis] > ct) {
      pj--;
      bj = (*pj)->get_bbox();
    }
    while (pi < pj && bi.centroid()[axis] <= ct) {
      pi++;
      bi = (*pi)->get_bbox();
    }
    if (pi == pj) {
      if (bi.centroid()[axis] <= ct)
        swap(*t, *pi);
      break;
    }
    swap(*pi, *pj);
  }

  if (abs((pi-start+left) - (end-pi-1+right)) <= 1)
    return (pi-start+left) > (end-1-pi+right) ?  pi-1 : pi;
  else if ((pi-start+left) > (end-1-pi+right))
    return quickSelect(start, pi, axis, left, right+end-pi);
  else
    return quickSelect(pi+1, end, axis, left+pi-start+1,right);

  return start;
}

BVHNode *BVHAccel::construct_bvh(std::vector<Primitive *>::iterator start,
                                 std::vector<Primitive *>::iterator end,
                                 size_t max_leaf_size) {
  // TODO (Part 2.1):
  // Construct a BVH from the given vector of primitives and maximum leaf
  // size configuration. The starter code build a BVH aggregate with a
  // single leaf node (which is also the root) that encloses all the
  // primitives.

  BBox bbox;
  
  // find the overall bounding box for this node
  for (auto p = start; p != end; p++) {
    BBox bb = (*p)->get_bbox();
    bbox.expand(bb);
  }

  BVHNode *node = new BVHNode(bbox);
  node->start = start;
  node->end = end;
  if (end - start <= max_leaf_size || PART < 2) {
    return node;
  }

  double longest = max(max(bbox.extent.x, bbox.extent.y), bbox.extent.z);

  // median is the (smaller one) of the median value of the centroid coordinate along some axis
  std::vector<Primitive *>::iterator median;
  if (longest == bbox.extent.x)
    median = quickSelect(start, end, 0);
  else if (longest == bbox.extent.y)
    median = quickSelect(start, end, 1);
  else if (longest == bbox.extent.z)
    median = quickSelect(start, end, 2);

  node->l = construct_bvh(start, median+1, max_leaf_size);

  node->r = construct_bvh(median+1, end, max_leaf_size);
  
  return node;
}

bool BVHAccel::has_intersection(const Ray &ray, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  // Take note that this function has a short-circuit that the
  // Intersection version cannot, since it returns as soon as it finds
  // a hit, it doesn't actually have to find the closest hit.

  if (node->bb.intersect(ray, ray.min_t, ray.max_t)) {
    if (node->isLeaf()) {
      for (auto p = node->start; p != node->end; p++)
        if ((*p)->has_intersection(ray))
          return true;
    }
    else {
      return has_intersection(ray, node->l) || has_intersection(ray, node->r);
    }
  }

  return false;
}

bool BVHAccel::intersect(const Ray &ray, Intersection *i, BVHNode *node) const {
  // TODO (Part 2.3):
  // Fill in the intersect function.
  if (node->bb.intersect(ray, ray.min_t, ray.max_t)) {
    total_isects++;
    if (node->isLeaf()) {
      bool hit = false;
      for (auto p = node->start; p != node->end; p++) {
        total_isects++;
        hit = (*p)->intersect(ray, i) || hit;
      }
      return hit;
    }
    else {
      bool lFlag = false, rFlag = false;
      if (node->l != NULL)
        lFlag = intersect(ray, i, node->l);
      if (node->r != NULL)
        rFlag = intersect(ray, i, node->r);
      return lFlag || rFlag;
      // return intersect(ray, i, node->l) || intersect(ray, i, node->r);
    }
  }

  return false;
}

} // namespace SceneObjects
} // namespace CGL
