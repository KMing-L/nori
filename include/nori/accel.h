/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/mesh.h>

NORI_NAMESPACE_BEGIN

struct AccelNode {
    uint32_t first_child = 0; // child index
    BoundingBox3f bbox;
    std::vector<std::pair<uint32_t, uint32_t>> indices;

    AccelNode() : bbox() {}

    explicit AccelNode(BoundingBox3f box) : bbox(std::move(box)) {}

    explicit AccelNode(BoundingBox3f box, uint32_t count)
        : bbox(std::move(box)), indices(count) {}
};

/**
 * \brief Acceleration data structure for ray intersection queries
 *
 * The current implementation falls back to a brute force loop
 * through the geometry.
 */
class Accel {
  public:
    /**
     * \brief Register a triangle mesh for inclusion in the acceleration
     * data structure
     *
     * This function can only be used before \ref build() is called
     */
    void addMesh(Mesh *mesh);

    /// Build the acceleration data structure (currently a no-op)
    void build();

    /// Return an axis-aligned box that bounds the scene
    const BoundingBox3f &getBoundingBox() const { return m_bbox; }

    /**
     * \brief Intersect a ray against all triangles stored in the scene and
     * return detailed intersection information
     *
     * \param ray
     *    A 3-dimensional ray data structure with minimum/maximum extent
     *    information
     *
     * \param its
     *    A detailed intersection record, which will be filled by the
     *    intersection query
     *
     * \param shadowRay
     *    \c true if this is a shadow ray query, i.e. a query that only aims to
     *    find out whether the ray is blocked or not without returning detailed
     *    intersection information.
     *
     * \return \c true if an intersection was found
     */
    bool rayIntersect(const Ray3f &ray, Intersection &its,
                      bool shadowRay) const;

    virtual void divide(uint32_t n, std::vector<AccelNode *> *children) = 0;

    virtual bool traverse(uint32_t n, Ray3f &ray, Intersection &its,
                          bool shadowRay, uint32_t &f) const = 0;

    virtual std::pair<uint32_t, uint32_t> getLimits() const = 0;

    uint32_t getTriCount() { return m_indices.size(); }

    virtual ~Accel() {
        for (int i = 0, n = m_meshes.size(); i < n; i++)
            delete m_meshes[i];
        for (int i = 0, n = m_tree.size(); i < n; i++)
            delete m_tree[i];
    }

  protected:
    std::vector<Mesh *> m_meshes; ///< Mesh (only a single one for now)
    BoundingBox3f m_bbox;         ///< Bounding box of the entire scene
    std::vector<std::pair<uint32_t, uint32_t>> m_indices;
    std::vector<AccelNode *> m_tree;

    uint32_t depth_curr_ = -1;
    uint32_t count_leaf_ = 0;

  private:
    uint32_t max_mesh_ = 10;
};

class OctTree : public Accel {
  public:
    void divide(uint32_t n, std::vector<AccelNode *> *children) override;

    bool traverse(uint32_t n, Ray3f &ray, Intersection &its, bool shadowRay,
                  uint32_t &f) const override;

    std::pair<uint32_t, uint32_t> getLimits() const override {
        return {COUNT_LEAF_TRI, LIMIT_DEPTH};
    };

  private:
    uint32_t COUNT_LEAF_TRI = 16;
    uint32_t LIMIT_DEPTH = 12;
};

class BVH : public Accel {
  public:
    void divide(uint32_t n, std::vector<AccelNode *> *children) override;

    bool traverse(uint32_t n, Ray3f &ray, Intersection &its, bool shadowRay,
                  uint32_t &f) const override;

    std::pair<uint32_t, uint32_t> getLimits() const override {
        return {COUNT_LEAF_TRI, LIMIT_DEPTH};
    };

    uint32_t getBlockCount() const { return COUNT_BLOCK; }

  private:
    uint32_t COUNT_LEAF_TRI = 16;
    uint32_t LIMIT_DEPTH = 32;
    uint32_t COUNT_BLOCK = 10;
};

NORI_NAMESPACE_END
