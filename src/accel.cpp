/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <Eigen/Geometry>
#include <chrono>
#include <nori/accel.h>
#include <queue>

NORI_NAMESPACE_BEGIN

void Accel::addMesh(Mesh *mesh) {
    if (m_meshes.size() >= max_mesh_)
        throw NoriException("Accel: only ten meshes is supported!");

    m_meshes.emplace_back(mesh);
    m_bbox.expandBy(mesh->getBoundingBox());

    uint32_t index = m_meshes.size() - 1;
    for (uint32_t i = 0, n = mesh->getTriangleCount(); i < n; i++) {
        m_indices.emplace_back(std::pair<uint32_t, uint32_t>(index, i));
    }
}

void Accel::build() {
    if (m_meshes.empty())
        return;

    cout << endl << ">>>> Begin building tree <<<<" << endl;

    auto start = std::chrono::high_resolution_clock::now();

    auto [COUNT_LEAF_TRI, LIMIT_DEPTH] = getLimits();

    AccelNode *root = new AccelNode(getBoundingBox(), getTriCount());
    root->indices = m_indices;
    m_tree.emplace_back(root);
    count_leaf_++;

    std::queue<uint32_t> q;
    q.push(0);

    std::vector<AccelNode *> children;
    while (!q.empty()) {
        depth_curr_++;
        for (uint32_t i = 0, n = q.size(); i < n; i++) {
            auto node = m_tree[q.front()];
            if (node->indices.size() > COUNT_LEAF_TRI &&
                this->depth_curr_ < LIMIT_DEPTH) {
                divide(q.front(), &children);
                node->first_child = m_tree.size();
                count_leaf_--;
                for (uint32_t i = 0, n = children.size(); i < n; i++) {
                    count_leaf_++;
                    q.push(m_tree.size());
                    m_tree.emplace_back(children[i]);
                }
            }
            q.pop();
            children.clear();
            children.shrink_to_fit();
        }
    }

    auto end = std::chrono::high_resolution_clock::now();
    auto time =
        std::chrono::duration_cast<std::chrono::milliseconds>(end - start)
            .count();

    cout << "[build time]: " << time << "ms" << endl;
    cout << "[max depth]: " << depth_curr_ << endl;
    cout << "[node count]: " << m_tree.size() << endl;
    cout << "[leaf count]: " << count_leaf_ << endl;
}

void OctTree::divide(uint32_t n, std::vector<AccelNode *> *children) {
    // cout << "divide tree[" << n
    //      << "] triangularCount = " << m_tree[n]->indices.size() << endl;
    auto node = m_tree[n];
    auto center = node->bbox.getCenter();
    for (int i = 0; i < 8; i++) {
        auto corner = node->bbox.getCorner(i);
        AccelNode *child = new AccelNode(BoundingBox3f(center));
        child->bbox.expandBy(corner);

        for (auto [index_mesh, index_indice] : node->indices) {
            if (child->bbox.overlaps(
                    m_meshes[index_mesh]->getBoundingBox(index_indice)))
                child->indices.emplace_back(
                    std::pair<uint32_t, uint32_t>(index_mesh, index_indice));
        }

        children->emplace_back(child);
    }
}

void BVH::divide(uint32_t n, std::vector<AccelNode *> *children) {
    
}

bool Accel::rayIntersect(const Ray3f &ray_, Intersection &its,
                         bool shadowRay) const {
    bool foundIntersection = false; // Was an intersection found so far?
    auto f = (uint32_t)-1; // Triangle index of the closest intersection
    Ray3f ray(ray_); /// Make a copy of the ray (we will need to update its
                     /// '.maxt' value)

    /* Brute force search through all triangles */
    // for (uint32_t idx = 0; idx < m_mesh->getTriangleCount(); ++idx) {
    //     float u, v, t;
    //     if (m_mesh->rayIntersect(idx, ray, u, v, t)) {
    //         /* An intersection was found! Can terminate
    //            immediately if this is a shadow ray query */
    //         if (shadowRay)
    //             return true;
    //         ray.maxt = its.t = t;
    //         its.uv = Point2f(u, v);
    //         its.mesh = m_mesh;
    //         f = idx;
    //         foundIntersection = true;
    //     }
    // }

    foundIntersection = this->traverse(0, ray, its, shadowRay, f);

    if (foundIntersection) {
        /* At this point, we now know that there is an intersection,
           and we know the triangle index of the closest such intersection.

           The following computes a number of additional properties which
           characterize the intersection (normals, texture coordinates, etc..)
        */

        /* Find the barycentric coordinates */
        Vector3f bary;
        bary << 1 - its.uv.sum(), its.uv;

        /* References to all relevant mesh buffers */
        const Mesh *mesh = its.mesh;
        const MatrixXf &V = mesh->getVertexPositions();
        const MatrixXf &N = mesh->getVertexNormals();
        const MatrixXf &UV = mesh->getVertexTexCoords();
        const MatrixXu &F = mesh->getIndices();

        /* Vertex indices of the triangle */
        uint32_t idx0 = F(0, f), idx1 = F(1, f), idx2 = F(2, f);

        Point3f p0 = V.col(idx0), p1 = V.col(idx1), p2 = V.col(idx2);

        /* Compute the intersection positon accurately
           using barycentric coordinates */
        its.p = bary.x() * p0 + bary.y() * p1 + bary.z() * p2;

        /* Compute proper texture coordinates if provided by the mesh */
        if (UV.size() > 0)
            its.uv = bary.x() * UV.col(idx0) + bary.y() * UV.col(idx1) +
                     bary.z() * UV.col(idx2);

        /* Compute the geometry frame */
        its.geoFrame = Frame((p1 - p0).cross(p2 - p0).normalized());

        if (N.size() > 0) {
            /* Compute the shading frame. Note that for simplicity,
               the current implementation doesn't attempt to provide
               tangents that are continuous across the surface. That
               means that this code will need to be modified to be able
               use anisotropic BRDFs, which need tangent continuity */

            its.shFrame =
                Frame((bary.x() * N.col(idx0) + bary.y() * N.col(idx1) +
                       bary.z() * N.col(idx2))
                          .normalized());
        } else {
            its.shFrame = its.geoFrame;
        }
    }

    return foundIntersection;
}

bool OctTree::traverse(uint32_t n, Ray3f &ray, Intersection &its,
                       bool shadowRay, uint32_t &f) const {
    auto node = m_tree[n];

    if (!node->bbox.rayIntersect(ray))
        return false;

    bool foundIntersection = false;
    if (!node->first_child) {
        for (auto [index_mesh, index_indice] : node->indices) {
            float u, v, t;
            if (m_meshes[index_mesh]->rayIntersect(index_indice, ray, u, v,
                                                   t)) {
                /* An intersection was found! Can terminate
                   immediately if this is a shadow ray query */
                if (shadowRay)
                    return foundIntersection;
                ray.maxt = its.t = t;
                its.uv = Point2f(u, v);
                its.mesh = m_meshes[index_mesh];
                f = index_indice;
                foundIntersection = true;
            }
        }
    } else {
        std::pair<uint32_t, float> dis[8];
        for (uint32_t i = 0; i < 8; i++) {
            dis[i] = std::pair<uint32_t, float>(
                i + node->first_child, this->m_tree[i]->bbox.distanceTo(ray.o));
        }
        std::sort(dis, dis + 8, [](const auto &x, const auto &y) {
            return x.second < y.second;
        });

        for (auto [childIDX, d] : dis) {
            foundIntersection |= traverse(childIDX, ray, its, shadowRay, f);
            if (shadowRay && foundIntersection) {
                return true;
            }
        }
    }

    return foundIntersection;
}

bool BVH::traverse(uint32_t n, Ray3f &ray, Intersection &its, bool shadowRay,
                   uint32_t &f) const {}

NORI_NAMESPACE_END
