/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <Eigen/Geometry>
#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/dpdf.h>
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

Mesh::Mesh() {}

Mesh::~Mesh() {
    delete m_bsdf;
    delete m_emitter;
}

void Mesh::activate() {
    if (!m_bsdf) {
        /* If no material was assigned, instantiate a diffuse BRDF */
        m_bsdf = static_cast<BSDF *>(
            NoriObjectFactory::createInstance("diffuse", PropertyList()));
    }

    m_dpdf = new DiscretePDF(getTriangleCount());
    for (uint32_t i = 0, n = getTriangleCount(); i < n; i++)
        m_dpdf->append(surfaceArea(i));
    m_dpdf->normalize();
}

float Mesh::surfaceArea(uint32_t index) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
}

bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v,
                        float &t) const {
    uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
    const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

    /* Find vectors for two edges sharing v[0] */
    Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

    /* Begin calculating determinant - also used to calculate U parameter */
    Vector3f pvec = ray.d.cross(edge2);

    /* If determinant is near zero, ray lies in plane of triangle */
    float det = edge1.dot(pvec);

    if (det > -1e-8f && det < 1e-8f)
        return false;
    float inv_det = 1.0f / det;

    /* Calculate distance from v[0] to ray origin */
    Vector3f tvec = ray.o - p0;

    /* Calculate U parameter and test bounds */
    u = tvec.dot(pvec) * inv_det;
    if (u < 0.0 || u > 1.0)
        return false;

    /* Prepare to test V parameter */
    Vector3f qvec = tvec.cross(edge1);

    /* Calculate V parameter and test bounds */
    v = ray.d.dot(qvec) * inv_det;
    if (v < 0.0 || u + v > 1.0)
        return false;

    /* Ray intersects triangle -> compute t */
    t = edge2.dot(qvec) * inv_det;

    return t >= ray.mint && t <= ray.maxt;
}

BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
    BoundingBox3f result(m_V.col(m_F(0, index)));
    result.expandBy(m_V.col(m_F(1, index)));
    result.expandBy(m_V.col(m_F(2, index)));
    return result;
}

Point3f Mesh::getCentroid(uint32_t index) const {
    return (1.0f / 3.0f) * (m_V.col(m_F(0, index)) + m_V.col(m_F(1, index)) +
                            m_V.col(m_F(2, index)));
}

void Mesh::addChild(NoriObject *obj) {
    switch (obj->getClassType()) {
    case EBSDF:
        if (m_bsdf)
            throw NoriException(
                "Mesh: tried to register multiple BSDF instances!");
        m_bsdf = static_cast<BSDF *>(obj);
        break;

    case EEmitter: {
        Emitter *emitter = static_cast<Emitter *>(obj);
        if (m_emitter)
            throw NoriException(
                "Mesh: tried to register multiple Emitter instances!");
        m_emitter = emitter;
    } break;

    default:
        throw NoriException("Mesh::addChild(<%s>) is not supported!",
                            classTypeName(obj->getClassType()));
    }
}

MeshSample Mesh::sample(Sampler *sampler) const {
    MeshSample sample;
    uint32_t idx = m_dpdf->sample(sampler->next1D());

    Point2f pointSample = sampler->next2D();
    float alpha = 1 - sqrt(1 - pointSample.x());
    float beta = pointSample.y() * (1 - pointSample.x());
    float gamma = 1 - alpha - beta;

    Point3f v0 = m_V.col(m_F(0, idx));
    Point3f v1 = m_V.col(m_F(1, idx));
    Point3f v2 = m_V.col(m_F(2, idx));

    sample.p = alpha * v0 + beta * v1 + gamma * v2;
    if (m_N.size() != 0) {
        sample.n = alpha * m_N.col(m_F(0, idx)) + beta * m_N.col(m_F(1, idx)) +
                   gamma * m_N.col(m_F(2, idx));
    } else {
        sample.n = (v1 - v0).cross(v2 - v1).normalized();
    }

    sample.pdf = m_dpdf->getNormalization();

    return sample;
}

std::string Mesh::toString() const {
    return tfm::format(
        "Mesh[\n"
        "  name = \"%s\",\n"
        "  vertexCount = %i,\n"
        "  triangleCount = %i,\n"
        "  bsdf = %s,\n"
        "  emitter = %s\n"
        "]",
        m_name, m_V.cols(), m_F.cols(),
        m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
        m_emitter ? indent(m_emitter->toString()) : std::string("null"));
}

std::string Intersection::toString() const {
    if (!mesh)
        return "Intersection[invalid]";

    return tfm::format("Intersection[\n"
                       "  p = %s,\n"
                       "  t = %f,\n"
                       "  uv = %s,\n"
                       "  shFrame = %s,\n"
                       "  geoFrame = %s,\n"
                       "  mesh = %s\n"
                       "]",
                       p.toString(), t, uv.toString(),
                       indent(shFrame.toString()), indent(geoFrame.toString()),
                       mesh ? mesh->toString() : std::string("null"));
}

NORI_NAMESPACE_END
