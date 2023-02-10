/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#pragma once

#include <nori/object.h>

NORI_NAMESPACE_BEGIN

struct EmitterQueryRecord {
    Point3f ref;
    Point3f p;
    Normal3f n;
    Vector3f wi;
    float pdf;

    EmitterQueryRecord(Point3f ref) : ref(ref) {}
    EmitterQueryRecord(Point3f ref, Point3f p, Normal3f n)
        : ref(ref), p(p), n(n) {}
};

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
  public:
    virtual Color3f sample(const Mesh *mesh, Sampler *sampler,
                           EmitterQueryRecord &eRec) const = 0;

    virtual Color3f eval(const EmitterQueryRecord &eRec) const = 0;

    virtual float pdf(const Mesh *mesh,
                      const EmitterQueryRecord &eRec) const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.)
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }
};

NORI_NAMESPACE_END
