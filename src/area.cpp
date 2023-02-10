#include <nori/emitter.h>
#include <nori/frame.h>
#include <nori/mesh.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AreaLight : public Emitter {
  public:
    AreaLight(const PropertyList &propList) {
        m_radiance = propList.getColor("radiance");
    }

    Color3f eval(const EmitterQueryRecord &eRec) const override {
        if (eRec.n.dot(-eRec.wi) <= 0.f)
            return Color3f(0.f);
        return m_radiance;
    }

    float pdf(const Mesh *mesh, const EmitterQueryRecord &eRec) const override {
        float cosTheta = eRec.n.dot(-eRec.wi);
        if (cosTheta <= 0.f)
            return 0.f;

        return mesh->getDPDF()->getNormalization() *
               (eRec.p - eRec.ref).squaredNorm() / cosTheta;
    }

    Color3f sample(const Mesh *mesh, Sampler *sampler,
                   EmitterQueryRecord &eRec) const override {
        auto result = mesh->sample(sampler);
        eRec.p = result.p;
        eRec.n = result.n;
        eRec.wi = (eRec.p - eRec.ref).normalized();
        eRec.pdf = pdf(mesh, eRec);
        
        if (eRec.pdf == 0.f)
            return Color3f(0.f);
        
        return eval(eRec) / eRec.pdf;
    }

    std::string toString() const {
        return tfm::format("AreaLight[\n"
                           "  radiance = %s\n"
                           "]",
                           m_radiance.toString());
    }

    EClassType getClassType() const { return EEmitter; }

  private:
    Color3f m_radiance;
};

NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END