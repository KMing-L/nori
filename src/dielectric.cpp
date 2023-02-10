/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob
*/

#include <nori/bsdf.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
  public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        Vector3f normalLocal{0., 0., 1.};
        float cosThetaI = Frame::cosTheta(bRec.wi);

        float etaI, etaT;
        if (cosThetaI > 0.f) {
            etaI = m_extIOR;
            etaT = m_intIOR;
        } else {
            etaI = m_intIOR;
            etaT = m_extIOR;
            normalLocal = -normalLocal;
            cosThetaI = -cosThetaI;
        }
        bRec.eta = etaI / etaT;

        float sin2ThetaI = Frame::sinTheta2(bRec.wi);
        float sin2ThetaT = bRec.eta * bRec.eta * sin2ThetaI;

        if (sin2ThetaT >= 1.) {
            bRec.wo = Vector3f{-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z()};
        } else {
            float fresnelTerm = nori::fresnel(cosThetaI, etaI, etaT);

            if (sample.x() < fresnelTerm) {
                bRec.wo = Vector3f{-bRec.wi.x(), -bRec.wi.y(), bRec.wi.z()};
            } else {
                float cosThetaT = std::sqrt(1. - sin2ThetaT);

                bRec.wo = bRec.eta * -bRec.wi +
                          (bRec.eta * cosThetaI - cosThetaT) * normalLocal;
            }
        }

        bRec.measure = EDiscrete;

        return Color3f(1.);
    }

    std::string toString() const {
        return tfm::format("Dielectric[\n"
                           "  intIOR = %f,\n"
                           "  extIOR = %f\n"
                           "]",
                           m_intIOR, m_extIOR);
    }

  private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
