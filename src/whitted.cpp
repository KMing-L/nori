#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
  public:
    WhittedIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Color3f Le(0.f), Lr(0.f);

        if (its.mesh->isEmitter()) {
            EmitterQueryRecord eRec(ray.o, its.p, its.shFrame.n);
            Le = its.mesh->getEmitter()->eval(eRec);
        }

        if (scene->getEmitterDpdf()) {
            const Mesh *light = scene->getEmitterMesh(
                scene->getEmitterDpdf()->sample(sampler->next1D()));
            float lightPdf = scene->getEmitterDpdf()->getNormalization() *
                        light->getDPDF()->getSum();

            EmitterQueryRecord eRec(its.p);
            Lr = light->getEmitter()->sample(light, sampler, eRec);

            if (scene->rayIntersect(
                    Ray3f(eRec.ref, eRec.wi, Epsilon,
                          (eRec.p - eRec.ref).norm() - Epsilon))) {
                Lr = 0.f;
            } else {
                float cosTheta = its.shFrame.n.dot(eRec.wi);
                if (cosTheta <= 0) {
                    Lr = 0.f;
                } else {
                    BSDFQueryRecord bRec(its.toLocal(-ray.d),
                                         its.toLocal(eRec.wi), ESolidAngle);
                    Color3f f = its.mesh->getBSDF()->eval(bRec);
                    Lr *= f * cosTheta / lightPdf;
                }
            }
        }

        return Le + Lr;
    }

    std::string toString() const { return "WhittedIntegrator[]"; }
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END