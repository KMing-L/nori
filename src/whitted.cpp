#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
  public:
    WhittedIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const {
        /* Find the surface that is visible in the requested direction */
        Ray3f ray = _ray;
        Color3f radiance(0.f);
        float coff = 1.f;

        while (true) {
            Intersection its;
            if (!scene->rayIntersect(ray, its))
                break;

            Color3f Le(0.f);
            if (its.mesh->isEmitter()) {
                EmitterQueryRecord eRec(ray.o, its.p, its.shFrame.n);
                Le = its.mesh->getEmitter()->eval(eRec);
            }
            radiance += Le;

            if (its.mesh->getBSDF()->isDiffuse()) {
                Color3f Lr(0.f);

                if (scene->getEmitterCount()) {
                    const Mesh *light =
                        scene->getRandomEmitter(sampler->next1D());

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
                                                 its.toLocal(eRec.wi),
                                                 ESolidAngle);
                            Color3f f = its.mesh->getBSDF()->eval(bRec);
                            Lr *= f * cosTheta * scene->getEmitterCount();
                        }
                    }
                }

                radiance += Lr;
                break;
            } else {
                BSDFQueryRecord bRec(its.toLocal(-ray.d));
                Color3f refColor =
                    its.mesh->getBSDF()->sample(bRec, sampler->next2D());

                if (refColor.isApprox(Color3f(0.f)) || sampler->next1D() > rr)
                    break;

                ray = Ray3f(its.p, its.toWorld(bRec.wo));

                coff *= rr;
            }
        }

        return radiance / coff;
    }

    std::string toString() const { return "WhittedIntegrator[]"; }

  private:
    float rr = 0.95;
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END