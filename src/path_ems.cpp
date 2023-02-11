#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PathEmsIntegrator : public Integrator {
  public:
    PathEmsIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const {
        /* Find the surface that is visible in the requested direction */
        Ray3f ray = _ray;
        Color3f Lo(0.f);
        Color3f coff(1.f);
        uint32_t depth = 0;

        bool is_last_emitter = false;

        while (true) {
            if (depth >= 3) {
                float probability = std::min(coff.maxCoeff(), 0.99f);
                if (sampler->next1D() > probability)
                    break;
                coff /= probability;
            }

            Intersection its;
            if (!scene->rayIntersect(ray, its))
                break;

            Color3f Le(0.f);
            if (its.mesh->isEmitter() && !is_last_emitter) {
                EmitterQueryRecord eRec(ray.o, its.p, its.shFrame.n);
                Le = its.mesh->getEmitter()->eval(eRec) * coff;
            }
            Lo += Le;
            is_last_emitter = false;

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

                            Lr *=
                                coff * f * cosTheta * scene->getEmitterCount();
                        }
                    }
                }

                Lo += Lr;
                is_last_emitter = true;
            }

            BSDFQueryRecord bRec(its.toLocal(-ray.d));
            coff *= its.mesh->getBSDF()->sample(bRec, sampler->next2D());

            ray = Ray3f(its.p, its.toWorld(bRec.wo));

            depth++;
        }

        return Lo;
    }

    std::string toString() const { return "PathEmsIntegrator[]"; }
};

NORI_REGISTER_CLASS(PathEmsIntegrator, "path_ems");
NORI_NAMESPACE_END