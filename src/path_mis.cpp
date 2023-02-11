#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PathMisIntegrator : public Integrator {
  public:
    PathMisIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const {
        Color3f Lo(0.f);
        Color3f coff(1.f);
        Ray3f ray(_ray);

        float w_mats = 1.f, w_ems = 1.f;

        uint32_t depth = 0;

        Intersection its;

        if (!scene->rayIntersect(ray, its))
            return Lo;

        while (true) {
            if (depth >= 3) {
                float probability = std::min(coff.maxCoeff(), 0.99f);
                if (sampler->next1D() > probability)
                    break;
                coff /= probability;
            }

            Color3f Le(0.f);
            if (its.mesh->isEmitter()) {
                EmitterQueryRecord eRec(ray.o, its.p, its.shFrame.n);
                Le = w_mats * its.mesh->getEmitter()->eval(eRec) * coff;
            }
            Lo += Le;

            if (its.mesh->getBSDF()->isDiffuse()) {
                Color3f Lr(0.f);

                if (scene->getEmitterCount()) {
                    const Mesh *light =
                        scene->getRandomEmitter(sampler->next1D());

                    EmitterQueryRecord eRec(its.p);
                    Lr = light->getEmitter()->sample(light, sampler, eRec);
                    float pdf_em = eRec.pdf;

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

                            float pdf_mat = its.mesh->getBSDF()->pdf(bRec);

                            w_ems = pdf_mat + pdf_em > 0.f
                                        ? pdf_em / (pdf_mat + pdf_em)
                                        : pdf_em;

                            Lr *= w_ems * coff * f * cosTheta *
                                  scene->getEmitterCount();
                        }
                    }
                }

                Lo += Lr;
            }

            BSDFQueryRecord bRec(its.toLocal(-ray.d));
            coff *= its.mesh->getBSDF()->sample(bRec, sampler->next2D());
            ray = Ray3f(its.p, its.toWorld(bRec.wo));

            float pdf_mat = its.mesh->getBSDF()->pdf(bRec);
            Point3f origin = its.p;

            if (!scene->rayIntersect(ray, its))
                break;

            if (its.mesh->isEmitter()) {
                EmitterQueryRecord eRec(origin, its.p, its.shFrame.n);

                float pdf_em = its.mesh->getEmitter()->pdf(its.mesh, eRec);
                w_mats = pdf_mat + pdf_em > 0.f ? pdf_mat / (pdf_mat + pdf_em)
                                                : pdf_mat;
            }

            if (bRec.measure == EDiscrete)
                w_mats = 1.f;

            depth++;
        }

        return Lo;
    };

    std::string toString() const { return "PathMisIntegrator[]"; }
};

NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis");
NORI_NAMESPACE_END