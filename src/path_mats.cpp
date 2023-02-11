#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class PathMatsIntegrator : public Integrator {
  public:
    PathMatsIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const {
        /* Find the surface that is visible in the requested direction */
        Ray3f ray = _ray;
        Color3f Lo(0.f);
        Color3f fr(1.f);
        uint32_t depth = 0;

        while (true) {
            if (depth >= 3) {
                float probability = std::min(fr.maxCoeff(), 0.99f);
                if (sampler->next1D() > probability)
                    break;
                fr /= probability;
            }

            Intersection its;
            if (!scene->rayIntersect(ray, its))
                break;

            Color3f Le(0.f);
            if (its.mesh->isEmitter()) {
                EmitterQueryRecord eRec(ray.o, its.p, its.shFrame.n);
                Le = its.mesh->getEmitter()->eval(eRec) * fr;
            }
            Lo += Le;

            BSDFQueryRecord bRec(its.toLocal(-ray.d));
            fr *= its.mesh->getBSDF()->sample(bRec, sampler->next2D());

            ray = Ray3f(its.p, its.toWorld(bRec.wo));

            depth++;
        }

        return Lo;
    }

    std::string toString() const { return "PathMatsIntegrator[]"; }
};

NORI_REGISTER_CLASS(PathMatsIntegrator, "path_mats");
NORI_NAMESPACE_END