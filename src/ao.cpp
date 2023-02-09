#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AoIntegrator : public Integrator {
  public:
    AoIntegrator(const PropertyList &props) {}

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        auto dir = Warp::squareToCosineHemisphere(sampler->next2D());
        auto pdf = Warp::squareToCosineHemispherePdf(dir);
        auto cos = its.shFrame.cosTheta(dir);
        dir = its.shFrame.toWorld(dir);

        if (scene->rayIntersect(Ray3f(its.p + dir * Epsilon, dir)))
            return Color3f(0.0f);

        return Color3f(1.0f) * INV_PI * cos / pdf;
    }

    std::string toString() const { return "AoIntegrator[]"; }
};

NORI_REGISTER_CLASS(AoIntegrator, "ao");
NORI_NAMESPACE_END