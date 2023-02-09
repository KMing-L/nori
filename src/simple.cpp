#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
  public:
    SimpleIntegrator(const PropertyList &props) {
        p = props.getPoint("position");
        Phi = props.getColor("energy");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Vector3f l = p - its.p;

        if (scene->rayIntersect(Ray3f(its.p + l * Epsilon, l)))
            return Color3f(0.0f);

        return 0.25f * pow(INV_PI, 2) * Phi *
               std::max(0.0f, its.shFrame.n.dot(l.normalized())) / l.dot(l);
    }

    std::string toString() const { return "SimpleIntegrator[]"; }

  private:
    Point3f p;
    Color3f Phi;
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END