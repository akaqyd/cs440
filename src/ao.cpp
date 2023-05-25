
#include <nori/integrator.h>
#include <nori/scene.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

class AOIntegrator : public Integrator {
   public:
    AOIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Vector3f vWorld = its.shFrame.toWorld(Warp::squareToCosineHemisphere(sampler->next2D()));
        if (scene->rayIntersect(Ray3f(its.p, vWorld)))
            return Color3f(0.0f);
        else
            return Color3f(1.0f);

    }

    std::string toString() const {
        return tfm::format(
            "AOIntegrator[]\n");
    }
};

NORI_REGISTER_CLASS(AOIntegrator, "ao");
NORI_NAMESPACE_END