
#include <nori/integrator.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class SimpleIntegrator : public Integrator {
   private:
    Point3f position;
    Color3f energy;

   public:
    SimpleIntegrator(const PropertyList &props) {
        position = props.getPoint("position");
        energy   = props.getColor("energy");
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        /* Find the surface that is visible in the requested direction */
        Intersection its;
        if (!scene->rayIntersect(ray, its))
            return Color3f(0.0f);

        Vector3f xp = position - its.p;
        Ray3f shadowRay(its.p, xp);
        if (scene->rayIntersect(shadowRay))
            return Color3f(0.0f);

        /* Return the component-wise absolute
           value of the shading normal as a color */
        Normal3f n     = its.shFrame.n;
        float xpNorm = xp.norm();
        float cosTheta = n.dot(xp / xpNorm);
        if (cosTheta < 0.0f)
            return Color3f(0.0f);

        return energy * INV_PI * INV_PI / 4.0 * cosTheta / (xpNorm * xpNorm);
    }

    std::string toString() const {
        return tfm::format(
            "SimpleIntegrator[\n"
            "  position = %s,\n"
            "  energy   = %s,\n"
            "]",
            position.toString(), energy.toString());
    }
};

NORI_REGISTER_CLASS(SimpleIntegrator, "simple");
NORI_NAMESPACE_END