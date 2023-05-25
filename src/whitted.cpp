
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegrator : public Integrator {
   private:
   public:
    WhittedIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        // Find the intersection point and consider the infinite light source
        Intersection its;
        if (!scene->rayIntersect(ray, its)) {
            Color3f L(0.f);
            for (Emitter* em : scene->getEmitters())
                if (em->isInfinite())
                    L += em->le(ray);
            return L;
        }

        // check if the intersection point is specular
        if (!its.mesh->getBSDF()->isDiffuse()) {
            BSDFQueryRecord bRec(its.toLocal(-1 * ray.d).normalized());
            if (sampler->next1D() < 0.95) {
                Color3f bsdf = its.mesh->getBSDF()->sample(bRec, sampler);
                if (bsdf.maxCoeff() <= 0.f)
                    return 0.f;
                else
                    return bsdf * Li(scene, sampler, Ray3f(its.p, its.toWorld(bRec.wo))) / 0.95;
            } else
                return Color3f(0.0f);
        }

        // Sample an emitter from the scene
        const Emitter *em = scene->sampleEmitter(sampler->next1D());


        // Sample an incident direction and l_i from the light source
        Vector3f wi;
        Vector3f wo = (-ray.d).normalized();
        Point3f p;
        Normal3f n;
        Color3f l_i = em->sampleLi(its.p, its.shFrame.n, wi, p, n, sampler, scene);
        l_i *= scene->getEmittersSize();

        // BSDF
        BSDFQueryRecord bRec(its.toLocal(wi), its.toLocal(wo), ESolidAngle);
        Color3f bsdfVal = its.mesh->getBSDF()->eval(bRec);

        // Compute l_e if intersection is on a emitter
        Color3f l_e(0.0f);
        if (its.mesh->isEmitter())
            l_e = its.mesh->getEmitter()->le(its.p, its.shFrame.n, wo);

        Color3f val = l_e + bsdfVal * l_i;
        return val;
    }

    std::string toString() const {
        return tfm::format(
            "WhittedIntegrator[]");
    }
};

NORI_REGISTER_CLASS(WhittedIntegrator, "whitted");
NORI_NAMESPACE_END
