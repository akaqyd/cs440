
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

NORI_NAMESPACE_BEGIN

class WhittedIntegratorHacker : public Integrator {
   private:
   public:
    WhittedIntegratorHacker(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {

        // Can only be used in the hack task of HW 4, assume there is just one emitter
        // Immediately return if hit an emitter
        Emitter* _em = scene->rayIntersectNonMeshEmitter(ray);
        if (_em != nullptr) 
            return _em->le(Point3f(0.f), Normal3f(0.f), Vector3f(0.f));

        Intersection its;
        if (!scene->rayIntersect(ray, its)) 
                return Color3f(0.0f);

        // check if the intersection point is specular
        if (!its.mesh->getBSDF()->isDiffuse()) {
            BSDFQueryRecord bRec(its.toLocal(-1 * ray.d));
            if (sampler->next1D() < 0.95) {
                Color3f brdf = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                return brdf * Li(scene, sampler, Ray3f(its.p, its.toWorld(bRec.wo))) / 0.95;
            } else
                return Color3f(0.0f);
        }

        // Sample an emitter from the scene
        const Emitter *em = scene->sampleEmitter(sampler->next1D());
        
        // Sample an incident direction and l_i from the light source
        Vector3f wi;
        Vector3f wo = (-ray.d).normalized();
        Color3f l_i = em->sampleLi(its.p, its.shFrame.n, wi, sampler, scene);

        // BSDF
        BSDFQueryRecord bRec(its.toLocal(wi), its.toLocal(wo), ESolidAngle);
        Color3f bsdfVal = its.mesh->getBSDF()->eval(bRec);

        return bsdfVal * l_i;
    }

    std::string toString() const {
        return tfm::format(
            "WhittedIntegratorHacker[]");
    }
};


NORI_REGISTER_CLASS(WhittedIntegratorHacker, "whitted-hacker");
NORI_NAMESPACE_END
