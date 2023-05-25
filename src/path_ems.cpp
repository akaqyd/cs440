

#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <algorithm>

#define ROULETTE_START_DEPTH 4
NORI_NAMESPACE_BEGIN

class PathEmsIntegrator : public Integrator {
   private:
    Color3f Ld(const Scene *scene, Sampler *sampler, const Intersection &its, const Vector3f &wo) const {
        // For specular surface
        if (!its.mesh->getBSDF()->isDiffuse()) {
            return Color3f(0.f);
        } else {  // for diffuse surface
            Vector3f wi;
            Point3f p;
            Normal3f n;
            const Emitter *em = scene->sampleEmitter(sampler->next1D());
            Color3f l_i       = em->sampleLi(its.p, its.shFrame.n, wi, p, n, sampler, scene);
            l_i *= scene->getEmittersSize();

            if (l_i.maxCoeff() <= 0.f) return Color3f(0.f);

            BSDFQueryRecord bRec(its.toLocal(wo), its.toLocal(wi), ESolidAngle);
            Color3f bsdfVal = its.mesh->getBSDF()->eval(bRec);

            return bsdfVal * l_i;
        }
    }

   public:
    PathEmsIntegrator(const PropertyList &props) {
    }

    Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
        Color3f li(0.f);
        Color3f throughput = 1.f;
        float eta = 1.000277f, prev_eta = 1.000277f, continue_p = 1.f;
        int depth            = 1;
        bool specular_bounce = true;
        Point3f prev_p       = ray.o;
        BSDFQueryRecord bRec(Vector3f(0.f));
        Intersection its;

        while (true) {
            // Find the intersection point
            if (depth == 1) {
                if (!scene->rayIntersect(ray, its)) {
                    for (Emitter *em : scene->getEmitters())
                        if (em->isInfinite())
                            li += em->le(ray);
                    break;
                }
            } else {
                Ray3f _ray(its.p, its.toWorld(bRec.wo));
                if (!scene->rayIntersect(_ray, its)) {
                    if (specular_bounce) {
                        for (Emitter *em : scene->getEmitters())
                            if (em->isInfinite())
                                li += em->le(_ray);
                        break;
                    }
                }
            }

            // Calculate the L_e contribution of current intersection point if it is an emitter
            Vector3f wo = (prev_p - its.p).normalized();
            if (its.mesh->isEmitter() && specular_bounce)
                li += throughput * its.mesh->getEmitter()->le(its.p, its.shFrame.n, wo);

            // Calculate the L_d contribution of current intersection (direct illumination)
            li += throughput * Ld(scene, sampler, its, wo);

            // Update stop probability
            continue_p = throughput.maxCoeff() * eta * eta;
            continue_p = continue_p < 0.99f ? continue_p : 0.99f;

            // Russian Roulette
            if (depth >= ROULETTE_START_DEPTH) {
                if (sampler->next1D() > continue_p)
                    break;
                else
                    throughput /= continue_p;
            }

            // Sample the BSDF to get a new direction
            bRec.wi         = its.shFrame.toLocal(wo);
            Color3f bsdfVal = its.mesh->getBSDF()->sample(bRec, sampler);
            if (bsdfVal.maxCoeff() == 0.f) break;  // stop if sampling failed

            // Update eta and throughput
            prev_eta /= bRec.eta;
            eta *= prev_eta;
            throughput *= bsdfVal;

            specular_bounce = !its.mesh->getBSDF()->isDiffuse();
            prev_p          = its.p;
            ++depth;
        }

        return li;
    }

    std::string toString() const {
        return tfm::format(
            "PathEmsIntegrator[]");
    }
};

NORI_REGISTER_CLASS(PathEmsIntegrator, "path_ems");
NORI_NAMESPACE_END
