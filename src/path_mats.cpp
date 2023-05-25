

#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <algorithm>

#define ROULETTE_START_DEPTH 4

NORI_NAMESPACE_BEGIN


    class PathMatsIntegrator : public Integrator {
    private:
    public:
        PathMatsIntegrator(const PropertyList &props) {
        }

        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
            Color3f li(0.f);
            Color3f throughput = 1.f;
            float eta = 1.000277f, prev_eta = 1.000277f, continue_p = 1.f;   // start with the eta of air
            int depth = 1;
            Point3f prev_p = ray.o;
            BSDFQueryRecord bRec(Vector3f(0.f));
            Intersection its;

            while (true) {
                // Find the intersection point
                if (depth == 1) {
                    if (!scene->rayIntersect(ray, its)) {
                        for (Emitter *em: scene->getEmitters())
                            if (em->isInfinite())
                                li += em->le(ray);
                        break;
                    }
                } else {
                    Ray3f _ray(its.p, its.toWorld(bRec.wo));
                    if (!scene->rayIntersect(_ray, its)) {
                        for (Emitter *em: scene->getEmitters())
                            if (em->isInfinite())
                                li += em->le(_ray);
                        break;
                    }
                }

                // Aggregate the contribution from current intersection point if it is an emitter
                if (its.mesh->isEmitter()) {
                    li += throughput * its.mesh->getEmitter()->le(its.p, its.shFrame.n, (prev_p - its.p).normalized());
                }

                // Update stop probability
                continue_p = throughput.maxCoeff() * eta * eta;
                continue_p = continue_p < 0.99f ? continue_p : 0.99f;

                // Using Russian Roulette after ROULETTE_START_DEPTH
                if (depth >= ROULETTE_START_DEPTH) {
                    if (sampler->next1D() > continue_p)
                        break;
                    else
                        throughput /= continue_p;
                }

                // Sample the BSDF to get a new direction
                bRec.wi = its.shFrame.toLocal((prev_p - its.p).normalized());
                Color3f bsdfVal = its.mesh->getBSDF()->sample(bRec, sampler->next2D());
                if (bsdfVal.maxCoeff() == 0.f) break;  // stop if sampling failed

                // Update eta and throughput
                prev_eta /= bRec.eta;
                eta *= prev_eta;
                throughput *= bsdfVal;

                prev_p = its.p;
                ++depth;
            }

            return li;
        }

        std::string toString() const {
            return tfm::format(
                    "PathMatsIntegrator[]");
        }
    };

    NORI_REGISTER_CLASS(PathMatsIntegrator, "path_mats");
NORI_NAMESPACE_END
