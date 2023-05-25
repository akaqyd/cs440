

#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <algorithm>

#define ROULETTE_START_DEPTH 4
#define MAX_DEPTH 16
NORI_NAMESPACE_BEGIN

    class PathVolMatsIntegrator : public Integrator {
    private:

    public:
        PathVolMatsIntegrator(const PropertyList &props) {
        }

        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const override {
            Color3f li(0.f);
            Color3f throughput = 1.f;
            float eta = 1.000277f, prev_eta = 1.000277f, continue_p = 1.f;
            int depth = 1;
//            bool specular_bounce = true;
            Ray3f ray(_ray);
            Point3f prev_p = ray.o;
            BSDFQueryRecord bRec(Vector3f(0.f));
            Intersection its;
//            MediumInteraction mi;

            while (true) {
                bool intersectSurface = scene->rayIntersect(ray, its);
                if (intersectSurface)
                    ray.maxt = (its.p - ray.o).norm();

                /* Ray has a corresponding medium. Sample the medium. */
                MediumInteraction mi;
                if (ray.m != nullptr)
                    throughput *= ray.m->sample(ray, sampler, &mi);

                if (throughput.maxCoeff() < Epsilon) break;

                /* Medium sampling has two cases */
                /* If it does sample an interaction at t < ray.maxt,
                 * the integral term along the interval is to be estimated,
                 * and the provided MediumInteraction should be initialized accordingly.
                 */
                if (mi.isInitialized()) {
                    if (depth > MAX_DEPTH) break;
//                    li += throughput *
//                          sampleDirectLight(scene, sampler, mi, -ray.d);   // account for direct illumination
                    Vector3f wo;
                    mi.mif.inside->sample_pf(-ray.d, wo, sampler);
                    mi.spawnRayToDir(ray, wo);
//                    specular_bounce = false;
                }
                    /* If Sample() doesn’t sample an interaction on the given ray interval [0, t_max],
                     * then the surface-related term Tr(p0→ p)Lo(p0, −ω) should be estimated.
                     * Falls back to original path_ems
                     */
                else {
                    if (intersectSurface) {
                        // Calculate emitted radiance L_e at its
//                        if (its.mesh->isEmitter() specular_bounce) {
                        if (its.mesh->isEmitter()) {
                            if (its.has_emitter_texture)
                                li += throughput *
                                      its.mesh->getEmitter()->le(its.p, its.shFrame.n, -ray.d, its.emitter_texture);
                            else
                                li += throughput * its.mesh->getEmitter()->le(its.p, its.shFrame.n, -ray.d);
                        }
                    } else {
//                        if (specular_bounce) {
                            // Calculate radiance from IBL
                            for (Emitter *em: scene->getEmitters())
                                if (em->isInfinite())
                                    li += throughput * em->le(ray);
                        break;
                    }

                    if (depth > MAX_DEPTH) break;

                    // If intersection is on a non-opaque surface that is used as the medium boundary,
                    // just skip over the boundary and continue
                    if (its.mesh->getBSDF()->isNull()) {
                        its.spawnRayToDir(ray, ray.d);
                        depth--;
                        continue;
                    }

                    // Account for direct illumination
//                    li += throughput * sampleDirectLight(scene, sampler, its, -ray.d);

                    // Sample the BSDF to get a new direction
                    bRec.wi = its.shFrame.toLocal(-ray.d);

                    // texture
                    its.setTexture(bRec);

                    if (its.mesh->getBSDF()->isDiffuse())
                        bRec.measure = ESolidAngle;
                    else
                        bRec.measure = EDiscrete;
                    Color3f bsdfVal = its.mesh->getBSDF()->sample(bRec, sampler);
                    if (bsdfVal.maxCoeff() == 0.f) break;  // stop if sampling failed
                    throughput *= bsdfVal;

                    prev_eta /= bRec.eta;
                    eta *= prev_eta;

//                    specular_bounce = !its.mesh->getBSDF()->isDiffuse();

                    its.spawnRayToDir(ray, its.toWorld(bRec.wo));
                }

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

                depth++;
            }

            return li;

        }


        [[nodiscard]] std::string toString() const override {
            return tfm::format(
                    "PathVolMatsIntegrator[]");
        }
    };

    NORI_REGISTER_CLASS(PathVolMatsIntegrator, "path_vol_mats");
NORI_NAMESPACE_END
