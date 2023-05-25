

#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <algorithm>

#define ROULETTE_START_DEPTH 4
NORI_NAMESPACE_BEGIN

    class PathMisIntegrator : public Integrator {
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
                Color3f l_i = em->sampleLi(its.p, its.shFrame.n, wi, p, n, sampler, scene);
                l_i *= scene->getEmittersSize();

                if (l_i.maxCoeff() <= 0.f) return Color3f(0.f);

                BSDFQueryRecord bRec(its.toLocal(wo), its.toLocal(wi), ESolidAngle);
                its.setTexture(bRec);

                Color3f bsdfVal = its.mesh->getBSDF()->eval(bRec);

                if (bsdfVal.maxCoeff() <= 0.f)
                    return Color3f(0.f);

                // compute MIS weight
                float pdf_bsdf = its.mesh->getBSDF()->pdf(bRec);         // wrt solid angle
                float pdf_light = em->pdfLi(its.p, its.shFrame.n, p, n);  // wrt solid angle
                pdf_light /= scene->getEmittersSize();                    // the probability to sample an emitter

                if (pdf_light + pdf_bsdf <= 0.f)
                    return 0.f;
                else
                    return bsdfVal * l_i * pdf_light / (pdf_light + pdf_bsdf);
            }
        }

    public:
        PathMisIntegrator(const PropertyList &props) {
        }

        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &ray) const {
            Color3f li(0.f);
            Color3f throughput = 1.f;
            Color3f bsdfVal(0.f);
            float eta = 1.000277f, prev_eta = 1.000277f, continue_p = 1.f;
            int depth = 1;
            bool specular_bounce = true;
            Point3f prev_p = ray.o;
            BSDFQueryRecord bRec(Vector3f(0.f));
            Intersection its, next_its;

            if (!scene->rayIntersect(ray, next_its)) {
                for (Emitter *em: scene->getEmitters())
                    if (em->isInfinite())
                        li += em->le(ray);
                return li;
            }

            while (true) {
                its = next_its;

                Vector3f wo = (prev_p - its.p).normalized();

                // Calculate emitted radiance L_e at its
                if (its.mesh->isEmitter() && specular_bounce) {
                    if (its.has_emitter_texture)
                        li += throughput *
                              its.mesh->getEmitter()->le(its.p, its.shFrame.n, wo, its.emitter_texture);
                    else
                        li += throughput * its.mesh->getEmitter()->le(its.p, its.shFrame.n, wo);
                }

                // Calculate the direct illumination at its via emitter MIS
                li += throughput * Ld(scene, sampler, its, wo);

                // Sample the BSDF to get a new direction
                bRec.wi = its.shFrame.toLocal(wo);
                if (its.mesh->getBSDF()->isDiffuse())
                    bRec.measure = ESolidAngle;
                else
                    bRec.measure = EDiscrete;

                its.setTexture(bRec);

                bsdfVal = its.mesh->getBSDF()->sample(bRec, sampler);
                if (bsdfVal.maxCoeff() == 0.f) break;  // stop if sampling failed

                throughput *= bsdfVal;

                // next intersection along the sampled direction
                Ray3f nextRay(its.p, its.toWorld(bRec.wo));
                if (!scene->rayIntersect(nextRay, next_its)) {
                    if (!its.mesh->getBSDF()->isDiffuse()) {
                        for (Emitter *em: scene->getEmitters())
                            if (em->isInfinite())
                                li += throughput * em->le(nextRay);
                    } else {
                        for (Emitter *em: scene->getEmitters()) {
                            if (em->isInfinite()) {
                                float pdf_bsdf, pdf_light;
                                pdf_bsdf = its.mesh->getBSDF()->pdf(
                                        bRec);                                                            // wrt. solid angle
                                pdf_light = em->pdfLi(its.p, its.shFrame.n, nextRay.d,
                                                      Normal3f(0.f));  // wrt solid angle
                                pdf_light /= scene->getEmittersSize();                                                                 // the probability to sample an emitter

                                if (pdf_light + pdf_bsdf <= 0.f)
                                    li += 0.f;
                                else {
                                    float weight = pdf_bsdf / (pdf_light + pdf_bsdf);
                                    li += throughput * weight *
                                          em->le(nextRay);
                                }

                            }
                        }
                    }
                    break;
                }

                // new direction hit an emitter, calculate the direct illumination at its via BSDF MIS
                if (next_its.mesh->isEmitter()) {
                    float pdf_bsdf, pdf_light;
                    pdf_bsdf = its.mesh->getBSDF()->pdf(
                            bRec);                                                            // wrt. solid angle
                    pdf_light = next_its.mesh->getEmitter()->pdfLi(its.p, its.shFrame.n, next_its.p,
                                                                   next_its.shFrame.n);  // wrt solid angle
                    pdf_light /= scene->getEmittersSize();                                                                 // the probability to sample an emitter

                    if (pdf_light + pdf_bsdf <= 0.f)
                        li += 0.f;
                    else {
                        if (next_its.has_emitter_texture)
                            li += throughput * pdf_bsdf / (pdf_light + pdf_bsdf) *
                                  next_its.mesh->getEmitter()->le(next_its.p, next_its.shFrame.n,
                                                                  (its.p - next_its.p).normalized(),
                                                                  next_its.emitter_texture);
//                        li += throughput *
//                              its.mesh->getEmitter()->le(its.p, its.shFrame.n, -ray.d, its.emitter_texture);
                        else
                            li += throughput * pdf_bsdf / (pdf_light + pdf_bsdf) *
                                  next_its.mesh->getEmitter()->le(next_its.p, next_its.shFrame.n,
                                                                  (its.p - next_its.p).normalized());
                    }
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

                // Update eta and throughput
                prev_eta /= bRec.eta;
                eta *= prev_eta;

                specular_bounce = !its.mesh->getBSDF()->isDiffuse();
                prev_p = its.p;
                ++depth;
            }

            return li;
        }

        std::string toString() const {
            return tfm::format(
                    "PathMisIntegrator[]");
        }
    };

    NORI_REGISTER_CLASS(PathMisIntegrator, "path_mis");
NORI_NAMESPACE_END
