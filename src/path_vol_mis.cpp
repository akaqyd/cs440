//
// Created by Qiyuan Dong on 21.05.22.
//


#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/scene.h>

#include <algorithm>

#define ROULETTE_START_DEPTH 4
#define MAX_DEPTH 16
NORI_NAMESPACE_BEGIN

    class PathVolMISIntegrator : public Integrator {
    private:
        static Color3f
        sampleDirectLight(const Scene *scene, Sampler *sampler, const Intersection &its, const Vector3f &wo) {
            Color3f l_d(0.f);

            // For specular surface
            if (!its.mesh->getBSDF()->isDiffuse()) {
                return Color3f(0.f);
            } else {  // for diffuse surface
                Vector3f wi;
                Point3f p;
                Normal3f n;
                Emitter *em = scene->sampleEmitter(sampler->next1D());
                Color3f l_i = em->sampleLi_no_occlusion(its.p, its.shFrame.n, wi, p, n, sampler, scene);
                l_i *= scene->getEmittersSize();

                float pdf_light = em->pdfLi(its.p, its.shFrame.n, p, n);  // wrt solid angle
                pdf_light /= scene->getEmittersSize();                    // the probability to sample an emitter

                if (l_i.maxCoeff() > 0 && pdf_light > 0) {
                    BSDFQueryRecord bRec(its.toLocal(wo), its.toLocal(wi), ESolidAngle);
                    // texture
                    its.setTexture(bRec);
                    Color3f bsdfVal = its.mesh->getBSDF()->eval(bRec);

                    if (bsdfVal.maxCoeff() > 0) {
                        // Handle medium between intersection point and light source
                        Ray3f ray;
                        its.spawnRayToPoint(ray, p);
                        ray.mint = Epsilon;
                        ray.maxt -= Epsilon;
                        l_i *= scene->Tr(ray, sampler);
                        if (l_i.maxCoeff() > 0) {

                            // MIS
                            float pdf_bsdf = its.mesh->getBSDF()->pdf(bRec);         // wrt solid angle
                            if (pdf_light + pdf_bsdf > 0.f) {
                                if (em->isDeltaLight())
                                    l_d += bsdfVal * l_i;
                                else
                                    l_d += bsdfVal * l_i * pdf_light / (pdf_light + pdf_bsdf);
                            }
                        }
                    }
                }
            }

            // MIS, BSDF
            BSDFQueryRecord bRec(its.shFrame.toLocal(wo));
            if (its.mesh->getBSDF()->isDiffuse())
                bRec.measure = ESolidAngle;
            else
                bRec.measure = EDiscrete;
            Color3f bsdfVal = its.mesh->getBSDF()->sample(bRec, sampler);
            float pdf_bsdf = its.mesh->getBSDF()->pdf(bRec);


            if (bsdfVal.maxCoeff() > 0 && pdf_bsdf > 0) {
                Intersection next_its;
                Ray3f nextRay(its.p, its.toWorld(bRec.wo));
                Color3f tr;
                if (scene->rayIntersectTr(nextRay, sampler, next_its, tr)) {
                    if (next_its.mesh->isEmitter() && !next_its.mesh->getEmitter()->isDeltaLight()) {
                        float pdf_light = next_its.mesh->getEmitter()->pdfLi(its.p, its.shFrame.n, next_its.p,
                                                                             next_its.shFrame.n);  // wrt solid angle
                        pdf_light /= scene->getEmittersSize();                                                                 // the probability to sample an emitter

                        if (pdf_light + pdf_bsdf <= 0.f)
                            l_d += 0.f;
                        else {
                            float weight = pdf_bsdf / (pdf_light + pdf_bsdf);
                            if (next_its.has_emitter_texture)
                                l_d += bsdfVal * weight * tr * next_its.mesh->getEmitter()->le(next_its.p, next_its.shFrame.n,
                                                                                     (its.p - next_its.p).normalized(),
                                                                                     next_its.emitter_texture);
                            else
                                l_d += bsdfVal * weight * tr * next_its.mesh->getEmitter()->le(next_its.p, next_its.shFrame.n,
                                                                                     (its.p - next_its.p).normalized());
                        }

                    }
                }
            }
            return l_d;
        }

        static Color3f
        sampleDirectLight(const Scene *scene, Sampler *sampler, const MediumInteraction &mi, const Vector3f &wo) {
            Color3f l_d(0.f);

            // Uniformly sample an emitter
            Emitter *em = scene->sampleEmitter(sampler->next1D());

            Vector3f wi;
            Point3f p;
            Normal3f n;
            Color3f l_i = em->sampleLi_no_occlusion(mi.p, Normal3f(0, 0, 0), wi, p, n, sampler, scene);
            l_i *= (float) scene->getEmittersSize();

            float pdf_light = em->pdfLi(mi.p, Normal3f(0.f), p, n);  // wrt solid angle
            pdf_light /= scene->getEmittersSize();                    // the probability to sample an emitter

            ///
            if (l_i.maxCoeff() > 0 && pdf_light > 0) {
                float pfVal = mi.mif.inside->eval_pf(wi, wo);
                float scattering_pdf = pfVal;
                if (pfVal > 0) {
                    // Handle medium between intersection point and light source
                    Ray3f ray;
                    mi.spawnRayToPoint(ray, p);
                    ray.mint = Epsilon;
                    ray.maxt -= Epsilon;
                    l_i *= scene->Tr(ray, sampler);
                    if (l_i.maxCoeff() > 0) {
                        // MIS
                        if (pdf_light + scattering_pdf > 0.f) {
                            if (em->isDeltaLight())
                                l_d += pfVal * l_i;
                            else
                                l_d += pfVal * l_i * pdf_light / (pdf_light + scattering_pdf);
                        }
                    }
                }
            }
            ///

            // MIS, phase func
            Vector3f scattering_dir;
            float pfVal = mi.mif.inside->sample_pf(wo, scattering_dir, sampler);
            float pdf_scattering = pfVal;


            if (pdf_scattering > 0) {
                Intersection next_its;
                Ray3f nextRay(mi.p, scattering_dir);
                Color3f tr;
                if (scene->rayIntersectTr(nextRay, sampler, next_its, tr)) {
                    if (next_its.mesh->isEmitter() && !next_its.mesh->getEmitter()->isDeltaLight()) {
                        pdf_light = next_its.mesh->getEmitter()->pdfLi(mi.p, Normal3f(0.f), next_its.p,
                                                                             next_its.shFrame.n);  // wrt solid angle
                        pdf_light /= scene->getEmittersSize();                                                                 // the probability to sample an emitter

                        if (pdf_light <= 0)
                            return l_d;

                        if (pdf_light + pdf_scattering <= 0.f)
                            return l_d;
                        else {
                            float weight = pdf_scattering / (pdf_light + pdf_scattering);
                            if (next_its.has_emitter_texture)
                                l_d += weight * tr * next_its.mesh->getEmitter()->le(next_its.p, next_its.shFrame.n,
                                                                                     (mi.p - next_its.p).normalized(),
                                                                                     next_its.emitter_texture);
                            else
                                l_d += weight * tr * next_its.mesh->getEmitter()->le(next_its.p, next_its.shFrame.n,
                                                                                     (mi.p - next_its.p).normalized());
                        }
                    }
                }
            }
            return l_d;
        }

    public:
        PathVolMISIntegrator(const PropertyList &props) {
        }

        Color3f Li(const Scene *scene, Sampler *sampler, const Ray3f &_ray) const override {
            Color3f li(0.f);
            Color3f throughput = 1.f;
            float eta = 1.000277f, prev_eta = 1.000277f, continue_p = 1.f;
            int depth = 1;
            bool specular_bounce = true;
            Ray3f ray(_ray);
            Point3f prev_p = ray.o;
            BSDFQueryRecord bRec(Vector3f(0.f));
            Intersection its;
//            MediumInteraction mi;

            while (true) {
                bool intersectSurface = scene->rayIntersect(ray, its);
                if (intersectSurface)
                    ray.maxt = (its.p - ray.o).norm() / ray.d.norm();

                /* Ray has a corresponding medium. Sample the medium. */
                MediumInteraction mi;
                if (ray.m != nullptr)
                    throughput *= ray.m->sample(ray, sampler, &mi);

                if (throughput.maxCoeff() < Epsilon)
                    break;

                /* Medium sampling has two cases */
                /* If it does sample an interaction at t < ray.maxt,
                 * the integral term along the interval is to be estimated,
                 * and the provided MediumInteraction should be initialized accordingly.
                 */
                if (mi.isInitialized()) {
                    if (depth > MAX_DEPTH)
                        break;
                    li += throughput *
                          sampleDirectLight(scene, sampler, mi, -ray.d);   // account for direct illumination
                    Vector3f wo;
                    mi.mif.inside->sample_pf(-ray.d, wo, sampler);
                    mi.spawnRayToDir(ray, wo);
                    specular_bounce = false;
                }
                    /* If Sample() doesn’t sample an interaction on the given ray interval [0, t_max],
                     * then the surface-related term Tr(p0→ p)Lo(p0, −ω) should be estimated.
                     * Falls back to original path_ems
                     */
                else {
                    if (intersectSurface) {
                        // Calculate emitted radiance L_e at its
                        if (its.mesh->isEmitter() && specular_bounce) {
                            if (its.has_emitter_texture)
                                li += throughput *
                                      its.mesh->getEmitter()->le(its.p, its.shFrame.n, -ray.d, its.emitter_texture);
                            else
                                li += throughput * its.mesh->getEmitter()->le(its.p, its.shFrame.n, -ray.d);
                        }
                    } else {
                        if (specular_bounce) {
                            // Calculate radiance from IBL
                            for (Emitter *em: scene->getEmitters())
                                if (em->isInfinite())
                                    li += throughput * em->le(ray);
                        }
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
                    li += throughput * sampleDirectLight(scene, sampler, its, -ray.d);

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

                    specular_bounce = !its.mesh->getBSDF()->isDiffuse();

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
                    "PathVolMISIntegrator[]");
        }
    };

    NORI_REGISTER_CLASS(PathVolMISIntegrator, "path_vol_mis");
NORI_NAMESPACE_END
