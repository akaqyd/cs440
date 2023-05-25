/*
    This file is part of Nori, a simple educational ray tracer

    Copyright (c) 2015 by Wenzel Jakob

    Nori is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Nori is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

    class Microfacet : public BSDF {
    public:
        Microfacet(const PropertyList &propList) {
            /* RMS surface roughness */
            m_alpha = propList.getFloat("alpha", 0.1f);

            /* Interior IOR (default: BK7 borosilicate optical glass) */
            m_intIOR = propList.getFloat("intIOR", 1.5046f);

            /* Exterior IOR (default: air) */
            m_extIOR = propList.getFloat("extIOR", 1.000277f);

            /* Albedo of the diffuse base material (a.k.a "kd") */
            m_kd = propList.getColor("kd", Color3f(0.5f));

            /* To ensure energy conservation, we must scale the
               specular component by 1-kd.

               While that is not a particularly realistic model of what
               happens in reality, this will greatly simplify the
               implementation. Please see the course staff if you're
               interested in implementing a more realistic version
               of this BRDF. */
            m_ks = 1 - m_kd.maxCoeff();
        }

        /// Evaluate the BRDF for the given pair of directions
        Color3f eval(const BSDFQueryRecord &bRec) const {
            // need to ensure that wi and wo in the bRec are normalized
            Color3f kd = m_kd;
            float ks = m_ks;
            float alpha = m_alpha;

            if (bRec.has_color_texture) {
                kd = bRec.color_texture;
                ks = 1 - kd.maxCoeff();
            }
            if (bRec.has_roughness_texture)
                alpha = bRec.roughness_texture;


            Color3f f = kd * INV_PI;

            Vector3f wh = (bRec.wi + bRec.wo).normalized();

            // check visibility first
            if (abs(bRec.wi.z()) < Epsilon || abs(bRec.wo.z()) < Epsilon ||
                bRec.wi.dot(wh) / bRec.wi.z() <= 0 ||
                bRec.wo.dot(wh) / bRec.wo.z() <= 0)
                return f;

            // Fresnel
            float Fr = fresnel(bRec.wi.dot(wh), m_extIOR, m_intIOR);

            // pdf of wh
//        float D_wh = Warp::squareToBeckmannPdf(wh, alpha);
            float D_wh = Warp::squareToGGXPdf(wh, alpha);

            f += ks * D_wh * Fr * G1(bRec.wi, wh, alpha) * G1(bRec.wo, wh, alpha) /
                 (4 * bRec.wi.z() * bRec.wo.z() * wh.z());

            return f;
        }

        /// Evaluate the sampling density of \ref sample() wrt. solid angles
        float pdf(const BSDFQueryRecord &bRec) const {
            Color3f kd = m_kd;
            float ks = m_ks;
            float alpha = m_alpha;

            if (bRec.has_color_texture) {
                kd = bRec.color_texture;
                ks = 1 - kd.maxCoeff();
            }
            if (bRec.has_roughness_texture)
                alpha = bRec.roughness_texture;

            if (bRec.measure != ESolidAngle || Frame::cosTheta(bRec.wi) <= 0 || Frame::cosTheta(bRec.wo) <= 0)
                return 0.0f;

            Vector3f wh = (bRec.wi + bRec.wo).normalized();
//        return m_ks * Warp::squareToBeckmannPdf(wh, alpha) / (4 * wh.dot(bRec.wo)) +
//               (1 - m_ks) * INV_PI * bRec.wo.z();
            return ks * Warp::squareToGGXPdf(wh, alpha) / (4 * wh.dot(bRec.wo)) +
                   (1 - ks) * INV_PI * bRec.wo.z();
        }

        /// Sample the BRDF
        Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
            if (bRec.wi.z() <= 0)
                return Color3f(0, 0, 0);

            Color3f kd = m_kd;
            float ks = m_ks;
            float alpha = m_alpha;

            if (bRec.has_color_texture) {
                kd = bRec.color_texture;
                ks = 1 - kd.maxCoeff();
            }
            if (bRec.has_roughness_texture)
                alpha = bRec.roughness_texture;

            // Sample reuse
            Point2f sample(_sample);

            if (_sample[0] < ks) {  // Specular
                sample[0] = sample[0] / ks;
//            Vector3f wh = Warp::squareToBeckmann(sample, alpha);
                Vector3f wh = Warp::squareToGGX(sample, alpha);
                bRec.wo = (2.0f * wh.dot(bRec.wi) * wh - bRec.wi).normalized();

            } else {  // Diffuse
                sample[0] = (sample[0] - ks) / (1.0f - ks);
                bRec.wo = Warp::squareToCosineHemisphere(sample);
            }

            bRec.measure = ESolidAngle;
            // bRec.eta     = 1.0f;

            if (bRec.wo.z() < 0.0f)
                return Color3f(0.0);

            return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
            // Note: Once you have implemented the part that computes the scattered
            // direction, the last part of this function should simply return the
            // BRDF value divided by the solid angle density and multiplied by the
            // cosine factor from the reflection equation, i.e.
            // return eval(bRec) * Frame::cosTheta(bRec.wo) / pdf(bRec);
        }

        bool isDiffuse() const {
            /* While microfacet BRDFs are not perfectly diffuse, they can be
               handled by sampling techniques for diffuse/non-specular materials,
               hence we return true here */
            return true;
        }

        std::string toString() const {
            return tfm::format(
                    "Microfacet[\n"
                    "  alpha = %f,\n"
                    "  intIOR = %f,\n"
                    "  extIOR = %f,\n"
                    "  kd = %s,\n"
                    "  ks = %f\n"
                    "]",
                    m_alpha,
                    m_intIOR,
                    m_extIOR,
                    m_kd.toString(),
                    m_ks);
        }

    private:
        float m_alpha;
        float m_intIOR, m_extIOR;
        float m_ks;
        Color3f m_kd;

        float G1(const Vector3f &wv, const Vector3f &wh, const float &alpha) const {
            float b = 1.0f / (alpha * Frame::tanTheta(wv));
            if (b < 1.6)
                return (3.535 * b + 2.181 * b * b) / (1 + 2.276 * b + 2.577 * b * b);
            else
                return 1.0f;
        }
    };

    NORI_REGISTER_CLASS(Microfacet, "microfacet");
NORI_NAMESPACE_END


