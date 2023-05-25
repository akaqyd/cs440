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

#include <iomanip>
#include <string>

NORI_NAMESPACE_BEGIN

class MicrofacetHacker : public BSDF {
   public:
    MicrofacetHacker(const PropertyList &propList) {
        /* RMS surface roughness */
        m_alpha = propList.getFloat("alpha", 0.1f);

        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &bRec) const {
        // need to ensure that wi and wo in the bRec are normalized
        float brdfVal = brdf(bRec.wi, bRec.wo, Normal3f(0, 0, 1));
        float btdfVal = btdf(bRec.wi, bRec.wo, Normal3f(0, 0, 1));

        return brdfVal + btdfVal;
    }

    /// Evaluate the sampling density of \ref sample() wrt. solid angles
    float pdf(const BSDFQueryRecord &bRec) const {
        throw NoriException("Full microfacet BSDF need a sampler, instead of just one sample!");
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &_sample) const {
        throw NoriException("Full microfacet BSDF need a sampler, instead of just one sample!");
    }

    /// Sample the BRDF
    Color3f sample(BSDFQueryRecord &bRec, Sampler *sampler) const {
        if (fabs(bRec.wi.z()) < Epsilon)
            return Color3f(0.f);

        // macrofacet normal n
        Normal3f n(0, 0, 1);

        // sample a microfacet normal m
        Normal3f m = Warp::squareToGGX(sampler->next2D(), m_alpha);

        // m need to be on the same side as wi
        if (bRec.wi.z() * m.z() < 0)
            m *= -1;

        // eta = idx_incident_side / idx_transmitted_side
        float etaI = m_extIOR, etaT = m_intIOR;
        if (bRec.wi.dot(n) < 0)
            std::swap(etaI, etaT);
        float eta = etaI / etaT;

        float F           = _F(bRec.wi, m);
        std::string event = "Reflection";

        float pdf = 0.f;
        if (sampler->next1D() < F) {  // Specular Reflection
            bRec.wo = (2.0f * fabs(m.dot(bRec.wi)) * m - bRec.wi).normalized();
            if (bRec.wo.z() * bRec.wi.z() < 0)  // wo must be in the same hemisphere as wi
                return 0.f;
            bRec.eta = 1.f;

        } else {  // Specular Transmission
            if (!refract(bRec.wi, m, eta, bRec.wo))
                return 0.f;                     // invalid wo
            if (bRec.wo.z() * bRec.wi.z() > 0)  // wo must be in the different hemisphere as wi
                return 0.f;
            bRec.eta = eta;
            event    = "Refraction";

        }

        bRec.measure = ESolidAngle;

        float reflectP = reflectionPdf(bRec.wi, bRec.wo);
        float refractP = refractionPdf(bRec.wi, bRec.wo);

        pdf = F * reflectP + (1 - F) * refractP;

        float reflectVal = brdf(bRec.wi, bRec.wo, Normal3f(0, 0, 1));
        float refractVal = btdf(bRec.wi, bRec.wo, Normal3f(0, 0, 1));
        float bsdf       = refractVal + reflectVal;

        float val = bsdf * fabs(bRec.wo.dot(n)) / pdf * bRec.eta * bRec.eta;
        if (isnan(val) || isinf(val))
            return 0.f;
        else
            return val;

        // // if (bsdf < Epsilon || pdf < 1e-4) {
        // if (val > 100) {
        //     cout << "=============================================================" << endl;
        //     cout << "Event: " << event << endl;
        //     cout << "wi: " << endl
        //          << bRec.wi << endl;
        //     cout << "wo: " << endl
        //          << bRec.wo << endl;
        //     cout << "m: " << endl
        //          << m << endl;
        //     cout << "Fresnel: " << F << endl;
        //     cout << "reflection_pdf: " << reflectP << endl;
        //     cout << "refraction_pdf: " << refractionPdf(bRec.wi, bRec.wo) << endl;
        //     cout << "pdf: " << pdf << endl;
        //     cout << "brdf: "
        //          << reflectVal << endl;
        //     cout << "btdf: "
        //          << btdf(bRec.wi, bRec.wo, Normal3f(0, 0, 1)) << endl;
        //     cout << "bsdf: "
        //          << bsdf << endl;
        //     cout << "val: "
        //          << val << endl;
        //     cout << "=============================================================" << endl;
        // }


    }

    bool isDiffuse() const {
        /* While microfacet BRDFs are not perfectly diffuse, they can be
           handled by sampling techniques for diffuse/non-specular materials,
           hence we return true here */
        return false;
    }

    std::string toString() const {
        return tfm::format(
            "Microfacet[\n"
            "  alpha = %f,\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_alpha,
            m_intIOR,
            m_extIOR);
    }

   private:
    float m_alpha;
    float m_intIOR, m_extIOR;

    float brdf(const Vector3f &wi, const Vector3f &wo, const Normal3f &n) const {
        Vector3f hr = (wi + wo).normalized();
        // if (wi.dot(n) < 0)
        if (hr.z() < 0)
            hr *= -1.f;

        float G = _G(wi, wo, hr, n);
        if (G == 0) return 0;

        float F = _F(wi, hr);

        float D = Warp::squareToGGXPdf(hr, m_alpha) / (hr.dot(n));  // compensate for the missing 1/cos(h)

        float val = F * G * D / (4.f * fabs(wi.dot(n)) * fabs(wo.dot(n)));

        return val;
    }

    float btdf(const Vector3f &wi, const Vector3f &wo, const Normal3f &n) const {

        /* implementation that takes reference to pbrt MicrofacetTransmission */
        if (wo.z() * wi.z() >= 0.f)  // wi and wo have to be in the different sides for transmission
            return 0.f;

        float cosThetaO = wo.dot(n);
        float cosThetaI = wi.dot(n);
        if (cosThetaI == 0 || cosThetaO == 0)
            return 0.f;

        float eta = wi.z() > 0 ? (m_intIOR / m_extIOR) : (m_extIOR / m_intIOR);

        // Compute the wh, note that wh lies in the same side as n
        Vector3f wh = (wi + wo * eta).normalized();
        if (wh.z() < 0)
            wh = -wh;

        if (wi.dot(wh) * wo.dot(wh) > 0)
            return 0.f;

        float F = _F(wi, wh);

        float sqrtDenom = wi.dot(wh) + eta * wo.dot(wh);
        // float factor    = 1.f; 
        float factor    = 1 / eta;
        float D = Warp::squareToGGXPdf(wh, m_alpha) / (wh.dot(n));  // compensate for the missing 1/cos(h)
        float G = _G(wi, wo, wh, n);

        return (1.f - F) *
               std::abs(D * G * eta * eta *
                        wi.dot(wh) * wo.dot(wh) * factor * factor /
                        (cosThetaI * cosThetaO * sqrtDenom * sqrtDenom));


        /* Implementation that takes reference to the Walter et al. paper */

        // float etaI = m_extIOR, etaO = m_intIOR;
        // if (wi.dot(n) < 0)
        //     std::swap(etaI, etaO);

        // Vector3f ht = -1.f * (etaI * wi + etaO * wo).normalized();
        // if (ht.z() < 0)
        //     ht *= -1;

        // float D = Warp::squareToGGXPdf(ht, m_alpha) / (ht.dot(n));  // compensate for the missing 1/cos(h)
        // if (D == 0) return 0.f;

        // float G = _G(wi, wo, ht, n);
        // if (G == 0) return 0.f;

        // float F = _F(wi, ht);

        // float val = fabs(wi.dot(ht)) * fabs(wo.dot(ht)) * etaO * etaO * (1 - F) * G * D;
        // val /= fabs(wi.dot(n)) * fabs(wo.dot(n)) * powf(etaI * wi.dot(ht) + etaO * wo.dot(ht), 2.0);

        // return val;
    }

    float _G(const Vector3f &wi, const Vector3f &wo, const Vector3f &m, const Normal3f &n) const {
        if (fabs(wi.z()) < Epsilon || fabs(wo.z()) < Epsilon)
            return 0.f;
        if ((wi.dot(m) * wi.dot(n)) <= 0)
            return 0.f;
        if ((wo.dot(m) * wo.dot(n)) <= 0)
            return 0.f;

        return G1(wi, m, n) * G1(wo, m, n);
    }

    float G1(const Vector3f &v, const Vector3f &m, const Normal3f &n) const {
        float tanT = Frame::tanTheta(v);
        return 2.f / (1.f + sqrt(1.f + m_alpha * m_alpha * tanT * tanT));
    }

    float _F(const Vector3f &wi, const Vector3f &m) const {
        return fresnel(wi.dot(m.z() > 0 ? m : -1 * m), m_extIOR, m_intIOR);
    }

    // from pbrt-v3 reflection.h line 97
    bool refract(const Vector3f &wi, const Normal3f &n, float eta, Vector3f &wt) const {
        // Compute $\cos \theta_\roman{t}$ using Snell's law
        float cosThetaI  = wi.dot(n);
        float sin2ThetaI = std::max(float(0), float(1 - cosThetaI * cosThetaI));
        float sin2ThetaT = eta * eta * sin2ThetaI;

        // Handle total internal reflection for transmission
        if (sin2ThetaT >= 1) return false;
        float cosThetaT = std::sqrt(1 - sin2ThetaT);
        wt              = eta * -wi + (eta * cosThetaI - cosThetaT) * Vector3f(n);
        wt.normalize();
        return true;
    }

    float refractionPdf(const Vector3f &wi, const Vector3f &wo) const {
        /* implementation that takes reference to pbrt MicrofacetTransmission */
        if (wi.z() * wo.z() > 0) return 0;

        float eta = wi.z() > 0 ? (m_intIOR / m_extIOR) : (m_extIOR / m_intIOR);


        Vector3f wh = (wi + wo * eta).normalized();

        if (wo.dot(wh) * wi.dot(wh) > 0) return 0;

        float sqrtDenom = wi.dot(wh) + eta * wo.dot(wh);
        float dwh_dwi   = std::abs((eta * eta * wo.dot(wh)) / (sqrtDenom * sqrtDenom));

        float D = Warp::squareToGGXPdf(wh.z() > 0 ? wh : -1 * wh, m_alpha);

        return D * dwh_dwi;
    }

    float reflectionPdf(const Vector3f &wi, const Vector3f &wo) const {
        if (wi.z() * wo.z() < 0) return 0;
        if ((wi + wo).norm() < Epsilon) return 0;
        Vector3f wh = (wo + wi).normalized();

        if (wh.z() < 0) wh *= -1;
        return Warp::squareToGGXPdf(wh, m_alpha) / fabs(4 * wo.dot(wh));
    }
};

NORI_REGISTER_CLASS(MicrofacetHacker, "microfacet-hacker");
NORI_NAMESPACE_END
