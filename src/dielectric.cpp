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
#include <nori/common.h>
#include <nori/frame.h>

NORI_NAMESPACE_BEGIN

/// Ideal dielectric BSDF
class Dielectric : public BSDF {
   public:
    Dielectric(const PropertyList &propList) {
        /* Interior IOR (default: BK7 borosilicate optical glass) */
        m_intIOR = propList.getFloat("intIOR", 1.5046f);

        /* Exterior IOR (default: air) */
        m_extIOR = propList.getFloat("extIOR", 1.000277f);
    }

    Color3f eval(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return Color3f(0.0f);
    }

    float pdf(const BSDFQueryRecord &) const {
        /* Discrete BRDFs always evaluate to zero in Nori */
        return 0.0f;
    }

    Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
        float cosThetaT;
        float Fr = fresnel(bRec.wi[2], m_extIOR, m_intIOR, cosThetaT);

        if (sample[0] <= Fr) {
            // Reflection in local coordinates
            bRec.wo = Vector3f(
                -bRec.wi.x(),
                -bRec.wi.y(),
                bRec.wi.z());
            bRec.measure = EDiscrete;

            bRec.eta = 1.0f;

            return Color3f(1.0f);
        } else {
            float sinThetaT = sqrt(1.0f - cosThetaT * cosThetaT);
            float sinTheta = sqrt(1.0f - bRec.wi.z() * bRec.wi.z());
            float wo_x = -1 * bRec.wi.x() * (sinThetaT / sinTheta);
            float wo_y = -1 * bRec.wi.y() * (sinThetaT / sinTheta);
            float wo_z = bRec.wi.z() < 0 ? cosThetaT : -cosThetaT;
            bRec.wo = Vector3f(wo_x, wo_y, wo_z);
            bRec.measure = EDiscrete;
            // bRec.eta = bRec.wi.z() > 0 ? (m_intIOR / m_extIOR) : (m_extIOR / m_intIOR);
            bRec.eta = bRec.wi.z() > 0 ? (m_extIOR / m_intIOR) : (m_intIOR / m_extIOR);

            // return Color3f(1.0f);
            return Color3f(1.0f) * bRec.eta * bRec.eta;     // multiply the correction term
        }
    }

    std::string toString() const {
        return tfm::format(
            "Dielectric[\n"
            "  intIOR = %f,\n"
            "  extIOR = %f\n"
            "]",
            m_intIOR, m_extIOR);
    }

   private:
    float m_intIOR, m_extIOR;
};

NORI_REGISTER_CLASS(Dielectric, "dielectric");
NORI_NAMESPACE_END
