//
// Created by Qiyuan Dong on 02.05.22.
//

#include <nori/medium.h>

NORI_NAMESPACE_BEGIN
    class Homogeneous : public Medium {
    private:
        Color3f m_sigma_a, m_sigma_s, m_sigma_t;
        float m_g;

    public:
        Homogeneous(const PropertyList &propList) {
            m_g = propList.getFloat("g", 0);
            float density_a = propList.getFloat("density_a", 1);
            float density_s = propList.getFloat("density_s", 1);
            m_sigma_a = propList.getColor("sigma_a", Color3f(0)) * density_a;
            m_sigma_s = propList.getColor("sigma_s", Color3f(0)) * density_s;
            m_sigma_t = m_sigma_a + m_sigma_s;
        }

        Color3f transmittance(const Ray3f &ray, Sampler *sampler) const override {
            float dist = ray.d.norm() * ray.maxt;
//            if (dist >dist= Scene_Boundary)
//                throw (NoriException("Homogeneous Transmittance computation: ray.d.norm() * ray.maxt >= SCENE_MAX"));
//            else
                return {std::exp(-m_sigma_t.x() * dist),
                        std::exp(-m_sigma_t.y() * dist),
                        std::exp(-m_sigma_t.z() * dist)};
        }

        Color3f sample(const Ray3f &ray, Sampler *sampler, MediumInteraction *mi) const override {
            /* Ensure that the passed ray.d is always normalized */

            Point2f sample = sampler->next2D();
            int c = sample.x() <= 0.3333 ? 0 : (sample.x() >= 0.6667 ? 2 : 1);
            float sampled_t = -1.f * std::log(1.f - sample.y()) / m_sigma_t(c);
            sampled_t = std::min(sampled_t, ray.maxt);    // sampled t must stay in the ray interval
            bool reachSurface = std::fabs(sampled_t - ray.maxt) <= Epsilon;

            Color3f tr{std::exp(-m_sigma_t.x() * sampled_t),
                       std::exp(-m_sigma_t.y() * sampled_t),
                       std::exp(-m_sigma_t.z() * sampled_t)};

            float pdf_surface = tr.sum() / 3.f;
            float pdf_medium = (m_sigma_t * tr).sum() / 3.f;

            if (reachSurface) {
                if (pdf_surface == 0.f)
                    throw NoriException("Homogeneous pdf_surface == 0.f");
                return tr / pdf_surface;
            } else {
                if (pdf_medium == 0.f)
                    throw NoriException("Homogeneous pdf_medium == 0.f");

                *mi = MediumInteraction(ray(sampled_t), -ray.d, this);

                return tr * m_sigma_s / pdf_medium;
            }
        }

        /* Use HG as the default phase function */
        float sample_pf(const Vector3f &wi, Vector3f &wo, Sampler *sampler) const override {
            Point2f sample = sampler->next2D();
            float cosTheta, sinTheta, phi;
            if (std::abs(m_g) < 0.001)
                cosTheta = 1.f - 3.f * sample[0];
            else {
                float t = (1 - m_g * m_g) / (1 + m_g - 2 * m_g * sample[0]);
                cosTheta = -(1 + m_g * m_g - t * t) / (2 * m_g);
            }
            sinTheta = std::sqrt(std::max(0.f, 1 - cosTheta * cosTheta));
            phi = 2 * M_PI * sample[1];

            Frame loc(wi);
            wo = sinTheta * std::cos(phi) * loc.s + sinTheta * std::sin(phi) * loc.t +
                 cosTheta * loc.n;

            float d = 1.f + m_g * m_g + 2.f * m_g * cosTheta;
            return INV_PI / 4.0 * (1.f - m_g * m_g) / (d * std::sqrt(d));
        }

        [[nodiscard]] float eval_pf(const Vector3f &wi, const Vector3f &wo) const override {
            float cosTheta = wi.dot(wo);
            float d = 1.f + m_g * m_g + 2.f * m_g * cosTheta;
            return INV_PI / 4.0 * (1.f - m_g * m_g) / (d * std::sqrt(d));
        }


        [[nodiscard]] std::string toString() const override {
            return tfm::format(
                    "Homogeneous [\n"
                    "  m_sigma_a = %s\n"
                    "  m_sigma_s = %s\n"
                    "  m_sigma_t = %s\n"
                    "  m_g = %f\n"
                    "]", m_sigma_a.toString(), m_sigma_s.toString(), m_sigma_t.toString(), m_g);
        }

    };

NORI_REGISTER_CLASS(Homogeneous, "homogeneous");

NORI_NAMESPACE_END
