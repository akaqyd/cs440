//
// Created by Qiyuan Dong on 02.05.22.
//

#include <nori/medium.h>
#include <nori/transform.h>
#include <iostream>
#include <nori/bbox.h>
#include <Eigen/LU>

NORI_NAMESPACE_BEGIN
/**
 * Implementation of image based lighting
 * My design and implementation of IBL make references to PBRT's grid.cpp.
 */

    class Heterogeneous : public Medium {
    private:
        Color3f m_sigma_a, m_sigma_s;
        float m_sigma_t;    // Implementation requires m_sigma_t is the same for R, G, and B
        float m_g;
        VDB m_vdb;
        Transform m_MediumToWorld;
        Transform m_WorldToMedium;

    public:
        Heterogeneous(const PropertyList &propList) {
            float density_a = propList.getFloat("density_a", 1);
            float density_s = propList.getFloat("density_s", 1);
            m_sigma_a = propList.getColor("sigma_a", Color3f(0)) * density_a;
            m_sigma_s = propList.getColor("sigma_s", Color3f(0)) * density_s;
            Color3f sigma_t = m_sigma_a + m_sigma_s;
            assert(sigma_t.x() == sigma_t.y() && sigma_t.x() == sigma_t.z());
            m_sigma_t = sigma_t.x();

            m_g = propList.getFloat("g", 0);

            std::string fp = propList.getString("vdb_path");
            m_vdb.reset(fp);

            m_MediumToWorld = propList.getTransform("MediumToWorld");
            m_WorldToMedium = Transform(m_MediumToWorld.getMatrix().inverse());

//            m_vdb.transform(m_MediumToWorld);
        }

        Color3f transmittance(const Ray3f &_ray, Sampler *sampler) const override {
            // Ensure _ray.d is normalized
            assert(std::fabs(_ray.d.norm() - 1.f) < Epsilon);

            // Convert ray from the world coordinates into medium local coordinates
            // The grid's extent is [0, 1]^3 in its local coordinates
            // Note that ray_loc.d is NOT normalized

//            Ray3f ray_loc(_ray);

            Ray3f ray_loc;
            ray_loc.o = m_WorldToMedium * _ray.o;
            ray_loc.d = m_WorldToMedium * _ray.d;

//            float scale = ray_loc.d.norm() / _ray.d.norm();
//            ray_loc.d /= scale;
            ray_loc.mint = Epsilon;
//            ray_loc.maxt = _ray.maxt * scale;
//            ray_loc.maxt = _ray.maxt / scale;
            ray_loc.maxt = _ray.maxt;
            ray_loc.m = _ray.m;
            ray_loc.update();

            float threshold = 0.1;
            Color3f Tr(1.f);
            float tMin = ray_loc.mint, tMax = ray_loc.maxt;
            float _tMin, _tMax;
            BoundingBox3f gridBbox = m_vdb.get_bbox();
            if (!gridBbox.rayIntersect(ray_loc, _tMin, _tMax))
                return Tr;

            tMin = fmax(tMin, _tMin);
            tMax = fmin(tMax, _tMax);

            float t = tMin;

            while (true) {
                Point2f sample = sampler->next2D();
                float dt = -1.f * std::log(1 - sample[0]) * m_vdb.get_inv_max_density() / m_sigma_t;
                t += dt;
                if (t >= tMax) break;   // exit the grid
                Tr *= 1.f - std::max(0.f, m_vdb.queryDensity(ray_loc(t)) * m_vdb.get_inv_max_density());

                if (Tr.minCoeff() < threshold) {
                    float q = 0.05;
                    if (1.f - Tr.minCoeff() > q)
                        q = 1.f - Tr.minCoeff();
                    if (sample[1] < q)
                        return 0.f;
                    Tr /= 1.f - q;
                }
            }
            return Tr;
        }

        Color3f sample(const Ray3f &_ray, Sampler *sampler, MediumInteraction *mi) const override {
            // Ensure _ray.d is normalized
            assert(std::fabs(_ray.d.norm() - 1.f) < Epsilon);

            // Convert ray from the world coordinates into medium local coordinates
            // The grid's extent is [0, 1]^3 in its local coordinates
            // Note that ray_loc.d is still normalized
//            Ray3f ray_loc(_ray);
            Ray3f ray_loc;
            ray_loc.o = m_WorldToMedium * _ray.o;
            ray_loc.d = m_WorldToMedium * _ray.d;

//            float scale = ray_loc.d.norm() / _ray.d.norm();
//            ray_loc.d /= scale;
            ray_loc.mint = Epsilon;
//            ray_loc.maxt = _ray.maxt * scale;
//            ray_loc.maxt = _ray.maxt / scale;
            ray_loc.maxt = _ray.maxt;
            ray_loc.m = _ray.m;
            ray_loc.update();

//            std::cout << "==========================" << std::endl;
//            std::cout << "Heterogeneous sample(), worldRay: " << _ray.toString() << ", localRay: " << ray_loc.toString() << std::endl;

            float t = 0, tMin = ray_loc.mint, tMax = ray_loc.maxt;
            float _tMin, _tMax;
            BoundingBox3f gridBbox = m_vdb.get_bbox();
            if (!gridBbox.rayIntersect(ray_loc, _tMin, _tMax)) {
//                std::cout << "No intersection with the grid bbox." << std::endl;
                return 1.f;
            }
            tMin = fmax(tMin, _tMin);
            tMax = fmin(tMax, _tMax);
//

//            std::cout << "Intersection with the grid bbox, tMin " << tMin << ", tMax " << tMax << std::endl;

            t = tMin;
            while (true) {
                Point2f sample = sampler->next2D();
                float dt = -1.f * std::log(1 - sample[0]) * m_vdb.get_inv_max_density() / m_sigma_t;
                t += dt;
//                std::cout << "dt: " << dt << ", t: " << t << std::endl;
                if (t >= tMax)
                    break;   // exit the grid
                if (m_vdb.queryDensity(ray_loc(t)) * m_vdb.get_inv_max_density() > sample[1]) {
                    *mi = MediumInteraction(_ray(t), -_ray.d, this);
//                    *mi = MediumInteraction(_ray(t / scale), -_ray.d, this);
//                    *mi = MediumInteraction(_ray(t * scale), -_ray.d, this);
//                    std::cout << "Sampled an interaction, t: " << t << ", value: " << (m_sigma_s / m_sigma_t) <<  std::endl;
                    return m_sigma_s / m_sigma_t;
                }
            }

            return 1.f;
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
                    "Heterogeneous [\n"
                    "  m_sigma_a = %s\n"
                    "  m_sigma_s = %s\n"
                    "  m_sigma_t = %f\n"
                    "  m_g = %f\n"
                    "  m_WorldToMedium = %s\n"
                    "  m_MediumToWorld= %s\n"
                    "]", m_sigma_a.toString(), m_sigma_s.toString(), m_sigma_t, m_g,
                    m_WorldToMedium.toString(), m_MediumToWorld.toString());
        }

    };

    NORI_REGISTER_CLASS(Heterogeneous, "heterogeneous");

NORI_NAMESPACE_END