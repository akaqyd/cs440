
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <nori/shape.h>
#include <math.h>

NORI_NAMESPACE_BEGIN

    class DirectionalLight : public Emitter {
    private:
        Vector3f m_dir;
        mutable Point3f m_scene_center;
        mutable float m_scene_radius = 0.f;
    public:
        DirectionalLight(const PropertyList &props) {
            m_radiance = props.getColor("radiance");
            float strength = props.getFloat("strength", 1.);
            m_radiance *= strength;

            float theta = props.getFloat("theta", -1);
            float phi = props.getFloat("phi", -1);
            if (theta != -1 && phi != -1) {
                theta = theta / 180.f * M_PI;
                phi = phi / 180.f * M_PI;
                m_dir(0) = std::cos(phi) * std::cos(theta);
                m_dir(1) = -1.f * std::sin(theta);
                m_dir(2) = std::sin(phi) * std::cos(theta);
                m_dir.normalize();
            }
            else
                m_dir = props.getVector("direction").normalized();
        }

        float pdfLi(const Point3f &ref, const Normal3f &ref_n, const Point3f &p, const Normal3f &n) const override {
            // from PBRT
//            return 0.f;

            return 1.f;
        }

        Color3f
        sampleLi(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Point3f &pp, Normal3f &nn, Sampler *sampler,
                 const Scene *scene) const override {
            /* For ref within medium, ref_n is Normal3f() */

            // not initialized
            if (m_scene_radius == 0.f) {
                BoundingBox3f scene_bbox = scene->getBoundingBox();
                m_scene_center = scene_bbox.getCenter();
                m_scene_radius = (m_scene_center - scene_bbox.min).norm();
            }

            wi = -1 * m_dir;

            pp = ref + 2.f * m_scene_radius * wi;
            nn = m_dir;

            // Test occlusion between p and ref for non-medium
            Ray3f shadowRay(ref, wi, Epsilon, 2.f * m_scene_radius);
            if (scene->rayIntersect(shadowRay))
                return Color3f(0, 0, 0);

            float cosTheta_ref = wi.dot(ref_n);

            if (cosTheta_ref <= 0)
                return Color3f(1, 0, 0);

            float pdf = 1.f;

            return m_radiance * abs(cosTheta_ref) / pdf;
        }

        Color3f
        sampleLi_no_occlusion(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Point3f &pp, Normal3f &nn,
                              Sampler *sampler, const Scene *scene) override {
            /* For ref within medium, ref_n is Normal3f() */

            // not initialized
            if (m_scene_radius == 0.f) {
                BoundingBox3f scene_bbox = scene->getBoundingBox();
                m_scene_center = scene_bbox.getCenter();
                m_scene_radius = (m_scene_center - scene_bbox.min).norm();
            }

            wi = -1 * m_dir;

            pp = ref + 2.f * m_scene_radius * wi;
            nn = m_dir;

            float cosTheta_ref;
            // for mediumInteraction, there is no ref_n
            if (ref_n.squaredNorm() == 0.f)
                cosTheta_ref = 1.f;
            else
                cosTheta_ref = wi.dot(ref_n);

            if (cosTheta_ref <= 0)
                return Color3f(0, 0, 0);

            float pdf = 1.f;

            return m_radiance * abs(cosTheta_ref) / pdf;
        }

        Color3f le(const Point3f &p, const Normal3f &n, const Vector3f &wo) const override {
            return 0.f;
        }

        bool isDeltaLight() const override {
            return true;
        }

        std::string toString() const override {
            return tfm::format(
                    "DirectionalLight[\n"
                    "  radiance = %s\n"
                    "  direction = %s\n"
                    "]",
                    m_radiance.toString(), m_dir.toString());
        }
    };

    NORI_REGISTER_CLASS(DirectionalLight, "directional");
NORI_NAMESPACE_END