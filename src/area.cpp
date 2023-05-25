
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <nori/shape.h>
#include <nori/texture.h>
#include <math.h>

NORI_NAMESPACE_BEGIN

    class AreaLight : public Emitter {
private:
    mutable Mesh* m_mesh = nullptr;
    public:
        AreaLight(const PropertyList &props) {
            m_radiance = props.getColor("radiance");
            float strength = props.getFloat("strength", 1.f);
            m_radiance *= strength;
        }

        float pdfLi(const Point3f &ref, const Normal3f &ref_n, const Point3f &p, const Normal3f &n) const {
            float pdf = m_shape->pdf();
            Vector3f dir = ref - p;  // from light source to the reference point
            Vector3f dir_u = dir.normalized();
            float cosTheta_p = dir_u.dot(n);

            if (cosTheta_p < 0.f)
                return 0;

            // return pdf is wrt. solid angle
            return pdf * dir.squaredNorm() / cosTheta_p;
        }

        Color3f
        sampleLi(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Point3f &pp, Normal3f &nn, Sampler *sampler,
                 const Scene *scene) const {
            /* For ref within medium, ref_n is Normal3f() */

            Point3f p;
            Normal3f n;
            Point2f uv;
            m_shape->samplePoint(sampler, p, n, uv);
            pp = p;
            nn = n;
            Vector3f dir = ref - p;  // from light source to the reference point
            float dir_norm = dir.norm();
            Vector3f dir_u = dir / dir_norm;
            wi = -1 * dir_u;

            // Test occlusion between p and ref for non-medium
            Ray3f shadowRay(p, dir_u, Epsilon, dir_norm - Epsilon);
            if (scene->rayIntersect(shadowRay))
                return Color3f(0, 0, 0);

            float cosTheta_p = dir_u.dot(n);
            float cosTheta_ref;

            cosTheta_ref = (-1 * dir_u).dot(ref_n);

            if (cosTheta_p <= 0 || cosTheta_ref <= 0)
                return Color3f(0, 0, 0);

            float pdf = 1.f / m_shape->getAreaSum();

            if (!m_mesh)
                m_mesh = static_cast<Mesh *>(m_shape);

            if (m_mesh->getEmitterTexture() != nullptr) {
                Vector3f c = m_mesh->getEmitterTexture()->queryUV(uv);
                return Color3f(c.x(), c.y(), c.z()) * abs(cosTheta_p * cosTheta_ref) / (dir.squaredNorm() * pdf);
            } else
                return m_radiance * abs(cosTheta_p * cosTheta_ref) / (dir.squaredNorm() * pdf);
        }

        Color3f
        sampleLi_no_occlusion(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Point3f &pp, Normal3f &nn,
                              Sampler *sampler, const Scene *scene) {
            /* For ref within medium, ref_n is Normal3f() */

            Point3f p;
            Normal3f n;
            Point2f uv;
            m_shape->samplePoint(sampler, p, n, uv);
            pp = p;
            nn = n;
            Vector3f dir = ref - p;  // from light source to the reference point
            float dir_norm = dir.norm();
            Vector3f dir_u = dir / dir_norm;
            wi = -1 * dir_u;


            float cosTheta_p = dir_u.dot(n);
            float cosTheta_ref;

            // for mediumInteraction, there is no ref_n
            if (ref_n.squaredNorm() == 0.f)
                cosTheta_ref = 1.f;
            else
                cosTheta_ref = (-1 * dir_u).dot(ref_n);

            if (cosTheta_p <= 0 || cosTheta_ref <= 0)
                return Color3f(0, 0, 0);

            float pdf = 1.f / m_shape->getAreaSum();

            if (!m_mesh)
                m_mesh = static_cast<Mesh *>(m_shape);

            if (m_mesh->getEmitterTexture() != nullptr) {
                Vector3f c = m_mesh->getEmitterTexture()->queryUV(uv);
                return Color3f(c.x(), c.y(), c.z()) * abs(cosTheta_p * cosTheta_ref) / (dir.squaredNorm() * pdf);
            } else
                return m_radiance * abs(cosTheta_p * cosTheta_ref) / (dir.squaredNorm() * pdf);
        }

        Color3f le(const Point3f &p, const Normal3f &n, const Vector3f &wo) const {

            if (n.dot(wo) > 0)
                return m_radiance;
            else
                return Color3f(0.f);
        }

        Color3f le(const Point3f &p, const Normal3f &n, const Vector3f &wo, const Color3f &texture) const {
            if (!m_mesh)
                m_mesh = static_cast<Mesh *>(m_shape);

            if (n.dot(wo) > 0) {
                if (m_mesh->getEmitterTexture() != nullptr) {
                    return texture;
                } else
                    return m_radiance;
            }
            else
                return Color3f(0.f);
        }

        bool isDeltaLight() const {
            return false;
        }

        std::string toString() const {
            return tfm::format(
                    "AreaLight[\n"
                    "  radiance = %s\n"
                    "]",
                    m_radiance.toString());
        }
    };

    NORI_REGISTER_CLASS(AreaLight, "area");
NORI_NAMESPACE_END