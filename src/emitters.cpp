
#include <nori/emitter.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <nori/shape.h>

#include "shapes.cpp"

NORI_NAMESPACE_BEGIN

bool better = false;
class SphereEmitter : public Emitter {
   public:
    // bool betterSampling = false;
    bool betterSampling = better;
    SphereEmitter(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
        m_shape    = new Sphere(props.getFloat("radius"), props.getPoint("position"));
    }

    Color3f sampleLi(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Sampler *sampler, const Scene *scene) const {
        Point3f p;
        Normal3f n;
        float pdf;
        bool over_area;
        if (betterSampling)
            pdf = m_shape->samplePoint(ref, ref_n, sampler, p, n, over_area);
        else
            pdf = m_shape->samplePoint(sampler, p, n);

        Vector3f dir   = ref - p;  // from light source to the reference point
        float dir_norm = dir.norm();
        Vector3f dir_u = dir / dir_norm;
        wi             = -1 * dir_u;

        // Test occlusion between p and ref
        Ray3f shadowRay(p, dir_u, Epsilon, dir_norm - Epsilon);
        if (scene->rayIntersect(shadowRay))
            // return Color3f(0.f, 1.f, 0.f);
            return Color3f(0.f);

        float cosTheta_p   = dir_u.dot(n);
        float cosTheta_ref = (-1 * dir_u).dot(ref_n);
        if (cosTheta_p <= 0 || cosTheta_ref <= 0)
            // return Color3f(1.f, 0, 0);
            return Color3f(0.0f);

        if (betterSampling) {
            return m_radiance * cosTheta_ref / pdf;
        } else {
            return m_radiance * abs(cosTheta_p * cosTheta_ref) / (dir.squaredNorm() * pdf);
        }
    }

    float pdfLi(const Point3f &ref, const Normal3f &ref_n, const Vector3f &wi) const {
        bool b;
        if (betterSampling)
            return m_shape->pdf(ref, ref_n, wi, b);  // wrt. solid angle
        else
            throw NoriException("Not implemented for Disk.");
        // return m_shape->pdf();  // wrt. area
    }
    Color3f le(const Point3f &p, const Normal3f &n, const Vector3f &wo) const {
        return m_radiance;
    }

    std::string toString() const {
        return tfm::format(
            "SphereEmitter[\n"
            "  radiance = %s\n"
            "]",
            m_radiance.toString());
    }
};

class DiskEmitter : public Emitter {
   public:
    // bool betterSampling = false;
    bool betterSampling = better;

    DiskEmitter(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
        m_shape    = new Disk(props.getFloat("radius"), props.getPoint("position"), props.getVector("normal"));
    }

    Color3f sampleLi(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Sampler *sampler, const Scene *scene) const {
        Point3f p;
        Normal3f n;
        float pdf = 0;
        bool b;
        if (betterSampling) {
            pdf = m_shape->samplePoint(ref, ref_n, sampler, p, n, b);
            while (pdf < 1e-2)
                pdf = m_shape->samplePoint(ref, ref_n, sampler, p, n, b);
        } else
            pdf = m_shape->samplePoint(sampler, p, n);

        Vector3f dir   = ref - p;  // from light source to the reference point
        float dir_norm = dir.norm();
        Vector3f dir_u = dir / dir_norm;
        wi             = -1 * dir_u;

        // Test occlusion between p and ref
        Ray3f shadowRay(p, dir_u, Epsilon, dir_norm - Epsilon);
        if (scene->rayIntersect(shadowRay))
            // return Color3f(0.f, 1.f, 0.f);
            return Color3f(0.f);

        float cosTheta_p   = dir_u.dot(n);
        float cosTheta_ref = (-1 * dir_u).dot(ref_n);
        // if (cosTheta_p <= 0 || cosTheta_ref <= 0)
        if (cosTheta_ref <= 0)
            // return Color3f(1.f, 0, 0);
            return Color3f(0.0f);

        return m_radiance * abs(cosTheta_p * cosTheta_ref) / (dir.squaredNorm() * pdf);
        // if (betterSampling) {
        //     return m_radiance * cosTheta_ref / pdf;
        // } else {
        //     return m_radiance * abs(cosTheta_p * cosTheta_ref) / (dir.squaredNorm() * pdf);
        // }
    }

    float pdfLi(const Point3f &ref, const Normal3f &ref_n, const Vector3f &wi) const {
        bool b;
        if (betterSampling)
            return m_shape->pdf(ref, ref_n, wi, b);  // wrt. solid angle
        else
            throw NoriException("Not implemented for Disk.");
        // return m_shape->pdf();  // wrt. area
    }

    Color3f le(const Point3f &p, const Normal3f &n, const Vector3f &wo) const {
        return m_radiance;
    }

    std::string toString() const {
        return tfm::format(
            "DiskEmitter[\n"
            "  radiance = %s\n"
            "]",
            m_radiance.toString());
    }
};

// class RectangleEmitter : public Emitter {
//    public:
//     // bool betterSampling = false;
//     bool betterSampling = better;

//     RectangleEmitter(const PropertyList &props) {
//         m_radiance = props.getColor("radiance");
//         m_shape    = new Rectangle(props.getPoint("point"), props.getVector("edge1"), props.getVector("edge2"), props.getVector("normal").normalized());
//     }

//     Color3f sampleLi(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Sampler *sampler, const Scene *scene) const {
//         Point3f p;
//         Normal3f n;
//         float pdf = 0;
//         bool b;
//         if (betterSampling) {
//             pdf = m_shape->samplePoint(ref, ref_n, sampler, p, n, b);
//             while (pdf < 1e-1)
//                 pdf = m_shape->samplePoint(ref, ref_n, sampler, p, n, b);
//         } else
//             pdf = m_shape->samplePoint(sampler, p, n);

//         Vector3f dir   = ref - p;  // from light source to the reference point
//         float dir_norm = dir.norm();
//         Vector3f dir_u = dir / dir_norm;
//         wi             = -1 * dir_u;

//         // Test occlusion between p and ref
//         Ray3f shadowRay(p, dir_u, Epsilon, dir_norm - Epsilon);
//         if (scene->rayIntersect(shadowRay))
//             return Color3f(0.f);

//         float cosTheta_p   = dir_u.dot(n);
//         float cosTheta_ref = (-1 * dir_u).dot(ref_n);
//         // if (cosTheta_p <= 0 || cosTheta_ref <= 0)
//         if (cosTheta_ref <= 0)
//             return Color3f(0.0f);

//         return m_radiance * abs(cosTheta_p * cosTheta_ref) / (dir.squaredNorm() * pdf);
//     }

//     float pdfLi(const Point3f &ref, const Normal3f &ref_n, const Vector3f &wi) const {
//         bool b;
//         if (betterSampling)
//             return m_shape->pdf(ref, ref_n, wi, b);  // wrt. solid angle
//         else
//             throw NoriException("Not implemented for Rectangle.");
//     }

//     Color3f le(const Point3f &p, const Normal3f &n, const Vector3f &wo) const {
//         return m_radiance;
//     }

//     std::string toString() const {
//         return tfm::format(
//             "RectangleEmitter[\n"
//             "  radiance = %s\n"
//             "]",
//             m_radiance.toString());
//     }
// };

class RectangleEmitter : public Emitter {
   public:
    // bool betterSampling = false;
    bool betterSampling = better;
    int nx              = 4;
    int ny              = 4;

    RectangleEmitter(const PropertyList &props) {
        m_radiance = props.getColor("radiance");
        nx         = props.getInteger("nx");
        ny         = props.getInteger("ny");
        m_shape    = new Rectangle(props.getPoint("point"), props.getVector("edge1"), props.getVector("edge2"), props.getVector("normal").normalized(), nx, ny);
    }

    Color3f sampleLi(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Sampler *sampler, const Scene *scene) const {
        Point3f p;
        Normal3f n;
        float pdf = 0;
        bool b;
        if (betterSampling) {
            Color3f res(0.f);
            for (int i = 0; i < nx * ny; i++) {
                pdf            = m_shape->samplePoint(ref, ref_n, sampler, p, n, b);
                Vector3f dir   = ref - p;  // from light source to the reference point
                float dir_norm = dir.norm();
                Vector3f dir_u = dir / dir_norm;
                wi             = -1 * dir_u;

                // Test occlusion between p and ref
                Ray3f shadowRay(p, dir_u, Epsilon, dir_norm - Epsilon);
                if (scene->rayIntersect(shadowRay))
                    res += Color3f(0.f);
                else {
                    float cosTheta_p   = dir_u.dot(n);
                    float cosTheta_ref = (-1 * dir_u).dot(ref_n);
                    if (cosTheta_ref <= 0)
                        res += Color3f(0.f);
                    else
                        res += m_radiance * abs(cosTheta_p * cosTheta_ref) / (dir.squaredNorm() * pdf);
                }
            }
            res /= float(nx * ny);
            return res;

        } else {
            pdf = m_shape->samplePoint(sampler, p, n);

            Vector3f dir   = ref - p;  // from light source to the reference point
            float dir_norm = dir.norm();
            Vector3f dir_u = dir / dir_norm;
            wi             = -1 * dir_u;

            // Test occlusion between p and ref
            Ray3f shadowRay(p, dir_u, Epsilon, dir_norm - Epsilon);
            if (scene->rayIntersect(shadowRay))
                return Color3f(0.f);

            float cosTheta_p   = dir_u.dot(n);
            float cosTheta_ref = (-1 * dir_u).dot(ref_n);
            // if (cosTheta_p <= 0 || cosTheta_ref <= 0)
            if (cosTheta_ref <= 0)
                return Color3f(0.0f);

            return m_radiance * abs(cosTheta_p * cosTheta_ref) / (dir.squaredNorm() * pdf);
        }
    }

    float pdfLi(const Point3f &ref, const Normal3f &ref_n, const Vector3f &wi) const {
        bool b;
        if (betterSampling)
            return m_shape->pdf(ref, ref_n, wi, b);  // wrt. solid angle
        else
            throw NoriException("Not implemented for Rectangle.");
    }

    Color3f le(const Point3f &p, const Normal3f &n, const Vector3f &wo) const {
        return m_radiance;
    }

    std::string toString() const {
        return tfm::format(
            "RectangleEmitter[\n"
            "  radiance = %s\n"
            "]",
            m_radiance.toString());
    }
};

NORI_REGISTER_CLASS(SphereEmitter, "sphere");
NORI_REGISTER_CLASS(DiskEmitter, "disk");
NORI_REGISTER_CLASS(RectangleEmitter, "rectangle");
NORI_NAMESPACE_END