
#include <nori/common.h>
#include <nori/frame.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <nori/shape.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

// Implementation follows PBRT pp. 841-845
class Sphere : public Shape {
   public:
    Sphere(float radius, Point3f center) : m_radius(radius), m_center(center){};

    float samplePoint(const Point3f &ref, const Normal3f &ref_n, Sampler *sampler, Point3f &p, Normal3f &n, bool &over_area) const {
        Vector3f dir = m_center - ref;
        float dir_sn = dir.squaredNorm();
        float dir_n  = std::sqrt(dir_sn);
        if (dir_sn <= m_radius * m_radius) {  // ref is inside the sphere or on the surface
            over_area = true;
            return samplePoint(sampler, p, n);
        }

        over_area          = false;
        Point2f sample     = sampler->next2D();
        float sinThetaMax2 = m_radius * m_radius / dir_sn;
        float cosThetaMax  = std::sqrt(std::max(0.f, 1.f - sinThetaMax2));
        float cosTheta     = (1.f - sample(0)) + sample(0) * cosThetaMax;
        float sinTheta     = std::sqrt(std::max(0.f, 1 - cosTheta * cosTheta));
        float phi          = 2 * M_PI * sample(1);

        float ds       = dir_n * cosTheta - std::sqrt(std::max(0.f, m_radius * m_radius - dir_n * dir_n * sinTheta * sinTheta));
        float cosAlpha = (dir_n * dir_n + m_radius * m_radius - ds * ds) / (2 * dir_n * m_radius);
        float sinAlpha = std::sqrt(std::max(0.f, 1 - cosAlpha * cosAlpha));

        Frame fr = Frame(dir / dir_n);  // z: ref -> sphere center
        p        = m_center +
            (sinAlpha * std::cos(phi) * -fr.s +
             sinAlpha * std::sin(phi) * -fr.t +
             cosAlpha * -fr.n) *
                m_radius;
        n = (p - m_center).normalized();

        return 1.f / (2.f * M_PI * (1.f - cosThetaMax));
    }

    float pdf(const Point3f &ref, const Normal3f &ref_n, const Vector3f &wi, bool &over_area) const {
        Vector3f dir = m_center - ref;
        float dir_sn = dir.squaredNorm();
        if (dir_sn <= m_radius * m_radius) {  // ref is inside the sphere or on the surface
            over_area = true;
            return pdf();
        }

        over_area          = false;
        float sinThetaMax2 = m_radius * m_radius / dir_sn;
        float cosThetaMax  = std::sqrt(std::max(0.f, 1.f - sinThetaMax2));
        return 1.f / (2.f * M_PI * (1.f - cosThetaMax));
    };

    float samplePoint(Sampler *sampler, Point3f &p, Normal3f &n) const {
        Vector3f dir = Warp::squareToUniformSphere(sampler->next2D());
        p            = m_center + m_radius * dir;
        n            = (p - m_center).normalized();
        return pdf();
    }
    float pdf() const {
        return INV_PI / 4.0 / m_radius / m_radius;
    }

    bool rayIntersect(const Ray3f &ray) const {
        Vector3f dir = m_center - ray.o;

        float test = ray.d.normalized().dot(dir);
        test       = test * test;
        test -= dir.squaredNorm() - m_radius * m_radius;

        return test >= 0;
    }

    float getAreaSum() const {
        throw NoriException("Not implemented for spheres.");
    }

    std::string toString() const {
        return tfm::format(
            "Sphere[\n"
            "  radius = %f\n"
            "  center = %s\n"
            "]",
            m_radius,
            m_center.toString());
    }

   private:
    float m_radius;
    Point3f m_center;
};

class Disk : public Shape {
   public:
    Disk(float radius, Point3f center, Normal3f normal) : m_radius(radius), m_center(center), m_normal(normal) {
        local = Frame(normal);
    };

    float samplePoint(const Point3f &ref, const Normal3f &ref_n, Sampler *sampler, Point3f &p, Normal3f &n, bool &over_area) const {
        over_area = true;

        Point3f ref_l = local.toLocal(ref - m_center);
        ref_l(2)      = 0.f;  // project the point onto the disk plane

        float theta_ref;

        // corner cases for atan
        if (fabs(ref_l.x()) <= Epsilon) {
            if (ref_l.y() > 0)
                theta_ref = M_PI / 2;
            else
                theta_ref = 1.5 * M_PI;
        }

        if (ref_l(0) >= 0 && ref_l(1) >= 0)
            theta_ref = std::atan(ref_l(1) / ref_l(0));
        else if (ref_l(0) < 0 && ref_l(1) >= 0)
            theta_ref = M_PI - std::atan(ref_l(1) / -ref_l(0));
        else if (ref_l(0) < 0 && ref_l(1) < 0)
            theta_ref = M_PI + std::atan(ref_l(1) / ref_l(0));
        else
            theta_ref = 2 * M_PI - std::atan(-ref_l(1) / ref_l(0));

        Point2f u   = sampler->next2D();
        float t     = std::asin(2 * u(0) - 1);
        float theta = std::fmod(2 * t + theta_ref + 2 * M_PI, 2 * M_PI);
        float r     = std::sqrt(u(1)) * m_radius;

        p = Point3f(r * std::cos(theta), r * std::sin(theta), 0);
        p = local.toWorld(p) + m_center;
        n = m_normal;
        if ((ref - m_center).dot(m_normal) < 0)
            n = -n;

        float prob_r     = 2 / (m_radius * m_radius);
        float prob_theta = 0.25 * std::cos(t);

        return prob_r * prob_theta;
    }

    float pdf(const Point3f &ref, const Normal3f &ref_n, const Vector3f &wi, bool &over_area) const {
        throw NoriException("Not implemented for Disk.");
    };

    float samplePoint(Sampler *sampler, Point3f &p, Normal3f &n) const {
        Point2f s = Warp::squareToUniformDisk(sampler->next2D());
        p         = m_center + m_radius * local.toWorld(Vector3f(s(0), s(1), 0));
        n         = m_normal;  // not correct because it is not known which side the ref point lies on. But it is fine because I discard the test of angle of light source for disk
        return pdf();
    }

    float pdf() const {
        return INV_PI / m_radius / m_radius;
    }

    bool rayIntersect(const Ray3f &ray) const {
        float denominator = ray.d.normalized().dot(m_normal);
        if (std::fabs(denominator) <= Epsilon) return false;
        float d   = (m_center - ray.o).dot(m_normal) / denominator;
        Point3f p = ray.o + ray.d * d;
        if ((p - m_center).squaredNorm() <= m_radius * m_radius)
            return true;
        else
            return false;
    }

    float getAreaSum() const {
        throw NoriException("Not implemented for Disk.");
    }

    std::string toString() const {
        return tfm::format(
            "Disk[\n"
            "  radius = %f\n"
            "  center = %s\n"
            "]",
            m_radius,
            m_center.toString());
    }

   private:
    float m_radius;
    Point3f m_center;
    Normal3f m_normal;
    Frame local;
};


class Rectangle : public Shape {
   private:
// from PBRT pp.437
    void StratifiedSample2D(Point2f *samp, int nx, int ny, Sampler *sampler) const {
        float dx = (float)1 / nx, dy = (float)1 / ny;
        for (int y = 0; y < ny; ++y)
            for (int x = 0; x < nx; ++x) {
                float jx  = sampler->next1D();
                float jy  = sampler->next1D();
                samp->x() = std::min((x + jx) * dx, 1.f - Epsilon);
                samp->y() = std::min((y + jy) * dy, 1.f - Epsilon);
                ++samp;
            }
    }

   public:
    Rectangle(Point3f p, Vector3f e1, Vector3f e2, Normal3f n, int _nx, int _ny) : m_point(p), m_edge1(e1), m_edge2(e2), m_normal(n), nx(_nx), ny(_ny) {
        m_frame      = Frame(n);
        m_edge1_loc  = m_frame.toLocal(e1);
        m_edge2_loc  = m_frame.toLocal(e2);
        m_stratified = new Point2f[_nx * _ny];
    }
    ~Rectangle() {
        delete m_stratified;
    }

    // stratified sampling
    float samplePoint(const Point3f &ref, const Normal3f &ref_n, Sampler *sampler, Point3f &p, Normal3f &n, bool &over_area) const {
        if (m_stratified_idx < 0 || m_stratified_idx >= nx * ny) {
            StratifiedSample2D(m_stratified, nx, ny, sampler);
            m_stratified_idx = 0;
        }

        Point2f *ptr = m_stratified + m_stratified_idx;
        p            = m_point + m_edge1 * ptr->x() + m_edge2 * ptr->y();
        n            = m_normal;
        m_stratified_idx++;
        return pdf();
    }

    float pdf(const Point3f &ref, const Normal3f &ref_n, const Vector3f &wi, bool &over_area) const {
        throw NoriException("Not implemented for Rectangle.");
    };

    float samplePoint(Sampler *sampler, Point3f &p, Normal3f &n) const {
        Point2f u = sampler->next2D();
        p         = m_point + m_edge1 * u(0) + m_edge2 * u(1);
        n         = m_normal;
        return pdf();
    }

    float pdf() const {
        return 1.f / (m_edge1.norm() * m_edge2.norm());
    }

    bool rayIntersect(const Ray3f &ray) const {
        float denominator = ray.d.normalized().dot(m_normal);
        if (std::fabs(denominator) <= Epsilon)
            return false;
        float d   = (m_point - ray.o).dot(m_normal) / denominator;
        Point3f p = ray.o + ray.d.normalized() * d;

        Vector3f v = p - m_point;

        float proj1 = v.dot(m_edge1) / m_edge1.squaredNorm();
        if (proj1 <= 0 || proj1 >= 1)
            return false;

        float proj2 = v.dot(m_edge2) / m_edge2.squaredNorm();
        if (proj2 <= 0 || proj2 >= 1)
            return false;

        return true;
    }

    float getAreaSum() const {
        throw NoriException("Not implemented for Rectangle.");
    }

    std::string toString() const {
        return tfm::format(
            "Rectangle[\n"
            "]");
    }

   private:
    Point3f m_point;
    Vector3f m_edge1;
    Vector3f m_edge2;
    Normal3f m_normal;
    Frame m_frame;
    Vector3f m_edge1_loc;
    Vector3f m_edge2_loc;
    Point2f *m_stratified = nullptr;
    int nx = 4, ny = 4;
    mutable int m_stratified_idx = -1;
};

NORI_NAMESPACE_END