#include <math.h>
#include <nori/bitmap.h>
#include <nori/common.h>
#include <nori/emitter.h>
#include <nori/frame.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <nori/shape.h>
#include <nori/warp.h>
#include <filesystem/resolver.h>
#include <nori/transform.h>
#include <Eigen/LU>

#include <string>

NORI_NAMESPACE_BEGIN

    using namespace std;

/**
 * Implementation of image based lighting
 * My design and implementation of IBL make some reference to PBRT's Infinite Light.
 */
    class IBLLight : public Emitter {
    public:
        IBLLight(const PropertyList &props) {
            Eigen::Matrix4f defaultTM;
            defaultTM <<
                      0, -1, 0, 0,
                    0, 0, 1, 0,
                    -1, 0, 0, 0,
                    0, 0, 0, 1;
            Transform defaultT(defaultTM, defaultTM.inverse());

            m_exr_path = getFileResolver()->resolve(props.getString("exr_path")).str();
            m_strength = props.getFloat("strength", 1.f);
            m_toWorld = props.getTransform("toWorld", defaultT);
            m_toLight = Transform(m_toWorld.getMatrix().inverse());

            Warp::buildMipMap(m_exr_path);
            m_bm = Bitmap(m_exr_path);

            // debug
            // m_bm.transposeInPlace();
            // cout << "IBL map 0,0 : " << endl
            //      << m_bm(0, 0) << endl;
            // cout << "IBL map 0,1 : " << endl
            //      << m_bm(0, 1) << endl;
            // cout << "IBL map 1,0 : " << endl
            //      << m_bm(1, 0) << endl;
            // cout << "IBL map 1,1 : " << endl
            //      << m_bm(1, 1) << endl;

            m_rows = m_bm.rows();
            m_cols = m_bm.cols();
        }

        Point2f sphere2uv(Vector3f wi) const {
            Point2f theta_phi = sphericalCoordinates(wi);
            float theta = theta_phi[0], phi = theta_phi[1];
            return Point2f(phi / 2 / M_PI, 1.f - theta / M_PI);
        }

        Vector3f uv2sphere(Point2f uv) const {
            float theta = (1.f - uv[1]) * M_PI, phi = uv[0] * 2 * M_PI;
            // return sphericalDirection(theta, phi);
            float cosTheta = std::cos(theta), sinTheta = std::sin(theta);
            float sinPhi = std::sin(phi), cosPhi = std::cos(phi);
            return Vector3f(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
        }

        /**
         * @brief Given a uv coordinates in [0, 1]^2, return the interpolated value from the bitmap
         *
         * @param bm The bitmap
         * @param uv Coordinates in [0, 1] ^ 2, where the origin is left bottom corner, uv[0] is x, uv[1] is y
         * @return Color3f Interpolated value.
         */
        Color3f lookup(const Bitmap &bm, const Point2f &uv) const {
            assert(uv[0] < 1 && uv[1] < 1);
            assert(uv[0] > 0 && uv[1] > 0);
            float r  = (1.f - uv[1]) * (float)m_rows;
            float c  = uv[0] * (float)m_cols;
            float r0, r1, c0, c1;

            r0 = std::floor(r - 0.5);
            c0 = std::floor(c - 0.5);
            r1 = r0 + 1.f;
            c1 = c0 + 1.f;

            // handle the border cases
            if (r < 0.5 || r > ((float)m_rows - 0.5)) {
                r0 = (float) m_rows - 1.f;
                r1 = 0;
            }
            if (c < 0.5 || c > ((float)m_cols - 0.5)) {
                c0 = (float) m_cols - 1.f;
                c1 = 0;
            }

            Color3f l00(0.f), l01(0.f), l10(0.f), l11(0.f);
            if (r0 >= 0 && r0 < (float)m_rows && c0 >= 0 && c0 < (float)m_cols) l00 = bm(r0, c0);
            if (r0 >= 0 && r0 < (float)m_rows && c1 >= 0 && c1 < (float)m_cols) l01 = bm(r0, c1);
            if (r1 >= 0 && r1 < (float)m_rows && c0 >= 0 && c0 < (float)m_cols) l10 = bm(r1, c0);
            if (r1 >= 0 && r1 < (float)m_rows && c1 >= 0 && c1 < (float)m_cols) l11 = bm(r1, c1);

            // bilinear interpolation
//        Color3f tmp_r0 = (c - c0) * l01 + (c1 - c) * l00;
//        Color3f tmp_r1 = (c - c0) * l11 + (c1 - c) * l10;

            Color3f tmp_r0 = fmod(c - c0 - 0.5 + (float)m_cols, (float)m_cols) * l01 + fmod(c1 - c + 0.5 + (float)m_cols, (float)m_cols) * l00;
            Color3f tmp_r1 = fmod(c - c0 - 0.5 + (float)m_cols, (float)m_cols) * l11 + fmod(c1 - c + 0.5 + (float)m_cols, (float)m_cols) * l10;

//        Color3f val = (r - r0) * tmp_r1 + (r1 - r) * tmp_r0;
            Color3f val = fmod(r - r0 - 0.5 + (float)m_rows, (float)m_rows) * tmp_r1 + fmod(r1 - r + 0.5 + (float)m_rows, (float)m_rows) * tmp_r0;

            // cout << "uv: " << endl << uv << endl;
            // cout << "rows: " << m_rows << ", cols: " << m_cols << endl;
            // cout << "r: " << r << ", c: " << c << ", r0: " << r0 << ", r1: " << r1 << ", c0: " << c0 << ", c1:" << c1 << endl;
            // cout << "l00: " << l00.x() << ", l01: " << l01.x() << ", l10: " << l10.x() << ", l11: " << l11.x() << endl;
            // cout << "tmp_r0: " << tmp_r0.x() << ", tmp_r1: " << tmp_r1.x() << endl;
            // cout << "val: " << val.x() <<  endl;

            return val * m_strength;
        }

        float pdfLi(const Point3f &ref, const Normal3f &ref_n, const Point3f &p, const Normal3f &n) const override {
            Vector3f wi    = p.normalized();
            Point2f uv     = sphere2uv(wi);
            float jacobian = 2 * M_PI * M_PI * Frame::sinTheta(wi);
//            float pdf = Warp::squareToHierarchicalPdf(uv) * m_rows * m_cols / jacobian;
            float pdf = Warp::squareToHierarchicalPdf(uv) / jacobian;
            if (pdf < Epsilon || isnan(pdf) || isnan(pdf))
                return 0.f;
            else
                return pdf;
        }

        Color3f sampleLi(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Point3f &pp, Normal3f &nn, Sampler *sampler, const Scene *scene) const override {
            Point2f uv = Warp::squareToHierarchical(sampler->next2D());
            Vector3f loc_wi         = uv2sphere(uv).normalized();
            wi = m_toWorld * loc_wi;
            pp         = wi * Scene_Boundary;
            nn         = -wi;

            float cosTheta = ref_n.dot(wi);
            if (cosTheta < 0.0f)
                return 0.f;

            // test visibility
            Ray3f shadowRay(ref, wi, Epsilon, Scene_Boundary);
            if (scene->rayIntersect(shadowRay))
                return 0.f;

            Color3f bm_val = lookup(m_bm, uv);
            float pdf      = pdfLi(ref, ref_n, loc_wi, nn);
            Color3f val    = bm_val * cosTheta / pdf;

            return isnan(val.x()) || isinf(val.x()) ? 0.f : val;
        }

        Color3f le(const Point3f &p, const Normal3f &n, const Vector3f &wo) const override{
            return 0.f;
        }

        Color3f le(const Ray3f &ray) const override {
            Vector3f wi =  ray.d.normalized();
            wi = m_toLight * wi;
            Point2f uv  = sphere2uv(wi);
            Color3f val = lookup(m_bm, uv);

            // cout << "----------------------------------------------------------------" << endl;
            // cout << "wi: " << wi.x() << " " << wi.y() << " " << wi.z() << endl;
            // cout << "uv: " << uv.x() << " " << uv.y() << endl;

            // cout << "Le: " << endl
            //      << "wi: " << endl
            //      << wi << endl
            //      << "uv: " << endl
            //      << uv << endl
            //      << "val: " << endl
            //      << val << endl;

            return val * m_strength;
        }

        Color3f sampleLi_no_occlusion(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Point3f &pp, Normal3f &nn,
                                      Sampler *sampler, const Scene *scene) override {
            Point2f uv = Warp::squareToHierarchical(sampler->next2D());
            Vector3f loc_wi         = uv2sphere(uv).normalized();
            wi = m_toWorld * loc_wi;
            pp         = wi * Scene_Boundary;
            nn         = -wi;

            float cosTheta;
            if (ref_n.norm() == 0.f)
                cosTheta = 1.f;
            else
                cosTheta = ref_n.dot(wi);

            if (cosTheta < 0.0f)
                return 0.f;

            Color3f bm_val = lookup(m_bm, uv);
            float pdf      = pdfLi(ref, ref_n, loc_wi, nn);
            Color3f val    = bm_val * cosTheta / pdf;

            return isnan(val.x()) || isinf(val.x()) ? 0.f : val;
        }

        bool isInfinite() const override { return true; }

        bool isDeltaLight() const override {
            return true;
        }

        [[nodiscard]] std::string toString() const override {
            return tfm::format(
                    "IBL[\n"
                    "  exr_path = %s\n"
                    "]",
                    m_exr_path);
        }

    private:
        std::string m_exr_path;
        float m_strength = 1.f;
        Bitmap m_bm;
        Transform m_toWorld;
        Transform m_toLight;
        int m_rows, m_cols;
    };

    NORI_REGISTER_CLASS(IBLLight, "IBL");
NORI_NAMESPACE_END