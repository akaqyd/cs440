//
// Created by Qiyuan Dong on 07.05.22.
//

#include <nori/frame.h>
#include <nori/warp.h>
#include <nori/texture.h>
#include <nori/bitmap.h>
#include <filesystem/resolver.h>

NORI_NAMESPACE_BEGIN

    class ColorTexture: public NoriTexture{
private:
    float m_scale = 1.f;
    Bitmap m_texture;
    std::string m_file_path;
    int m_rows, m_cols;

    public:
        ColorTexture(const PropertyList &propList) {
            m_scale = propList.getFloat("scale", 1);
            std::string path = propList.getString("texture_path");
            float strength = propList.getFloat("strength", 1.f);
            m_file_path = getFileResolver()->resolve(path).str();
            m_texture = Bitmap(m_file_path);
            m_texture *= strength;

            m_rows = m_texture.rows();
            m_cols = m_texture.cols();

            // normalize
            for (int r = 0; r < m_rows; r++) {
                for (int c = 0; c < m_cols; c++) {
                    m_texture(r, c).x() = std::fmax(m_texture(r,c).x(), 0);
                    m_texture(r, c).y() = std::fmax(m_texture(r,c).y(), 0);
                    m_texture(r, c).z() = std::fmax(m_texture(r,c).z(), 0);
                }
            }
        }

        Vector3f queryUV(const Point2f &uv) const override {
//            assert(uv[0] < 1 && uv[1] < 1);
//            assert(uv[0] > 0 && uv[1] > 0);

            float u = std::fmod(uv[0] + 1.f, 1.f);
            float v = std::fmod(uv[1] + 1.f, 1.f);

            /* Enable scale parameter */
            float r = (1.f - v) * (float)m_rows / m_scale;
            float c = u * (float)m_cols / m_scale;
            r = std::fmod(r, m_rows);
            c = std::fmod(c, m_cols);

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
            if (r0 >= 0 && r0 < (float)m_rows && c0 >= 0 && c0 < (float)m_cols) l00 = m_texture(r0, c0);
            if (r0 >= 0 && r0 < (float)m_rows && c1 >= 0 && c1 < (float)m_cols) l01 = m_texture(r0, c1);
            if (r1 >= 0 && r1 < (float)m_rows && c0 >= 0 && c0 < (float)m_cols) l10 = m_texture(r1, c0);
            if (r1 >= 0 && r1 < (float)m_rows && c1 >= 0 && c1 < (float)m_cols) l11 = m_texture(r1, c1);


            Color3f tmp_r0 = fmod(c - c0 - 0.5 + (float)m_cols, (float)m_cols) * l01 + fmod(c1 - c + 0.5 + (float)m_cols, (float)m_cols) * l00;
            Color3f tmp_r1 = fmod(c - c0 - 0.5 + (float)m_cols, (float)m_cols) * l11 + fmod(c1 - c + 0.5 + (float)m_cols, (float)m_cols) * l10;

            Color3f val = fmod(r - r0 - 0.5 + (float)m_rows, (float)m_rows) * tmp_r1 + fmod(r1 - r + 0.5 + (float)m_rows, (float)m_rows) * tmp_r0;

            return Vector3f(val.x(), val.y(), val.z());

//            val = val * 2.f - 1.f;
//            Normal3f n(val.x(), val.y(), val.z());
//            n.normalize();
//            return n;

             /* Old version, buggy, produces ugly seams */
//            float u = uv.x(), v = uv.y();
//            assert(u >= 0 && u <= 1 && v >= 0 && v <= 1);
//
//            float r = (1.f - v) * (float)m_rows / m_scale;
//            float c = u * (float)m_cols / m_scale;
//            r = std::fmod(r, m_rows);
//            c = std::fmod(c, m_cols);
//
//            int r0 = std::floor(r), c0 = std::floor(c);
//            int r1 = r0 + 1, c1 = c0 + 1;
//
//            Color3f v00(0), v01(0), v10(0), v11(0), v_top(0), v_bottom(0), val(0);
//            if (r0 >= 0 && c0 >= 0 && r0 < m_rows && c0 < m_cols)
//                v00 = m_texture(r0, c0);
//            if (r0 >= 0 && c1 >= 0 && r0 < m_rows && c1 < m_cols)
//                v01 = m_texture(r0, c1);
//            if (r1 >= 0 && c0 >= 0 && r1 < m_rows && c0 < m_cols)
//                v10 = m_texture(r1, c0);
//            if (r1 >= 0 && c1 >= 0 && r1 < m_rows && c1 < m_cols)
//                v11 = m_texture(r1, c1);
//
//            float a = c - (float)c0;
//            float b = r - (float)r0;
//            v_top = (1.f - a) * v00 + a * v01;
//            v_bottom = (1.f - a) * v10 + a * v11;
//            val = (1.f - b) * v_top + b * v_bottom;

//            val = (val / 255.f) * 2.f - 1.f;

//            val = val * 2.f - 1.f;
//
//            Normal3f n(val.x(), val.y(), val.z());
//            n.normalize();
//            return n;
        }

        bool isColorTexture() override {
            return true;
        }

        /// Return a human-readable summary
        std::string toString() const override {
            return tfm::format(
                    "ColorTexture[\n"
                    "  m_file_path= %s\n"
                    "  m_scale= %f\n"
                    "]", m_file_path, m_scale);
        }

        EClassType getClassType() const override { return ENoriTexture; }

    };

    NORI_REGISTER_CLASS(ColorTexture, "color_texture");
NORI_NAMESPACE_END
