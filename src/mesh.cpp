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

#include <nori/bbox.h>
#include <nori/bsdf.h>
#include <nori/emitter.h>
#include <nori/mesh.h>
#include <nori/warp.h>
#include <nori/shape.h>
#include <nori/texture.h>
#include <nori/common.h>

#include <Eigen/Geometry>

NORI_NAMESPACE_BEGIN

    Mesh::Mesh() {}

    Mesh::~Mesh() {
        delete m_bsdf;
        delete m_emitter;
    }

    void Mesh::activate() {
        if (!m_bsdf) {
            /* If no material was assigned, instantiate a null BRDF */
            m_bsdf = static_cast<BSDF *>(
                    NoriObjectFactory::createInstance("null", PropertyList()));

            // It is ok to be nullptr for the boundary of medium
//            m_bsdf = nullptr;
//            std::cout << "Mesh " << this->m_name
//                      << " has no BSDF. This is only allowed for meshes used as medium boundary." << std::endl;
        }

        // preprocess the mesh if used as area light source
        if (m_emitter) {
            uint32_t triCnt = getTriangleCount();
            dpdf.reserve(triCnt);
            for (uint32_t i = 0; i < triCnt; ++i)
                dpdf.append(surfaceArea(i));
            dpdf.normalize();
            area_sum = dpdf.getSum();
        }

    }

    float Mesh::samplePoint(Sampler *sampler, Point3f &p, Normal3f &n) const {
        uint32_t idx = dpdf.sample(sampler->next1D());
        Point2f sample2 = sampler->next2D();
        float alpha = 1 - sqrt(1 - sample2[0]);
        float beta = sample2[1] * sqrt(1 - sample2[0]);

        uint32_t v0 = m_F(0, idx), v1 = m_F(1, idx), v2 = m_F(2, idx);
        Point3f p0 = m_V.col(v0), p1 = m_V.col(v1), p2 = m_V.col(v2);
        p = alpha * p0 + beta * p1 + (1.0 - alpha - beta) * p2;

        if (m_N.size() > 0) {
            n = alpha * m_N.col(v0) + beta * m_N.col(v1) + (1.0 - alpha - beta) * m_N.col(v2);
        } else {
            n = (p1 - p0).cross(p2 - p0).normalized();
        }

        return 1 / area_sum;
    }

    float Mesh::samplePoint(Sampler *sampler, Point3f &p, Normal3f &n, Point2f &uv) const {
        uint32_t idx = dpdf.sample(sampler->next1D());
        Point2f sample2 = sampler->next2D();
        float alpha = 1 - sqrt(1 - sample2[0]);
        float beta = sample2[1] * sqrt(1 - sample2[0]);

        uint32_t v0 = m_F(0, idx), v1 = m_F(1, idx), v2 = m_F(2, idx);
        Point3f p0 = m_V.col(v0), p1 = m_V.col(v1), p2 = m_V.col(v2);
        p = alpha * p0 + beta * p1 + (1.0 - alpha - beta) * p2;

        if (m_N.size() > 0) {
            n = alpha * m_N.col(v0) + beta * m_N.col(v1) + (1.0 - alpha - beta) * m_N.col(v2);
        } else {
            n = (p1 - p0).cross(p2 - p0).normalized();
        }

        if (m_UV.size() > 0) {
            uv = alpha * m_UV.col(v0) + beta * m_UV.col(v1) + (1.0 - alpha - beta) * m_UV.col(v2);
        } else {
            uv = Point2f(-1.f);
        }

        return 1 / area_sum;
    }

    float Mesh::surfaceArea(uint32_t index) const {
        uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);

        const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

        return 0.5f * Vector3f((p1 - p0).cross(p2 - p0)).norm();
    }

    bool Mesh::rayIntersect(uint32_t index, const Ray3f &ray, float &u, float &v, float &t) const {
        uint32_t i0 = m_F(0, index), i1 = m_F(1, index), i2 = m_F(2, index);
        const Point3f p0 = m_V.col(i0), p1 = m_V.col(i1), p2 = m_V.col(i2);

        /* Find vectors for two edges sharing v[0] */
        Vector3f edge1 = p1 - p0, edge2 = p2 - p0;

        /* Begin calculating determinant - also used to calculate U parameter */
        Vector3f pvec = ray.d.cross(edge2);

        /* If determinant is near zero, ray lies in plane of triangle */
        float det = edge1.dot(pvec);

        if (det > -1e-8f && det < 1e-8f)
            return false;
        float inv_det = 1.0f / det;

        /* Calculate distance from v[0] to ray origin */
        Vector3f tvec = ray.o - p0;

        /* Calculate U parameter and test bounds */
        u = tvec.dot(pvec) * inv_det;
        if (u < 0.0 || u > 1.0)
            return false;

        /* Prepare to test V parameter */
        Vector3f qvec = tvec.cross(edge1);

        /* Calculate V parameter and test bounds */
        v = ray.d.dot(qvec) * inv_det;
        if (v < 0.0 || u + v > 1.0)
            return false;

        /* Ray intersects triangle -> compute t */
        t = edge2.dot(qvec) * inv_det;

        return t >= ray.mint && t <= ray.maxt;
    }

    BoundingBox3f Mesh::getBoundingBox(uint32_t index) const {
        BoundingBox3f result(m_V.col(m_F(0, index)));
        result.expandBy(m_V.col(m_F(1, index)));
        result.expandBy(m_V.col(m_F(2, index)));
        return result;
    }

    Point3f Mesh::getCentroid(uint32_t index) const {
        return (1.0f / 3.0f) *
               (m_V.col(m_F(0, index)) +
                m_V.col(m_F(1, index)) +
                m_V.col(m_F(2, index)));
    }

    void Mesh::addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case EBSDF:
                if (m_bsdf)
                    throw NoriException(
                            "Mesh: tried to register multiple BSDF instances!");
                m_bsdf = static_cast<BSDF *>(obj);
//                if (m_bsdf->isNull())
//                    m_bsdf = nullptr;
                break;

            case EEmitter: {
                Emitter *emitter = static_cast<Emitter *>(obj);
                if (m_emitter)
                    throw NoriException(
                            "Mesh: tried to register multiple Emitter instances!");
                m_emitter = emitter;
                m_emitter->setShape(this);
            }
                break;

            case EMedium: {
                Medium *m = static_cast<Medium *>(obj);
                if (m_mif.inside == nullptr && m_mif.outside == nullptr)
                    m_mif = MediumInterface(m, nullptr);
                else if (m_mif.inside != nullptr && m_mif.outside == nullptr)
                    m_mif.outside = m;
                else
                    throw NoriException(
                            "Mesh: tried to register multiple Medium instances!");
            }
                break;

            case ENoriTexture: {
                NoriTexture *m = static_cast<NoriTexture *>(obj);
                if (m->isNormalTexture()) {
                    if (m_normal_texture != nullptr)
                        throw NoriException(
                                "Mesh: tried to register multiple normal texture instances!");
                    m_normal_texture = m;
                } else if (m->isColorTexture()) {
                    if (m_color_texture != nullptr)
                        throw NoriException(
                                "Mesh: tried to register multiple color texture instances!");
                    m_color_texture = m;
                } else if (m->isRoughnessTexture()) {
                    if (m_roughness_texture != nullptr)
                        throw NoriException(
                                "Mesh: tried to register multiple roughness texture instances!");
                    m_roughness_texture = m;
                } else if (m->isEmitterTexture()) {
                    if (m_emitter_texture != nullptr)
                        throw NoriException(
                                "Mesh: tried to register multiple emitter texture instances!");
                    m_emitter_texture = m;
                } else {
                    throw NoriException(
                            "Mesh: tried to register an unknown texture instance!");
                }

            }
                break;

            default:
                throw NoriException("Mesh::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }

    std::string Mesh::toString() const {
        return tfm::format(
                "Mesh[\n"
                "  name = \"%s\",\n"
                "  vertexCount = %i,\n"
                "  triangleCount = %i,\n"
                "  bsdf = %s,\n"
                "  emitter = %s\n"
                "  medium_inside = %s\n"
                "  medium_outside = %s\n"
                "]",
                m_name,
                m_V.cols(),
                m_F.cols(),
                m_bsdf ? indent(m_bsdf->toString()) : std::string("null"),
                m_emitter ? indent(m_emitter->toString()) : std::string("null"),
                m_mif.inside ? indent(m_mif.inside->toString()) : std::string("null"),
                m_mif.outside ? indent(m_mif.outside->toString()) : std::string("null")
        );
    }

    std::string Intersection::toString() const {
        if (!mesh)
            return "Intersection[invalid]";

        return tfm::format(
                "Intersection[\n"
                "  p = %s,\n"
                "  t = %f,\n"
                "  uv = %s,\n"
                "  shFrame = %s,\n"
                "  geoFrame = %s,\n"
                "  mesh = %s\n"
                "]",
                p.toString(),
                t,
                uv.toString(),
                indent(shFrame.toString()),
                indent(geoFrame.toString()),
                mesh ? mesh->toString() : std::string("null"));
    }

NORI_NAMESPACE_END
