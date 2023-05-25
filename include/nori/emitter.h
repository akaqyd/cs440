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

#pragma once

#include <nori/object.h>
#include <nori/shape.h>
#include <pcg32.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Superclass of all emitters
 */
class Emitter : public NoriObject {
   public:
    Shape *getShape() const { return m_shape; }

    void setShape(Shape *s) { m_shape = s; }

    /**
     * @brief Samples an incident direction at a point in the scene along which radiance from the light may be arriving. Note that the cosine factor in the reflection equation, the pdf used to sample this incident direction and occlusion are already taken into account.
     *
     * @param ref The intersection point in the scene at which the incident direction is sampled
     * @param ref_n The corresponding normal of p
     * @param wi To store the sampled incident direction along which radiance from the light may be arriving. Normalized.
     * @param p To store sampled point on the light source
     * @param n To store  normal of p
     * @param sampler Used to sample.
     * @param scene Pointer to the scene, which will be used to test occlusion
     * @return Color3f The illumination from the light source along the sampled direction.
     */
    virtual Color3f sampleLi(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Point3f &p, Normal3f &n, Sampler *sampler, const Scene *scene) const = 0;

    virtual Color3f sampleLi_no_occlusion(const Point3f &ref, const Normal3f &ref_n, Vector3f &wi, Point3f &p, Normal3f &n, Sampler *sampler, const Scene *scene) = 0;

    /**
     * @brief Returns the radiance from a point on the surface of the light source along a given direction.
     *
     * @param p The point on the surface of the light source
     * @param n The corresponding normal of p
     * @param wo The outgoing direction from p to be evaluated
     * @return Color3f Radiance emitted by the light source at point p to the direction wo.
     */
    virtual Color3f le(const Point3f &p, const Normal3f &n, const Vector3f &wo) const = 0;

    /**
     * UV support
     */
    virtual Color3f le(const Point3f &p, const Normal3f &n, const Vector3f &wo, const Color3f &texture) const {
        return le(p, n, wo);
    }

    /**
     * @brief Used for infinite light source when the ray hits the scene boundary.
     *
     * @param ray 
     * @return Color3f Radiance emitted by the infinite light source for the direction ray.d
     */
    virtual Color3f le(const Ray3f &ray) const {
        return 0.f;
    }

    virtual bool isInfinite() const { return false;}

    /**
     * @brief Returns the probability density with respect to ** solid angle ** for the light’s Sample_Li() method to sample the point p on the surface of the light source
     *
     * @param ref The intersection point in the scene at which the incident direction is sampled
     * @param ref_n The corresponding normal of ref
     * @param p The point on light source
     * @param n Normal of p
     *
     * @return float The probability density with respect to solid angle for the light’s Sample_Li() method to sample the direction wi
     */
    virtual float pdfLi(const Point3f &ref, const Normal3f &ref_n, const Point3f &p, const Normal3f &n) const = 0;

    virtual bool isDeltaLight() const = 0;

    /**
     * \brief Return the type of object (i.e. Mesh/Emitter/etc.)
     * provided by this instance
     * */
    EClassType getClassType() const { return EEmitter; }

   protected:
    Color3f m_radiance;
    Shape *m_shape;
};

NORI_NAMESPACE_END
