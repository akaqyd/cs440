#pragma once

#include <nori/common.h>
#include <nori/object.h>

NORI_NAMESPACE_BEGIN
class Shape : public NoriObject {
public:

    /**
     * @brief Uniformly samples a point on the surface of the shape and returns the local geometric information about the sampled point
     *
     * @param sampler Random variable used to sample a triangle and a point inside the triangle
     * @param p To store the sampled point
     * @param n To store the corresponding normal of p
     * @return float The probability density for this sample
     */
    virtual float samplePoint(Sampler* sampler, Point3f &p, Normal3f &n) const = 0;


    /**
     * UV support
     */
    virtual float samplePoint(Sampler* sampler, Point3f &p, Normal3f &n, Point2f &uv) const {
        uv = Point2f(-1.f);
        return samplePoint(sampler, p, n);
    }

    /**
     * @brief Samples a point on the surface but also takes the point in the scene from which the surface of the shape is being integrated over as a parameter. Try to only sample the portion of the shape that is potentially visible from that point. Fall back to original sampling method if this is not possible.
     *
     * @param ref The point in the scene from which the surface of the shape is being integrated over as a parameter
     * @param ref_n Corresponding normal of ref
     * @param sampler Random variable used to sample a point
     * @param p To store the sampled point
     * @param n To store the corresponding normal of p
     * @param over_area To store if the probability density returned is defined over area (over solid angel otherwise)
     * @return float The probability density for this sample
     */
    virtual float samplePoint(const Point3f &ref, const Normal3f &ref_n, Sampler* sampler, Point3f &p, Normal3f &n, bool &over_area) const = 0;


    /**
     * @brief Returns the probability density corresponding to the first sampling approach, which uniformly sample a point on the surface
     *
     * @return float Probability density defined over AREA
     */
    virtual float pdf() const = 0;

    /**
     * @brief Returns the probability density to sample the direction wi with respect to solid angle from the reference point ref
     *
     * @param ref The point in the scene from which the surface of the shape is being integrated over as a parameter
     * @param ref_n Corresponding normal of ref
     * @param wi The incident direction at ref
     * @param over_area To store if the probability density returned is defined over area (over solid angel otherwise)
     * @return float the probability density to sample the direction wi with respect to SOLID ANGLE from the reference point ref
     */
    virtual float pdf(const Point3f &ref, const Normal3f &ref_n, const Vector3f &wi, bool &over_area) const = 0;

    virtual float getAreaSum() const = 0;

    virtual bool rayIntersect(const Ray3f &ray) const = 0;

    EClassType getClassType() const { return EEmitter; }

};

NORI_NAMESPACE_END
