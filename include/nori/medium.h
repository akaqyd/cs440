//
// Created by Qiyuan Dong on 02.05.22.
//

#pragma once

#include <nori/object.h>
#include <nori/sampler.h>
#include <nori/color.h>
#include <nori/frame.h>
#include <nori/common.h>
#include <nori/bbox.h>
#include <nori/transform.h>
#include <filesystem/resolver.h>

#include <utility>
#include <openvdb/openvdb.h>
#include <openvdb/tools/Interpolation.h>
#include <openvdb/tools/GridTransformer.h>

NORI_NAMESPACE_BEGIN
    struct MediumInterface {
        const Medium *inside = nullptr, *outside = nullptr;

        // Constructors
        MediumInterface() = default;

        explicit MediumInterface(const Medium *medium) : inside(medium), outside(medium) {}

        MediumInterface(const Medium *inside, const Medium *outside)
                : inside(inside), outside(outside) {}

        [[nodiscard]] bool IsMediumTransition() const { return inside != outside; }
    };

    struct MediumInteraction {
        Point3f p;
        Normal3f n;
        Vector3f wi;
        MediumInterface mif;
        bool initialized = false;

        // Constructors
        MediumInteraction() = default;;

        MediumInteraction(Point3f p, Vector3f wi, const Medium *m) :
                p(std::move(p)), wi(std::move(wi)), mif(m, m), initialized(true) {}

        [[nodiscard]] bool isInitialized() const { return initialized; }

        void spawnRayToDir(Ray3f &ray, const Vector3f &d) const {
            ray.o = p;
            ray.d = d.normalized();
            ray.mint = Epsilon;
            ray.maxt = Scene_Boundary;
            ray.m = mif.inside;
            ray.update();
        }

        void spawnRayToPoint(Ray3f &ray, const Point3f &p1) const {
            Vector3f dir = p1 - p;
            ray.o = p;
            ray.d = dir.normalized();
            ray.mint = Epsilon;
            ray.maxt = dir.norm();
            ray.m = mif.inside;
            ray.update();
        }


    };

    class VDB {
    private:
        std::string file_path;
        float inv_max_density = -1.f;
        float max_density = -1.f;
        int nx, ny, nz;
        openvdb::FloatGrid::Ptr m_grid;
        openvdb::FloatGrid::Ptr m_grid_tmp;
        openvdb::CoordBBox m_bbox;
        int x_min, y_min, z_min;
        int x_max, y_max, z_max;

    public:

        VDB() = default;

        void reset(std::string &fp) {
            openvdb::initialize();
            file_path = getFileResolver()->resolve(fp).str();

            openvdb::initialize();
            openvdb::io::File file(file_path);
            file.open();

            openvdb::GridBase::Ptr baseGrid;
            openvdb::Coord dim;

            std::cout << "================================" << std::endl;

            // Iterate over all named grids
            for (openvdb::io::File::NameIterator nameIter = file.beginName();
                 nameIter != file.endName(); ++nameIter) {                // Read in only the grid we are interested in.

                if (nameIter.gridName() == "density") {
                    baseGrid = file.readGrid(nameIter.gridName());
                    std::cout << "Find grid: " << nameIter.gridName() << std::endl;
                } else
                    std::cout << "Skip grid: " << nameIter.gridName() << std::endl;
            }

            m_grid = openvdb::gridPtrCast<openvdb::FloatGrid>(baseGrid);


            std::cout << "Target grid has ";
            if (m_grid->hasUniformVoxels())
                std::cout << "uniform voxels,";
            else
                std::cout << "non-uniform voxels,";

            dim = m_grid->evalActiveVoxelDim();
            std::cout << " dimensions: " << dim << std::endl;
            nx = dim.x();
            ny = dim.y();
            nz = dim.z();

            std::cout << "Translation: " << m_grid->transform() << std::endl;

            // compute maximum density
            for (openvdb::FloatGrid::ValueOnCIter iter = m_grid->cbeginValueOn(); iter; ++iter)
                max_density = std::max(max_density, *iter);
            inv_max_density = 1.f / max_density;
            std::cout << "Maximum density: " << max_density << std::endl;

            // bounding box
            m_bbox = m_grid->evalActiveVoxelBoundingBox();
            x_min = m_bbox.getStart().x();
            y_min = m_bbox.getStart().y();
            z_min = m_bbox.getStart().z();
            x_max = m_bbox.getEnd().x();
            y_max = m_bbox.getEnd().y();
            z_max = m_bbox.getEnd().z();
            std::cout << "Bounding box: " << m_bbox << std::endl;

            std::cout << "================================" << std::endl;

        }

//        void transform(const Transform &mediumToWorld) {
//            m_grid_tmp = openvdb::FloatGrid::create(0);
//
//            Eigen::Matrix4f m = mediumToWorld.getMatrix();
//            openvdb::Mat4R xform(m(0, 0), m(0, 1), m(0, 2), m(0, 3),
//                                 m(1, 0), m(1, 1), m(1, 2), m(1, 3),
//                                 m(2, 0), m(2, 1), m(2, 2), m(2, 3),
//                                 m(3, 0), m(3, 1), m(3, 2), m(3, 3)
//            );
//
//            openvdb::tools::GridTransformer transformer(xform);
//            transformer.setThreaded(false);
//
//            // Resample using trilinear interpolation.
//            transformer.transformGrid<openvdb::tools::BoxSampler, openvdb::FloatGrid>(
//                    *m_grid, *m_grid_tmp);
//
//            m_grid = m_grid_tmp;
//
//            for (openvdb::FloatGrid::ValueOnCIter iter = m_grid->cbeginValueOn(); iter; ++iter)
//                max_density = std::max(max_density, *iter);
//            inv_max_density = 1.f / max_density;
//            std::cout << "Maximum density after transformation: " << max_density << std::endl;
//
//
//            // update bbox
////            m_bbox = m_grid->evalActiveVoxelBoundingBox();
////            x_min = m_bbox.getStart().x();
////            y_min = m_bbox.getStart().y();
////            z_min = m_bbox.getStart().z();
////            x_max = m_bbox.getEnd().x();
////            y_max = m_bbox.getEnd().y();
////            z_max = m_bbox.getEnd().z();
//        }

        float get_inv_max_density() const {
            return inv_max_density;
        }

        BoundingBox3f get_bbox() const {
//            return BoundingBox3f(Point3f(x_min, y_min, z_min),
//                                 Point3f(x_max, y_max, z_max));
            const openvdb::Vec3R min(x_min, y_min, z_min);
            const openvdb::Vec3R max(x_max, y_max, z_max);
            const openvdb::Vec3R min_world = m_grid->indexToWorld(min);
            const openvdb::Vec3R max_world = m_grid->indexToWorld(max);
            return BoundingBox3f(Point3f(min_world.x(), min_world.y(), min_world.z()),
                                 Point3f(max_world.x(), max_world.y(), max_world.z()));
        }

        /**
         * Query the density using the VDB grid for a given point
         * @param p Point at which the density is to be queried.
         *          Need to be transformed into medium local coordinates.
         *          NOT normalized.
         * @return The interpolated density.
         */
        float queryDensity(const Point3f &p) const {
//            const openvdb::Vec3R ijk(p.x() * (float)nx - .5f, p.y() * (float)ny - .5f, p.z() * (float)nz - .5f);
            const openvdb::Vec3R ijk(p.x(), p.y(), p.z());
            const openvdb::Vec3R ijk_index = m_grid->worldToIndex(ijk);
            float v = openvdb::tools::BoxSampler::sample(m_grid->tree(), ijk_index);

            return v;
        }

    };

    class Medium : public NoriObject {
    public:
        /**
         * Compute the beam transmittance along a given ray
         * @param ray
         * @param sampler
         * @return Transmittance on the interval between the ray.o and the point at a distance of Ray::tMax
         */
        [[nodiscard]] virtual Color3f transmittance(const Ray3f &ray, Sampler *sampler) const = 0;

        /**
         * Given a ray (p, ω), sample a medium scattering interaction along it.
         * The input ray will generally have **already been intersected** against the scene geometry.
         * Thus, implementations of this method shouldn’t ever sample a medium interaction at a point on the ray beyond ray.t_max.
         * Here, the algorithm neglects the effect of medium emission and assumes directionally constant medium properties
         *
         * Two cases can occur:
         *      - If Sample() doesn’t sample an interaction on the given ray interval [0, t_max], then the surface-related term Tr(p0→ p)Lo(p0, −ω) should be estimated
         *      - If it does sample an interaction, the integral term along the interval is to be estimated, and the provided MediumInteraction should be initialized accordingly.
         *
         * @param ray The ray along which the scattering interaction is sampled. Assumes that there is always a surface at ray.t_max. Particularly, p0= ray.o + t_max*ray.d is the point on the surface. Need to make sure that ray.d is *** normalized ***
         * @param sampler
         * @param mi MediumInteraction to be initialized if it does sample an interaction
         * @return Function value divided by the PDF at the sampled position
         */
        virtual Color3f sample(const Ray3f &ray, Sampler *sampler, MediumInteraction *mi) const = 0;

        /**
         * Given an incident direction, sample the phase function to get a exiting direction.
         * @param wi Incident normalized direction (point outward)
         * @param wo Exiting normalized direction (point outward)
         * @param sampler
         * @return pdf of this sample
         */
        virtual float sample_pf(const Vector3f &wi, Vector3f &wo, Sampler *sampler) const = 0;

        /**
         * Evaluate phase function value for a given pair of w_i and w_o.
         * @param wi Incident normalized direction (point outward)
         * @param wo Exiting normalized direction (point outward)
         * @return phase function value
         */
        [[nodiscard]] virtual float eval_pf(const Vector3f &wi, const Vector3f &wo) const = 0;


        [[nodiscard]] EClassType getClassType() const override {
            return EMedium;
        }

    };

    class PhaseFunction : public NoriObject {
    public:
        /**
         * Evaluate phase function value for a given pair of w_i and w_o.
         * @param wi Incident normalized direction (point outward)
         * @param wo Exiting normalized direction (point outward)
         * @return phase function value
         */
        [[nodiscard]] virtual float eval(const Vector3f &wi, const Vector3f &wo) const = 0;

        /**
         * Given an incident direction, sample a exiting direction.
         * @param wi Incident normalized direction (point outward)
         * @param wo Exiting normalized direction (point outward)
         * @param sampler
         * @return pdf of this sample
         */
        virtual float sample(const Vector3f &wi, Vector3f &wo, Sampler *sampler) const = 0;

        [[nodiscard]] EClassType getClassType() const override {
            return EPhaseFunction;
        }

    };

    class HGPhaseFunc : public PhaseFunction {
    private:
        const float g;
    public:
        explicit HGPhaseFunc(float g) : g(g) {}

        [[nodiscard]] float HG(float cosTheta) const {
            float d = 1 + g * g + 2 * g * cosTheta;
            return INV_PI / 4.0 * (1 - g * g) / (d * std::sqrt(d));
        }

        [[nodiscard]] float eval(const Vector3f &wi, const Vector3f &wo) const override {
            return HG(wi.dot(wo));
        }

        float sample(const Vector3f &wi, Vector3f &wo, Sampler *sampler) const override {
            Point2f sample = sampler->next2D();
            float cosTheta, sinTheta, phi;
            if (std::abs(g) < 0.001)
                cosTheta = 1.f - 3.f * sample[0];
            else {
                float t = (1 - g * g) / (1 + g - 2 * g * sample[0]);
                cosTheta = -(1 + g * g - t * t) / (2 * g);
            }
            sinTheta = std::sqrt(std::max(0.f, 1 - cosTheta * cosTheta));
            phi = 2 * M_PI * sample[1];

            Frame loc(wi);
            wo = sinTheta * std::cos(phi) * loc.s + sinTheta * std::sin(phi) * loc.t +
                 cosTheta * loc.n;

            return HG(cosTheta);
        }

        [[nodiscard]] std::string toString() const override {
            return tfm::format(
                    "[ HenyeyGreenstein g: %f ]", g);
        }

    };


NORI_NAMESPACE_END