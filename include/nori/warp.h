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

#include <nori/common.h>
#include <nori/sampler.h>

#include <string>

NORI_NAMESPACE_BEGIN

struct MipMapNode {
    float l00, l01, l10, l11;
    float lsum;
    float xmin, xmax, ymin, ymax;
    MipMapNode *children[4];

    MipMapNode(float _l00, float _l01, float _l10, float _l11) {
        l00  = _l00;
        l01  = _l01;
        l10  = _l10;
        l11  = _l11;
        lsum = l00 + l01 + l10 + l11;
        for (int i = 0; i < 4; ++i)
            children[i] = nullptr;
    }

    void setRange(float _xmin, float _xmax, float _ymin, float _ymax) {
        xmin = _xmin;
        xmax = _xmax;
        ymin = _ymin;
        ymax = _ymax;
    }

    void setChildren(MipMapNode *c1, MipMapNode *c2, MipMapNode *c3, MipMapNode *c4) {
        children[0] = c1;
        children[1] = c2;
        children[2] = c3;
        children[3] = c4;
    }

    bool isLeaf() {
        return children[0] == nullptr;
    }

    /**
     * Warp a sample point according to the parameters stored in this node.
     * @param p Reference of the point to be warped in global coordinates.
     * @return The (x, y) local coordinates of the warped point.
     **/
    Point2f warpPoint(Point2f &p) {
        assert(p.x() < 1 && p.y() < 1);
        assert(p.x() > 0 && p.y() > 0);

        float x = p.x();
        float y = p.y();

        float xloc = (x - xmin) / (xmax - xmin);
        float yloc = (y - ymin) / (ymax - ymin);

        if (xloc <= (l00 + l10) / lsum) {
            xloc = xloc / ((l00 + l10) / lsum) * 0.5;
            if (yloc <= (l10 / (l10 + l00)))
                yloc = yloc / (l10 / (l10 + l00)) * 0.5;
            else
                yloc = 0.5 + (yloc - (l10 / (l10 + l00))) / (l00 / (l00 + l10)) * 0.5;
        } else {
            xloc = 0.5 + (xloc - ((l00 + l10) / lsum)) / ((l01 + l11) / lsum) * 0.5;
            if (yloc <= (l11 / (l11 + l01)))
                yloc = yloc / (l11 / (l11 + l01)) * 0.5;
            else
                yloc = 0.5 + (yloc - (l11 / (l11 + l01))) / (l01 / (l01 + l11)) * 0.5;
        }
        p(0) = xmin + xloc * (xmax - xmin);
        p(1) = ymin + yloc * (ymax - ymin);

        return Point2f(xloc, yloc);
    }

    /**
     *
     * @brief Return the pointer to the mipmap node in the next level
     *
     * @param loc (x, y) in the local coordinate after being warped in the current node
     * @return MipMapNode* Pointer the the next node.
     */
    MipMapNode *nextNode(Point2f loc) {
        if (isLeaf())
            return nullptr;

        float x = loc.x(), y = loc.y();
        if (x <= 0.5)
            if (y <= 0.5)
                return children[2];
            else
                return children[0];
        else if (y <= 0.5)
            return children[3];
        else
            return children[1];
    }
};

/// A collection of useful warping functions for importance sampling
class Warp {
   public:
    /// Dummy warping function: takes uniformly distributed points in a square and just returns them
    static Point2f squareToUniformSquare(const Point2f &sample);

    /// Probability density of \ref squareToUniformSquare()
    static float squareToUniformSquarePdf(const Point2f &p);

    /// Sample a 2D tent distribution
    static Point2f squareToTent(const Point2f &sample);

    /// Probability density of \ref squareToTent()
    static float squareToTentPdf(const Point2f &p);

    /// Uniformly sample a vector on a 2D disk with radius 1, centered around the origin
    static Point2f squareToUniformDisk(const Point2f &sample);

    /// Probability density of \ref squareToUniformDisk()
    static float squareToUniformDiskPdf(const Point2f &p);

    /// Uniformly sample a vector on the unit sphere with respect to solid angles
    static Vector3f squareToUniformSphere(const Point2f &sample);

    /// Probability density of \ref squareToUniformSphere()
    static float squareToUniformSpherePdf(const Vector3f &v);

    /// Uniformly sample a vector on the unit hemisphere around the pole (0,0,1) with respect to solid angles
    static Vector3f squareToUniformHemisphere(const Point2f &sample);

    /// Probability density of \ref squareToUniformHemisphere()
    static float squareToUniformHemispherePdf(const Vector3f &v);

    /// Uniformly sample a vector on the unit hemisphere around the pole (0,0,1) with respect to projected solid angles
    static Vector3f squareToCosineHemisphere(const Point2f &sample);

    /// Probability density of \ref squareToCosineHemisphere()
    static float squareToCosineHemispherePdf(const Vector3f &v);

    /// Warp a uniformly distributed square sample to a Beckmann distribution * cosine for the given 'alpha' parameter
    static Vector3f squareToBeckmann(const Point2f &sample, float alpha);

    /// Probability density of \ref squareToBeckmann()
    static float squareToBeckmannPdf(const Vector3f &m, float alpha);


    // Below is for hierarchical warping
    static Point2f squareToHierarchical(const Point2f &sample);
    static float squareToHierarchicalPdf(const Point2f &v);
    static void buildMipMap(std::string EXR_filename);
    static void setFilename(std::string fn) {EXR_filename = fn;}


    // Below is for GGX used in full microfacet BSDF
    /// Warp a uniformly distributed square sample to a GGX distribution * cosine for the given 'alpha' parameter
    static Vector3f squareToGGX(const Point2f &sample, float alpha);

    /// Probability density of \ref squareToGGX()
    static float squareToGGXPdf(const Vector3f &m, float alpha);

   private:
    static std::string EXR_filename;
    static MipMapNode *root;
    static MatrixXf luminance;
    static float luminance_sum;
};

NORI_NAMESPACE_END
