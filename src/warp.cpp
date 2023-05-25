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

#include <nori/bitmap.h>
#include <nori/frame.h>
#include <nori/vector.h>
#include <nori/warp.h>

#include <cassert>
#include <cmath>
#include <vector>

NORI_NAMESPACE_BEGIN

Point2f Warp::squareToUniformSquare(const Point2f &sample) {
    return sample;
}

float Warp::squareToUniformSquarePdf(const Point2f &sample) {
    return ((sample.array() >= 0).all() && (sample.array() <= 1).all()) ? 1.0f : 0.0f;
}

float tent(float t) {
    if (t <= 1 && t >= -1)
        return 1 - fabs(t);
    else
        return 0;
}

float tentCDFInv(float y) {
    if (0 <= y && y < 0.5)
        return sqrt(2 * y) - 1;
    else if (0.5 <= y && y <= 1)
        return 1 - sqrt(2 - 2 * y);
    else
        throw NoriException("Input is out of range [0, 1)!");
}

Point2f Warp::squareToTent(const Point2f &sample) {
    return Point2f(
        tentCDFInv(sample[0]), tentCDFInv(sample[1]));
}

float Warp::squareToTentPdf(const Point2f &p) {
    return tent(p[0]) * tent(p[1]);
}

Point2f Warp::squareToUniformDisk(const Point2f &sample) {
    float r     = sqrt(sample[0]);
    float theta = 2 * M_PI * sample[1];
    return Point2f(
        r * cos(theta), r * sin(theta));
}


float Warp::squareToUniformDiskPdf(const Point2f &p) {
    float r2 = p[0] * p[0] + p[1] * p[1];
    if (r2 <= 1 && r2 >= 0)
        return INV_PI;
    else
        return 0;
}

Vector3f Warp::squareToUniformSphere(const Point2f &sample) {
    float u = sample[0];
    float v = sample[1];

    float sinTheta = sqrt(-4 * u * u + 4 * u);
    float phi      = 2 * M_PI * v;

    return Vector3f(
        sinTheta * cos(phi),
        sinTheta * sin(phi),
        1 - 2 * u);
}

float Warp::squareToUniformSpherePdf(const Vector3f &v) {
    if (fabs(v.norm() - 1.0) < 1e-6)
        return INV_PI / 4.0;
    else
        return 0;
}

Vector3f Warp::squareToUniformHemisphere(const Point2f &sample) {
    float u = sample[0];
    float v = sample[1];

    float sinTheta = sqrt(-1 * u * u + 2 * u);
    float phi      = 2 * M_PI * v;

    return Vector3f(
        sinTheta * cos(phi),
        sinTheta * sin(phi),
        1 - u);
}

float Warp::squareToUniformHemispherePdf(const Vector3f &v) {
    if (fabs(v.norm() - 1.0) < 1e-6 && v[2] >= 0.0)
        return INV_PI / 2.0;
    else
        return 0;
}

Vector3f Warp::squareToCosineHemisphere(const Point2f &sample) {
    // float u = sample[0];
    // float v = sample[1];

    // float theta = acos(1 - 2 * u) / 2.0;
    // float phi      = 2 * M_PI * v;

    // return Vector3f(
    //     sin(theta) * cos(phi),
    //     sin(theta) * sin(phi),
    //     cos(theta));

    // Malley’s method from PBRT pp.779
    Point2f d = squareToUniformDisk(sample);
    float z   = sqrt(fmax((float)0, 1 - d[0] * d[0] - d[1] * d[1]));
    return Vector3f(d[0], d[1], z);
}

float Warp::squareToCosineHemispherePdf(const Vector3f &v) {
    if (v[2] < 0.0 || fabs(v.norm() - 1.0f) > 1e-6)
        return 0;

    float cosTheta = v[2];
    // float theta = acos(cosTheta);
    // return INV_PI * sin(theta) * cosTheta;
    return INV_PI * cosTheta;
}

Vector3f Warp::squareToBeckmann(const Point2f &sample, float alpha) {
    // CDF(w) = 1 - exp(-tan(w)^2 / alpha^2)
    float phi    = 2.0f * M_PI * sample[0];
    float sinPhi = sin(phi), cosPhi = cos(phi);
    float tan2Theta = -1.0f * alpha * alpha * log(sample[1]);
    float cosTheta  = sqrt(1.0f / (1.0f + tan2Theta));
    float sinTheta  = sqrt(tan2Theta / (1.0f + tan2Theta));
    return Vector3f(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
}

float Warp::squareToBeckmannPdf(const Vector3f &m, float alpha) {
    if (m.z() <= 0)
        return 0.0f;
    float cosT = Frame::cosTheta(m);
    float tanT = Frame::tanTheta(m);
    return INV_PI * exp(-1 * tanT * tanT / (alpha * alpha)) / (alpha * alpha * cosT * cosT * cosT);
}

//----------------------------------------------------------------
// Below is for Hierarchical Warping
//----------------------------------------------------------------
std::string Warp::EXR_filename = "scene/pa3/4by4.exr";
MipMapNode *Warp::root         = nullptr;
float Warp::luminance_sum      = 0.0f;
MatrixXf Warp::luminance       = MatrixXf(0, 0);

void Warp::buildMipMap(std::string EXR_filename) {
    Bitmap bm(EXR_filename);
    luminance.resize(bm.rows(), bm.cols());

    std::vector<std::vector<MipMapNode *> > prev(bm.rows(), std::vector<MipMapNode *>(bm.cols()));

    float dr = 1.0f / bm.rows(), dc = 1.0f / bm.cols();
    for (int r = 0; r < bm.rows(); ++r) {
        for (int c = 0; c < bm.cols(); ++c) {
            float l         = bm(r, c).getLuminance();
            luminance(r, c) = l;
            luminance_sum += l;
            prev[r][c]       = new MipMapNode(0, 0, 0, 0);
            prev[r][c]->lsum = l;
            prev[r][c]->setRange(c * dc, (c + 1.f) * dc, (bm.rows() - 1.f - r) * dr, (bm.rows() - r) * dr);
        }
    }

    int rows = bm.rows() / 2, cols = bm.cols() / 2;
    while (rows > 0 && cols > 0) {
        std::vector<std::vector<MipMapNode *> > curr(rows, std::vector<MipMapNode *>(cols));
        for (int r = 0; r < rows; ++r) {
            for (int c = 0; c < cols; ++c) {
                int rr = r * 2, cc = c * 2;
                curr[r][c] = new MipMapNode(prev[rr][cc]->lsum, prev[rr][cc + 1]->lsum, prev[rr + 1][cc]->lsum, prev[rr + 1][cc + 1]->lsum);
                curr[r][c]->setRange(prev[rr][cc]->xmin, prev[rr][cc + 1]->xmax, prev[rr + 1][cc]->ymin, prev[rr][cc]->ymax);
                if (rows != bm.rows() / 2)
                    curr[r][c]->setChildren(prev[rr][cc], prev[rr][cc + 1], prev[rr + 1][cc], prev[rr + 1][cc + 1]);
            }
        }

        rows /= 2;
        cols /= 2;
        prev = curr;
    }
    root = prev[0][0];
};

Point2f Warp::squareToHierarchical(const Point2f &sample) {
    if (root == nullptr) 
        throw NoriException("Must call Warp::buildMipMap() first.");

    MipMapNode *currNode = root;
    Point2f p(sample);
    while (currNode != nullptr) {
        Point2f loc = currNode->warpPoint(p);
        currNode = currNode->nextNode(loc);
    }
    return p;
};

float Warp::squareToHierarchicalPdf(const Point2f &p) {
    if (root == nullptr) 
        throw NoriException("Must call Warp::buildMipMap() first.");

    if ((p.array() > 0).all() && (p.array() < 1).all()) {
        int r = luminance.rows() * (1 - p(1));
        int c = luminance.cols() * p(0);
        r     = fmin(luminance.rows(), r);
        r     = fmax(0, r);
        c     = fmin(luminance.cols(), c);
        c     = fmax(0, c);

        return luminance(r, c) * luminance.rows() * luminance.cols() / luminance_sum;
    } else {
        return 0;
    }
};

Vector3f Warp::squareToGGX(const Point2f &sample, float alpha) {
    float phi    = 2.0f * M_PI * sample[0];
    float sinPhi = sin(phi), cosPhi = cos(phi);

    float tan2Theta = alpha * alpha * sample[1] / (1 - sample[1]);
    float cosTheta  = sqrt(1.0f / (1.0f + tan2Theta));
    float sinTheta  = sqrt(tan2Theta / (1.0f + tan2Theta));

    return Vector3f(sinTheta * cosPhi, sinTheta * sinPhi, cosTheta);
}

float Warp::squareToGGXPdf(const Vector3f &m, float alpha) {
    if (m.z() <= 0)
        return 0.f;
    
    float cosT = Frame::cosTheta(m);
    float tanT = Frame::tanTheta(m);

    float pdf = alpha * alpha;
    // pdf /= M_PI * powf(cosT, 4) * powf(alpha * alpha + tanT * tanT, 2);
    pdf /= M_PI * powf(cosT, 3) * powf(alpha * alpha + tanT * tanT, 2);

    return pdf;
}



NORI_NAMESPACE_END
