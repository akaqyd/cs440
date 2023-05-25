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

#include <nori/bsdf.h>
#include <nori/frame.h>
#include <nori/warp.h>

NORI_NAMESPACE_BEGIN

/**
 * \brief Diffuse / Lambertian BRDF model
 */
    class NullDiffuse : public BSDF {
    public:
        NullDiffuse(const PropertyList &propList) {
        }

        /// Evaluate the BRDF model
        Color3f eval(const BSDFQueryRecord &bRec) const {
            throw NoriException("NullBSDF!");
        }

        /// Compute the density of \ref sample() wrt. solid angles
        float pdf(const BSDFQueryRecord &bRec) const {
            throw NoriException("NullBSDF!");
        }

        /// Draw a a sample from the BRDF model
        Color3f sample(BSDFQueryRecord &bRec, const Point2f &sample) const {
            throw NoriException("NullBSDF!");

        }

        bool isNull() const {
            return true;
        }

        bool isDiffuse() const {
            return false;
        }

        /// Return a human-readable summary
        std::string toString() const {
            return tfm::format(
                    "NullDiffuse[\n"
                    "]");
        }

        EClassType getClassType() const { return EBSDF; }

    private:
    };

    NORI_REGISTER_CLASS(NullDiffuse, "null");
NORI_NAMESPACE_END
