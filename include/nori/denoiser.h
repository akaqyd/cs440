//
// Created by Qiyuan Dong on 12.05.22.
//

#pragma once

#include <nori/object.h>
#include <nori/common.h>
#include <string>

NORI_NAMESPACE_BEGIN
    class Denoiser : public NoriObject {
    private:


    public:
        int m_rows, m_cols;
        int m_spp;

        Denoiser() = default;

        virtual void run(Scene *scene, const std::string &filename, bool gui, int threadCount) = 0;

        virtual void setConstants(int rows, int cols, int spp) {
            m_rows = rows;
            m_cols = cols;
            m_spp = spp;
        }

        [[nodiscard]] EClassType getClassType() const override {
            return EDenoiser;
        }
    };

NORI_NAMESPACE_END