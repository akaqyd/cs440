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

#include <nori/scene.h>
#include <nori/bitmap.h>
#include <nori/integrator.h>
#include <nori/sampler.h>
#include <nori/camera.h>
#include <nori/emitter.h>
#include <nori/medium.h>

NORI_NAMESPACE_BEGIN

    Scene::Scene(const PropertyList &) {
        m_accel = new Accel();
    }

    Scene::~Scene() {
        delete m_accel;
        delete m_sampler;
        delete m_camera;
        delete m_integrator;
    }

    void Scene::activate() {
        m_accel->build();
        sampleEmitterInit();

        if (!m_integrator)
            throw NoriException("No integrator was specified!");
        if (!m_camera)
            throw NoriException("No camera was specified!");

        if (!m_sampler) {
            /* Create a default (independent) sampler */
            m_sampler = static_cast<Sampler *>(
                    NoriObjectFactory::createInstance("independent", PropertyList()));
        }

        cout << endl;
        cout << "Configuration: " << toString() << endl;
        cout << endl;
    }

    void Scene::addChild(NoriObject *obj) {
        switch (obj->getClassType()) {
            case EMesh: {
                Mesh *mesh = static_cast<Mesh *>(obj);
                m_accel->addMesh(mesh);
                m_meshes.push_back(mesh);
                if (mesh->isEmitter())
                    m_emitters.push_back(mesh->getEmitter());
            }
                break;

            case EEmitter: {
                Emitter *em = static_cast<Emitter *>(obj);
                m_emitters.push_back(em);
                //Emitter *emitter = static_cast<Emitter *>(obj);
                /* TBD */
                // throw NoriException("Scene::addChild(): You need to implement this for emitters");
            }
                break;

            case ESampler:
                if (m_sampler)
                    throw NoriException("There can only be one sampler per scene!");
                m_sampler = static_cast<Sampler *>(obj);
                break;

            case ECamera:
                if (m_camera)
                    throw NoriException("There can only be one camera per scene!");
                m_camera = static_cast<Camera *>(obj);
                break;

            case EIntegrator:
                if (m_integrator)
                    throw NoriException("There can only be one integrator per scene!");
                m_integrator = static_cast<Integrator *>(obj);
                break;

            case EMedium: {
                Medium *m = static_cast<Medium *>(obj);
                if (m_medium != nullptr)
                    throw NoriException(
                            "Scene: tried to register multiple Medium instances!");
                m_medium = m;
            }
                break;

            case EDenoiser: {
                Denoiser *m = static_cast<Denoiser *>(obj);
                if (m_denoiser != nullptr)
                    throw NoriException(
                            "Scene: tried to register multiple Denoiser instances!");
                m_denoiser = m;
                if (m_camera && m_sampler)
                    m_denoiser->setConstants(m_camera->getOutputSize().y(),
                                             m_camera->getOutputSize().x(),
                                             m_sampler->getSampleCount()
                    );
            }
                break;

            default:
                throw NoriException("Scene::addChild(<%s>) is not supported!",
                                    classTypeName(obj->getClassType()));
        }
    }

    std::string Scene::toString() const {
        std::string meshes;
        for (size_t i = 0; i < m_meshes.size(); ++i) {
            meshes += std::string("  ") + indent(m_meshes[i]->toString(), 2);
            if (i + 1 < m_meshes.size())
                meshes += ",";
            meshes += "\n";
        }

        return tfm::format(
                "Scene[\n"
                "  integrator = %s,\n"
                "  sampler = %s\n"
                "  camera = %s,\n"
                "  meshes = {\n"
                "  %s  }\n"
                "]",
                indent(m_integrator->toString()),
                indent(m_sampler->toString()),
                indent(m_camera->toString()),
                indent(meshes, 2)
        );
    }

    void Scene::sampleEmitterInit() {
        for (int i = 0; i < (int) m_emitters.size(); i++) {
            emitter_dpdf.append(1);
        }
        emitter_dpdf.normalize();
    }

    Emitter *Scene::sampleEmitter(float x) const {
        size_t idx = emitter_dpdf.sample(x);
        return m_emitters[idx];
    }

    Color3f Scene::Tr(const Ray3f &_ray, Sampler *sampler) const {
        Ray3f ray(_ray);

        Point3f p1 = ray.o + ray.d * ray.maxt;

        Color3f Tr(1.f);
        while (true) {
            Intersection its;
            bool hit = this->rayIntersect(ray, its);
            if (hit)
                ray.maxt = (its.p - ray.o).norm();

            // Find non-opaque surface along the ray (should be rare!)
            if (hit && !its.mesh->getBSDF()->isNull())
                return Color3f(0.f);

            if (ray.m) Tr *= ray.m->transmittance(ray, sampler);

            if (!hit) break;
            its.spawnRayToPoint(ray, p1);
        }
        return Tr;

    }

    bool Scene::rayIntersectTr(const Ray3f &_ray, Sampler *sampler, Intersection &its, Color3f &Tr) const {
        Ray3f ray(_ray);
        Tr = 1.f;
        while (true) {
            bool hit = this->rayIntersect(ray, its);
            if (ray.m) Tr *= ray.m->transmittance(ray, sampler);

            if (!hit)
                return false;

            // Find non-opaque surface along the ray (should be rare!)
            if (!its.mesh->getBSDF()->isNull())
                return true;

            its.spawnRayToPoint(ray, ray.d);
        }

    }


    NORI_REGISTER_CLASS(Scene, "scene");
NORI_NAMESPACE_END
