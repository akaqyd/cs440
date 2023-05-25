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

#include <filesystem/resolver.h>
#include <nori/bitmap.h>
#include <nori/block.h>
#include <nori/camera.h>
#include <nori/gui.h>
#include <nori/integrator.h>
#include <nori/parser.h>
#include <nori/sampler.h>
#include <nori/scene.h>
#include <nori/timer.h>
#include <tbb/blocked_range.h>
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/concurrent_vector.h>
#include <openvdb/openvdb.h>
#include <openvdb/tools/LevelSetSphere.h>

#include <thread>

using namespace nori;

static int threadCount = -1;
static bool gui = true;
static tbb::concurrent_vector<double> variances;    // store variances of each pixel
static bool keepVariances = false;


static void renderBlock(const Scene *scene, Sampler *sampler, ImageBlock &block) {
    const Camera *camera = scene->getCamera();
    const Integrator *integrator = scene->getIntegrator();

    Point2i offset = block.getOffset();
    Vector2i size = block.getSize();

    /* Clear the block contents */
    block.clear();

    /* For each pixel and pixel sample sample */
    for (int y = 0; y < size.y(); ++y) {
        for (int x = 0; x < size.x(); ++x) {
            Eigen::MatrixXf vars(sampler->getSampleCount(), 3);

            for (uint32_t i = 0; i < sampler->getSampleCount(); ++i) {
                Point2f pixelSample = Point2f((float) (x + offset.x()), (float) (y + offset.y())) + sampler->next2D();
                Point2f apertureSample = sampler->next2D();

                /* Sample a ray from the camera */
                Ray3f ray;
                Color3f value = camera->sampleRay(ray, pixelSample, apertureSample);

                /* Qiyuan: For volumetric rendering, set the limit of ray.maxt */
                /* Note that ray.d is normalized */
//                ray.maxt = std::min(ray.maxt, Scene_Boundary);

                if (scene->getMedium() != nullptr)
                    ray.m = scene->getMedium();

                /* Compute the incident radiance */
                value *= integrator->Li(scene, sampler, ray);

                /* Store in the image block */
                block.put(pixelSample, value);

                /* Keep records of variance*/
                if (keepVariances)
                    vars.row(i) = value;
            }

            if (keepVariances) {
                Eigen::RowVector3f mean = vars.colwise().mean();
                vars = vars.rowwise() - mean;
                Eigen::VectorXf norms = vars.rowwise().squaredNorm();
                variances.push_back(sqrt(norms.mean()));
            }

        }
    }
}

static void render(Scene *scene, const std::string &filename) {
    const Camera *camera = scene->getCamera();
    Vector2i outputSize = camera->getOutputSize();
    scene->getIntegrator()->preprocess(scene);

    /* Create a block generator (i.e. a work scheduler) */
    BlockGenerator blockGenerator(outputSize, NORI_BLOCK_SIZE);

    /* Allocate memory for the entire output image and clear it */
    ImageBlock result(outputSize, camera->getReconstructionFilter());
    result.clear();

    /* Create a window that visualizes the partially rendered result */
    NoriScreen *screen = nullptr;
    if (gui) {
        nanogui::init();
        screen = new NoriScreen(result);
    }

    /* Do the following in parallel and asynchronously */
    std::thread render_thread = std::thread([&] {
        tbb::task_scheduler_init init(threadCount);

        cout << "Rendering .. " << endl;
        cout.flush();
        Timer timer;

        tbb::blocked_range<int> range(0, blockGenerator.getBlockCount());

        int allBlocks = blockGenerator.getBlockCount();

        auto map = [&](const tbb::blocked_range<int> &range) {
            /* Allocate memory for a small image block to be rendered
               by the current thread */
            ImageBlock block(Vector2i(NORI_BLOCK_SIZE),
                             camera->getReconstructionFilter());

            /* Create a clone of the sampler for the current thread */
            std::unique_ptr<Sampler> sampler(scene->getSampler()->clone());

            for (int i = range.begin(); i < range.end(); ++i) {
                /* Request an image block from the block generator */
                blockGenerator.next(block);


                /* Inform the sampler about the block to be rendered */
                sampler->prepare(block);

                /* Render all contained pixels */
                renderBlock(scene, sampler.get(), block);

                /* The image block has been processed. Now add it to
                   the "big" block that represents the entire image */
                result.put(block);

                int remains = blockGenerator.getBlockCount();
                if (remains % (allBlocks / 10) == 0 && remains != 0) {
                    float finished = (allBlocks - remains) * 100.f / allBlocks;
                    std::cout << finished << "% (took "
                              << timer.elapsedString() << "), " << remains << " blocks to be rendered ..." << std::endl;
                }
            }
        };

        /// Default: parallel rendering
        tbb::parallel_for(range, map);

        /// (equivalent to the following single-threaded call)
        // map(range);

        cout << "done. (took " << timer.elapsedString() << ")" << endl;
    });

    pthread_t phandle = render_thread.native_handle();
    render_thread.detach();

    /* Enter the application main loop */
    if (gui) {
        nanogui::mainloop(50.f);
        pthread_cancel(phandle);
        delete screen;
        nanogui::shutdown();
    } else {
        throw NoriException("Must use GUI for now!");
    }

    /* Shut down the user interface */
//    render_thread.join();

//    if (gui) {
//        delete screen;
//        nanogui::shutdown();
//    }

    /* Now turn the rendered image block into
       a properly normalized bitmap */
    std::unique_ptr<Bitmap> bitmap(result.toBitmap());

    /* Determine the filename of the output bitmap */
    std::string outputName = filename;
    size_t lastdot = outputName.find_last_of(".");
    if (lastdot != std::string::npos)
        outputName.erase(lastdot, std::string::npos);

    /* Save using the OpenEXR format */
    bitmap->saveEXR(outputName);

    /* Save tonemapped (sRGB) output using the PNG format */
    bitmap->savePNG(outputName);
}

int main(int argc, char **argv) {
    if (argc < 2) {
        cerr << "Syntax: " << argv[0] << " <scene.xml> [--no-gui] [--threads N]" << endl;
        return -1;
    }

    std::string sceneName = "";
    std::string exrName = "";

    for (int i = 1; i < argc; ++i) {
        std::string token(argv[i]);
        if (token == "-t" || token == "--threads") {
            if (i + 1 >= argc) {
                cerr << "\"--threads\" argument expects a positive integer following it." << endl;
                return -1;
            }
            threadCount = atoi(argv[i + 1]);
            i++;
            if (threadCount <= 0) {
                cerr << "\"--threads\" argument expects a positive integer following it." << endl;
                return -1;
            }

            continue;
        } else if (token == "--no-gui") {
            gui = false;
            continue;
        }

        filesystem::path path(argv[i]);

        try {
            if (path.extension() == "xml") {
                sceneName = argv[i];

                /* Add the parent directory of the scene file to the
                   file resolver. That way, the XML file can reference
                   resources (OBJ files, textures) using relative paths */
                getFileResolver()->prepend(path.parent_path());
            } else if (path.extension() == "exr") {
                /* Alternatively, provide a basic OpenEXR image viewer */
                exrName = argv[i];
            } else {
                cerr << "Fatal error: unknown file \"" << argv[i]
                     << "\", expected an extension of type .xml or .exr" << endl;
            }
        } catch (const std::exception &e) {
            cerr << "Fatal error: " << e.what() << endl;
            return -1;
        }
    }

    if (exrName != "" && sceneName != "") {
        cerr << "Both .xml and .exr files were provided. Please only provide one of them." << endl;
        return -1;
    } else if (exrName == "" && sceneName == "") {
        cerr << "Please provide the path to a .xml (or .exr) file." << endl;
        return -1;
    } else if (exrName != "") {
        if (!gui) {
            cerr << "Flag --no-gui was set. Please remove it to display the EXR file." << endl;
            return -1;
        }
        try {
            Bitmap bitmap(exrName);
            ImageBlock block(Vector2i((int)bitmap.cols(), (int) bitmap.rows()), nullptr);
            block.fromBitmap(bitmap);
            nanogui::init();
            NoriScreen *screen = new NoriScreen(block);
            nanogui::mainloop(50.f);
            delete screen;
            nanogui::shutdown();
        } catch (const std::exception &e) {
            cerr << e.what() << endl;
            return -1;
        }
    } else {  // sceneName != ""
        if (threadCount < 0) {
            threadCount = tbb::task_scheduler_init::automatic;
        }
        try {
            std::unique_ptr<NoriObject> root(loadFromXML(sceneName));
            /* When the XML root object is a scene, start rendering it .. */
//            if (root->getClassType() == NoriObject::EScene)
//                render(static_cast<Scene *>(root.get()), sceneName);

            if (root->getClassType() == NoriObject::EScene) {
                Scene *s = static_cast<Scene *>(root.get());
                if (s->getDenoiser() != nullptr)
                    s->getDenoiser()->run(s, sceneName, gui, threadCount);
                else
                    render(static_cast<Scene *>(root.get()), sceneName);
            }


        } catch (const std::exception &e) {
            cerr << e.what() << endl;
            return -1;
        }
    }

    // Create a 2*2 EXR
    // Bitmap bm(Vector2i(2,2));
    // bm(0, 0) = Color3f(60, 60, 60);
    // bm(0, 1) = Color3f(20, 20, 20);
    // bm(1, 0) = Color3f(4, 4, 4);
    // bm(1, 1) = Color3f(16, 16, 16);
    // bm.saveEXR("2by2");

    // Create a 2*2 EXR
    // Bitmap bm(Vector2i(4,4));
    // for (int r = 0; r < 4; ++r) { 
    //     for (int c = 0; c < 4; ++c) {
    //         bm(r, c) = Color3f(r * r + c * c + 10);
    //     }
    // }
    // bm.saveEXR("4by4");

//    if (keepVariances) {
//        double sum = 0.0;
//        for (double i: variances)
//            sum += i;
//        sum /= variances.size();
//        cout << "The average variance (RMSE) is " << sum << endl;
//    }
    return 0;
}
