//
// Created by Qiyuan Dong on 12.05.22.
//

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
#include <nori/object.h>
#include <nori/common.h>
#include <nori/denoiser.h>
#include <thread>
#include <vector>
#include <unordered_map>

NORI_NAMESPACE_BEGIN
    typedef Eigen::Array<Color3f, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXc;

    class AdaptiveNLM : public Denoiser {
    private:
//        int m_rows, m_cols;
//        int m_spp;
        int m_its;
        long m_total_samples, m_it_samples;

        int m_it_r, m_final_r;
        int m_it_f, m_final_f;
        float m_it_alpha, m_final_alpha;
        float m_it_k, m_final_k;

        MatrixXc m_buffer_A, m_buffer_B;
        MatrixXc m_filtered_A, m_filtered_B;

        MatrixXc m_empvar_A, m_empvar_B, m_empvar, m_empvar_var;
        MatrixXc m_sum, m_sum_A, m_sum_B;
        MatrixXc m_sum2, m_sum2_A, m_sum2_B;
        MatrixXc m_var, m_var_cache;
        MatrixXc m_result;

        MatrixXf m_NLM_cache;
        MatrixXf m_error;
        MatrixXf m_error_w;
        MatrixXf m_sampling_map;

        MatrixXu m_sampling_cnt;
        MatrixXu m_sampling_cnt_A;
        MatrixXu m_sampling_cnt_B;

        MatrixXf m_gaussian_k;

        int m_threadCnt;
        int m_strength;
        bool initialized = false;

    public:

        explicit AdaptiveNLM(const PropertyList &propList) {
//            m_cols = propList.getInteger("width");
//            m_rows = propList.getInteger("height");
            m_strength = propList.getInteger("strength", 2);
            m_its = propList.getInteger("its");
//            m_spp = propList.getInteger("spp");
            m_total_samples = (long) m_spp * m_rows * m_cols;
            m_it_samples = m_total_samples / m_its;

            m_it_r = propList.getInteger("it_r", 7);
            m_final_r = propList.getInteger("final_r", 10);
            m_it_f = propList.getInteger("it_f", 3);
            m_final_f = propList.getInteger("final_f", 3);

            m_it_alpha = propList.getFloat("it_alpha", 0.5);
            m_final_alpha = propList.getFloat("final_alpha", 1.f);
            m_it_k = propList.getFloat("it_k", 0.45);
            m_final_k = propList.getFloat("final_k", 0.45);

        }


        void setConstants(int rows, int cols, int spp) override {
            m_rows = rows;
            m_cols = cols;
            m_spp = spp;

            m_buffer_A.resize(m_rows, m_cols);
            m_buffer_B.resize(m_rows, m_cols);
            m_filtered_A.resize(m_rows, m_cols);
            m_filtered_B.resize(m_rows, m_cols);

            m_empvar_A.resize(m_rows, m_cols);
            m_empvar_B.resize(m_rows, m_cols);
            m_empvar.resize(m_rows, m_cols);
            m_empvar_var.resize(m_rows, m_cols);
            m_var.resize(m_rows, m_cols);
            m_var_cache.resize(m_rows, m_cols);

            m_sum.resize(m_rows, m_cols);
            m_sum_A.resize(m_rows, m_cols);
            m_sum_B.resize(m_rows, m_cols);
            m_sum2.resize(m_rows, m_cols);
            m_sum2_A.resize(m_rows, m_cols);
            m_sum2_B.resize(m_rows, m_cols);

            m_NLM_cache.resize(m_rows, m_cols);
            m_error.resize(m_rows, m_cols);
            m_error_w.resize(m_rows, m_cols);
            m_sampling_map.resize(m_rows, m_cols);
            m_sampling_cnt.resize(m_rows, m_cols);
            m_sampling_cnt_A.resize(m_rows, m_cols);
            m_sampling_cnt_B.resize(m_rows, m_cols);

            m_result.resize(m_rows, m_cols);

            // initialize sampling map to uniform sampling
            m_sampling_map.array() = (float) m_spp / (m_its * 2.0);

            // initialize gaussian kernel (5*5, sigma=0.8)
            m_gaussian_k.resize(5, 5);
            m_gaussian_k <<
                         0.00048091, 0.00501119, 0.01094545, 0.00501119, 0.00048091,
                    0.00501119, 0.0522178, 0.11405416, 0.0522178, 0.00501119,
                    0.01094545, 0.11405416, 0.2491172, 0.11405416, 0.01094545,
                    0.00501119, 0.0522178, 0.11405416, 0.0522178, 0.00501119,
                    0.00048091, 0.00501119, 0.01094545, 0.00501119, 0.00048091;
           initialized = true;
        }



        inline bool isValidCoord(const int &r, const int &c) {
            return (r >= 0 && r < m_rows && c >= 0 && c < m_cols);
        }

        inline Color3f cwiseMin(const Color3f &c1, const Color3f &c2) {
            return Color3f(
                    c1(0) > c2(0) ? c2(0) : c1(0),
                    c1(1) > c2(1) ? c2(1) : c1(1),
                    c1(2) > c2(2) ? c2(2) : c1(2)
            );
        }

        void
        NLM_filtering_block(const MatrixXc &src, MatrixXc &dst, const MatrixXc &weights, const MatrixXc &var, int r,
                            int f,
                            float alpha, float k, ImageBlock &block) {

            Point2i offset = block.getOffset();
            Vector2i size = block.getSize();
            std::unordered_map<unsigned long, float> hashTable;
            hashTable.reserve(size.x() * size.y() * size.x() * size.y() / 4);


            // for each pixel p(r, c)
//            for (int rp = 0; rp < m_rows; rp++) {
//                for (int cp = 0; cp < m_cols; cp++) {
            for (int br = 0; br < size.y(); ++br) {
                int rp = offset.y() + br;

                for (int bc = 0; bc < size.x(); ++bc) {
                    int cp = offset.x() + bc;

                    // for each pixel q in Neighbor(p)
                    for (int rq = rp - r; rq <= rp + r; rq++) {
                        for (int cq = cp - r; cq <= cp + r; cq++) {
                            if (!isValidCoord(rq, cq)) continue;


                            // look up the hash table first
                            unsigned long idx_p = rp * m_rows + cp;
                            unsigned long idx_q = rq * m_rows + cq;

                            float w = 0.f;
                            auto search = hashTable.find((idx_q << 32) | idx_p);
                            if (search != hashTable.end()) {
                                w = search->second;
                            } else {
                                // compute the distance between two patches
                                float d2 = 0;
                                int cnt = 0;
                                // for each offset
                                for (int dr = -f; dr <= f; dr++) {
                                    for (int dc = -f; dc <= f; dc++) {
                                        int rrp = rp + dr, ccp = cp + dc;
                                        int rrq = rq + dr, ccq = cq + dc;
                                        if (!isValidCoord(rrp, ccp) || !isValidCoord(rrq, ccq)) continue;

//                                        Color3f d = (weights(rrp, ccp) - weights(rrq, ccq)).square() -
//                                                    alpha * (var(rrp, ccp) + cwiseMin(var(rrp, ccp), var(rrq, ccq)));

                                        // modify the calculation of d so that d(p,q) = d(q,p)
                                        Color3f d = (weights(rrp, ccp) - weights(rrq, ccq)).square() -
                                                    alpha * ((var(rrp, ccp) + var(rrq, ccq)) / 2.f +
                                                             cwiseMin(var(rrp, ccp), var(rrq, ccq)));

                                        d /= 1e-10 + k * k * (var(rrp, ccp) + var(rrq, ccq));

                                        d2 += d.sum() / 3.f;

                                        cnt++;
                                    }
                                }
                                d2 /= (float) cnt * 3.f;
                                w = std::exp(-1.f * std::max(0.f, d2));

                                // preserve details in noisy area
                                if (w < 0.05)
                                    w = 0;

                                hashTable[(idx_p << 32) | idx_q] = w;

                                if (w == 0) continue;
                            }

                            // record w in error_w, will be used in computing sampling map
                            m_error_w(rp, cp) += w;

                            // for each offset of p, store weighted contribution
                            for (int dr = -f; dr <= f; dr++) {
                                for (int dc = -f; dc <= f; dc++) {
                                    int rrp = rp + dr, ccp = cp + dc;
                                    int rrq = rq + dr, ccq = cq + dc;
                                    if (!isValidCoord(rrp, ccp) || !isValidCoord(rrq, ccq)) continue;
                                    m_NLM_cache(rrp, ccp) += w;
                                    dst(rrp, ccp) += w * src(rrq, ccq);
                                }
                            }
                        }
                    }
                }
            }

//            for (int br = 0; br < size.y(); ++br) {
//                int rp = offset.y() + br;
//                for (int bc = 0; bc < size.x(); ++bc) {
//                    int cp = offset.x() + bc;
//
//                    dst(rp, cp) = dst(rp, cp) / m_NLM_cache(rp, cp);
//                }
//            }
        }


        void
        NLM_filtering(const MatrixXc &src, MatrixXc &dst, const MatrixXc &weights, const MatrixXc &var, int r, int f,
                      float alpha, float k) {
            dst.setZero();
            m_NLM_cache.setZero();
            BlockGenerator blockGenerator(Vector2i(m_cols, m_rows), NORI_BLOCK_SIZE);
            tbb::task_scheduler_init init(m_threadCnt);

            tbb::blocked_range<int> range(0, blockGenerator.getBlockCount());

            auto map = [&](const tbb::blocked_range<int> &range) {
                ImageBlock block(Vector2i(NORI_BLOCK_SIZE * 2), nullptr);

                for (int i = range.begin(); i < range.end(); ++i) {
                    /* Request an image block from the block generator */
                    blockGenerator.next(block);
                    NLM_filtering_block(src, dst, weights, var, r, f, alpha, k, block);
                }
            };

            /// Default: parallel rendering
            tbb::parallel_for(range, map);

            for (int rp = 0; rp < m_rows; rp++)
                for (int cp = 0; cp < m_cols; cp++)
                    dst(rp, cp) = dst(rp, cp) / m_NLM_cache(rp, cp);
        }

        void compute_var_from_record(const MatrixXu &cnt, const MatrixXc &sum, const MatrixXc &sum2, MatrixXc &var) {
            for (int r = 0; r < m_rows; r++) {
                for (int c = 0; c < m_cols; c++) {
                    float n = (float) cnt(r, c);
                    Color3f mean = sum(r, c) / n;
                    var(r, c) = sum2(r, c) / n - mean * mean;
                }
            }
        }

        void compute_var() {
            for (int r = 0; r < m_rows; r++)
                for (int c = 0; c < m_cols; c++) {
                    /* Obtain an initial estimate ∆[p] of the pixel variance Var [p]
                     * using the squared difference between the two buffers */
                    m_var(r, c) = (m_buffer_A(r, c) - m_buffer_B(r, c)).square() / 2.f;

                    /* To obtain the filter weights we also need the variance of Σ[p],
                     * which we compute as the difference between Σ[p] from buffers A and B */
                    m_empvar_var(r, c) = (m_empvar_A(r, c) - m_empvar_B(r, c)).abs();
                }

            /* Smooth the initial estimates ∆[p] by cross filtering them with the empirical sample variances Σ[p],
             * that is, we compute NL-Means filter weights using Σ[p] and apply them to filter ∆[p] */
            NLM_filtering(m_var, m_var_cache, m_empvar, m_empvar_var, 1, 3, 4, 0.45);
            m_var = m_var_cache;

            /* Since the filtering tends to increase the variance estimate of pixels neighboring noisy regions,
             * we also clamp the filtering value of ∆[p] so that it does not exceed Σ[p] */
            for (int r = 0; r < m_rows; r++)
                for (int c = 0; c < m_cols; c++)
                    m_var(r, c) = cwiseMin(m_var(r, c), m_empvar(r, c));
        }

        void compute_sampling_map() {
            for (int r = 0; r < m_rows; r++) {
                for (int c = 0; c < m_cols; c++) {
                    m_error(r, c) = (float) ((m_buffer_A(r, c) - m_buffer_B(r, c)).square().sum())
                                    / (1e-3 + m_buffer_A(r, c).square().sum());

                    m_error(r, c) *= m_error_w(r, c) / (1.f + (float) m_sampling_cnt(r, c));

                }
            }


            m_sampling_map.setZero();
            // gaussian filter
            int k_radius = 2;
            double sum = 0;
            for (int r = 0; r < m_rows; r++) {
                for (int c = 0; c < m_cols; c++) {
                    // do not filter border
                    if (r < k_radius || r >= m_rows - k_radius || c < k_radius || c >= m_cols - k_radius) {
                        m_sampling_map(r, c) = m_error(r, c);
                    } else {
                        for (int kr = -k_radius; kr <= k_radius; kr++)
                            for (int kc = -k_radius; kc <= k_radius; kc++)
                                m_sampling_map(r, c) +=
                                        m_error(r + kr, c + kc) * m_gaussian_k(kr + k_radius, kc + k_radius);
                    }

                    sum += m_sampling_map(r, c);
                }
            }

            // sampling map from error
            float samples = (float) m_it_samples / 2.f;
            m_sampling_map *= samples / sum;
            m_sampling_map = m_sampling_map.array().ceil();

            // clamp
            int it_spp = m_spp / (2 * m_its);
            m_sampling_map.array().min((float) it_spp / m_strength);
            m_sampling_map.array().max((float) it_spp * m_strength);

            // normalized
            sum = m_sampling_map.array().sum();
            m_sampling_map *= samples / sum;
            m_sampling_map = m_sampling_map.array().ceil();

        }

        void update_result() {
            for (int r = 0; r < m_rows; r++)
                for (int c = 0; c < m_cols; c++)
                    m_result(r, c) = (m_filtered_A(r, c) + m_filtered_B(r, c)) / 2.f;
        }

        void run(Scene *scene, const std::string &filename, bool gui, int threadCount) override {
            assert(initialized);

            m_threadCnt = threadCount;

            /* Determine the filename of the output bitmap */
            std::string outputName = filename;
            size_t lastdot = outputName.find_last_of(".");
            if (lastdot != std::string::npos)
                outputName.erase(lastdot, std::string::npos);

            const Camera *camera = scene->getCamera();
            Vector2i outputSize = camera->getOutputSize();
            scene->getIntegrator()->preprocess(scene);

            /* Allocate memory for the entire output image and clear it */
            ImageBlock result(outputSize, camera->getReconstructionFilter());
            ImageBlock buffer_A(outputSize, camera->getReconstructionFilter());
            ImageBlock buffer_B(outputSize, camera->getReconstructionFilter());
            result.clear();
            buffer_A.clear();
            buffer_B.clear();

            m_sampling_cnt.setZero();
            m_sampling_cnt_A.setZero();
            m_sampling_cnt_B.setZero();
            m_buffer_A.setZero();
            m_buffer_B.setZero();
            m_empvar.setZero();
            m_empvar_A.setZero();
            m_empvar_B.setZero();
            m_sum.setZero();
            m_sum_A.setZero();
            m_sum_B.setZero();
            m_sum2.setZero();
            m_sum2_A.setZero();
            m_sum2_B.setZero();

            /* Create a window that visualizes the partially rendered result */
            NoriScreen *screen = nullptr;

            if (gui) {
                nanogui::init();
                screen = new NoriScreen(result);
            }

            Timer t;
            std::thread it_thread(
                    [&] {
                        for (int it = 1; it <= m_its; it++) {

                            Timer timer;
                            cout << "Adaptive rendering with NLM, iter #" << it << endl;
                            render(scene, result, buffer_A, buffer_B);

                            Timer timer1;
                            cout << "Computing variance ...";
                            compute_var();
                            cout << "done (took " << timer1.elapsedString() << ")" << endl;

                            m_error_w.setZero();

                            Timer timer2;
                            cout << "NLM filtering ...";
                            if (it != m_its) {
                                NLM_filtering(m_buffer_A, m_filtered_A, m_buffer_B, m_var, m_it_r, m_it_f, m_it_alpha,
                                              m_it_k);
                                NLM_filtering(m_buffer_B, m_filtered_B, m_buffer_A, m_var, m_it_r, m_it_f, m_it_alpha,
                                              m_it_k);
                            } else {
                                NLM_filtering(m_buffer_A, m_filtered_A, m_buffer_B, m_var, m_final_r, m_final_f,
                                              m_final_alpha,
                                              m_final_k);
                                NLM_filtering(m_buffer_B, m_filtered_B, m_buffer_A, m_var, m_final_r, m_final_f,
                                              m_final_alpha,
                                              m_final_k);
                            }
                            cout << "done (took " << timer2.elapsedString() << ")" << endl;

                            if (it != m_its) {
                                Timer timer3;
                                cout << "Computing sampling map ...";
                                compute_sampling_map();
                                cout << "done (took " << timer3.elapsedString() << ")" << endl;

                                // save sampling map to a .exr file
                                Bitmap bm(Vector2i(m_cols, m_rows));
                                for (int r = 0; r < m_rows; r++)
                                    for (int c = 0; c < m_cols; c++)
                                        bm.coeffRef(r, c) = m_sampling_map(r, c);
                                std::string sm_outputName = outputName + "-sm-" + std::to_string(it);
                                bm.saveEXR(sm_outputName);
                            }

                            cout << "Iter #" << it << " done (took " << timer.elapsedString() << ")" << endl;
                        }
                    });

            if (gui)
                nanogui::mainloop(50.f);

            it_thread.join();

            if (gui) {
                delete screen;
                nanogui::shutdown();
            }


            cout << "Adaptive rendering with NLM took " << t.elapsedString() << endl;

//            /* Enter the application main loop */
//            if (gui)
//                nanogui::mainloop(50.f);
//
//            /* Shut down the user interface */
//            render_thread.join();
//
//            if (gui) {
//                delete screen;
//                nanogui::shutdown();
//            }
            /* Now turn the rendered image block into a properly normalized bitmap */

            // filtered final result
            update_result();    // merge m_filtered_A and m_filtered_B
            Bitmap bitmap(Vector2i(m_cols, m_rows)

            );
            for (
                    int r = 0;
                    r < m_rows;
                    r++)
                for (
                        int c = 0;
                        c < m_cols;
                        c++)
                    bitmap.
                            coeffRef(r, c
                    ) =
                            m_result(r, c
                            );
            /* Save using the OpenEXR format */
            bitmap.
                    saveEXR(outputName
                            + "-nlm");
            /* Save tonemapped (sRGB) output using the PNG format */
            bitmap.
                    savePNG(outputName
                            + "-nlm");

            // unfiltered final result
            result.
                    saveToMatrix(m_result);
            for (
                    int r = 0;
                    r < m_rows;
                    r++)
                for (
                        int c = 0;
                        c < m_cols;
                        c++)
                    bitmap.
                            coeffRef(r, c
                    ) =
                            m_result(r, c
                            );
            /* Save using the OpenEXR format */
            bitmap.
                    saveEXR(outputName
                            + "-ori");
            /* Save tonemapped (sRGB) output using the PNG format */
            bitmap.
                    savePNG(outputName
                            + "-ori");

        }

        void renderBlock(const Scene *scene, Sampler *sampler, ImageBlock &block, const char &buffer) {
            assert(buffer == 'A' || buffer == 'B');

            const Camera *camera = scene->getCamera();
            const Integrator *integrator = scene->getIntegrator();

            Point2i offset = block.getOffset();
            Vector2i size = block.getSize();


            /* Clear the block contents */
            block.clear();

            /* For each pixel and pixel sample sample */
            for (int y = 0; y < size.y(); ++y) {
                for (int x = 0; x < size.x(); ++x) {

                    // Maybe problematic about row, col, x, y
                    int pr = offset.y() + y;
                    int pc = offset.x() + x;
                    uint32_t samples = (int) m_sampling_map(pr, pc);

                    for (uint32_t i = 0; i < samples; ++i) {
                        Point2f pixelSample =
                                Point2f((float) (x + offset.x()), (float) (y + offset.y())) + sampler->next2D();
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

                        /* Keep records of mean, variance, count*/
                        m_sampling_cnt(pr, pc) += 1;
                        if (buffer == 'A') {
                            m_sampling_cnt_A(pr, pc) += 1;
                            m_sum_A(pr, pc) += value;
                            m_sum2_A(pr, pc) += value.square();
                        } else {
                            m_sampling_cnt_B(pr, pc) += 1;
                            m_sum_B(pr, pc) += value;
                            m_sum2_B(pr, pc) += value.square();
                        }
                        m_sampling_cnt(pr, pc) += 1;
                        m_sum(pr, pc) += value;
                        m_sum2(pr, pc) += value.square();
                    }
//                    if (keepVariances) {
//                        Eigen::RowVector3f mean = vars.colwise().mean();
//                        vars = vars.rowwise() - mean;
//                        Eigen::VectorXf norms = vars.rowwise().squaredNorm();
//                        variances.push_back(sqrt(norms.mean()));
//                    }

                }
            }
        }

        void render(Scene *scene, ImageBlock &result, ImageBlock &buffer_A, ImageBlock &buffer_B) {
            const Camera *camera = scene->getCamera();
            Vector2i outputSize = camera->getOutputSize();
//            scene->getIntegrator()->preprocess(scene);

            /* Create a block generator (i.e. a work scheduler) */
            BlockGenerator blockGenerator(outputSize, NORI_BLOCK_SIZE);

            /* Do the following in parallel and asynchronously */
            std::thread render_thread([&] {
                tbb::task_scheduler_init init(m_threadCnt);

                cout << "Rendering ... ";
                cout.flush();
                Timer timer;

                tbb::blocked_range<int> range(0, blockGenerator.getBlockCount());

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
                        renderBlock(scene, sampler.get(), block, 'A');

                        /* The image block has been processed. Now add it to
                           the "big" block that represents the entire image */
                        result.put(block);
                        buffer_A.put(block);

//                        sampler->prepare(block);
                        renderBlock(scene, sampler.get(), block, 'B');
                        result.put(block);
                        buffer_B.put(block);
                    }
                };

                /// Default: parallel rendering
                tbb::parallel_for(range, map);

                /// (equivalent to the following single-threaded call)
                // map(range);
                cout << "done (took " << timer.elapsedString() << ")" << endl;
            });

//            /* Enter the application main loop */
//            if (gui)
//                nanogui::mainloop(50.f);
//
//            /* Shut down the user interface */
            render_thread.join();


            // save to buffer matrix
            buffer_A.saveToMatrix(m_buffer_A);
            buffer_B.saveToMatrix(m_buffer_B);

            // compute variance
            compute_var_from_record(m_sampling_cnt_A, m_sum_A, m_sum2_A, m_empvar_A);
            compute_var_from_record(m_sampling_cnt_B, m_sum_B, m_sum2_B, m_empvar_B);
            compute_var_from_record(m_sampling_cnt, m_sum, m_sum2, m_empvar);
//
//            if (gui) {
//                delete screen;
//                nanogui::shutdown();
//            }
//
//            /* Now turn the rendered image block into
//               a properly normalized bitmap */
//            std::unique_ptr<Bitmap> bitmap(result.toBitmap());
//
//            /* Determine the filename of the output bitmap */
//            std::string outputName = filename;
//            size_t lastdot = outputName.find_last_of(".");
//            if (lastdot != std::string::npos)
//                outputName.erase(lastdot, std::string::npos);
//
//            /* Save using the OpenEXR format */
//            bitmap->saveEXR(outputName);
//
//            /* Save tonemapped (sRGB) output using the PNG format */
//            bitmap->savePNG(outputName);
        }


        [[nodiscard]] std::string toString() const

        override {
            return tfm::format(
                    "AdaptiveNLM [\n"
                    "]");
        }


    };

    NORI_REGISTER_CLASS(AdaptiveNLM, "adaptive_NLM");

NORI_NAMESPACE_END