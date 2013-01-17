#pragma warning(disable:4996)

#include <iostream>
#include <opencv2/opencv.hpp>
#include "gauss_mixture.h"

/*
int main()
{
        float *bg_model = NULL;
        unsigned char *bg_model_used = NULL;
        int width = 704;
        int height = 576;

        gauss_mixture_initialize(width, height, &bg_model, &bg_model_used);

        cv::Mat cv_img = cv::imread("./test_resource/001.jpg");
        gaussmix_image_t *gauss_img = 
                gauss_mixture_create_image(width, height, 3);
        memcpy(gauss_img->gi_ucdata, cv_img.data, sizeof(unsigned char) * 
                width * height * cv_img.channels());
        gaussmix_image_t *mask = 
                gauss_mixture_create_image(width, height, 1);

        gauss_mixture_update(gauss_img, mask, bg_model, bg_model_used);

        cv::Mat cv_img2 = cv::imread("./test_resource/002.jpg");
        memcpy(gauss_img->gi_ucdata, cv_img2.data, sizeof(unsigned char) * 
                width * height * cv_img.channels());
        gauss_mixture_update(gauss_img, mask, bg_model, bg_model_used);
        

        gauss_mixture_release_image(&gauss_img);
        gauss_mixture_final(&bg_model, &bg_model_used);
        
        return 0;
}
*/

static __inline unsigned __int64 GetCycleCount()
{
        __asm _emit 0x0F
        __asm _emit 0x31
}

int main()
{
        float *bg_model = NULL;
        unsigned char *bg_model_used = NULL;
        int width = 704;
        int height = 576;

        gauss_mixture_initialize(width, height, &bg_model, &bg_model_used);
        gaussmix_image_t *gauss_img = 
                gauss_mixture_create_image(width, height, 3);
        gaussmix_image_t *gauss_mask = 
                gauss_mixture_create_image(width, height, 1);

        cv::VideoCapture video("./test_resource/shadow.avi");
        cv::Mat frame;
        cv::Mat fg_frame(height, width, CV_8U);
        int nframe = 0;
        unsigned __int64 total_time = 0;
        char szpath[256];
        memset(szpath, 0, sizeof(szpath));

        while (true) {
                video >> frame;
                if (!frame.data) break;

                memcpy(gauss_img->gi_ucdata, frame.data, sizeof(unsigned char) 
                        * width * height * frame.channels());


                unsigned __int64 start = GetCycleCount();
                gauss_mixture_update(gauss_img, gauss_mask, bg_model, 
                                     bg_model_used);
                unsigned __int64 end = GetCycleCount();
                total_time += (end - start);

                memcpy(fg_frame.data, gauss_mask->gi_ucdata, 
                        sizeof(unsigned char) * width * height);

                cv::imshow("foreground", fg_frame);
                cv::imshow("show", frame);
                if (cv::waitKey(1) == 27) break;

                sprintf(szpath, "./test_resource/without_shadow_detect/%d.jpg", nframe);
                cv::imwrite(szpath, fg_frame);

                nframe++;
        }

        gauss_mixture_release_image(&gauss_img);
        gauss_mixture_release_image(&gauss_mask);
        gauss_mixture_final(&bg_model, &bg_model_used);

        std::cout<<"average time is: "<<total_time * 1.0 / (2.67 * 1000000) / nframe
                 <<std::endl;

        return 0;
}



/*
int main()
{
        cv::VideoCapture video("./test_resource/11.10.8.13.36.avi");
        cv::Mat frame;

        video >> frame;
        std::cout<<frame.rows<<" "<<frame.cols<<std::endl;

        cv::BackgroundSubtractorMOG2 bg_model;
        while (true) {
                video >> frame;
                if (!frame.data) break;

                cv::imshow("show", frame);
                if (cv::waitKey(30) == 27) break;
        }

        return 0;
}
*/