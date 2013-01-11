#pragma warning(disable:4996)

#include <iostream>
#include <opencv2/opencv.hpp>
#include "gauss_mixture.h"

int main()
{
        float *bg_model = NULL;
        unsigned char *bg_model_used = NULL;
        int width = 704;
        int height = 576;

        gauss_mixture_initialize(width, height, &bg_model, &bg_model_used);

        cv::Mat cv_img = cv::imread("./test_resource/001.jpg");
        gaussmix_image_t *gauss_img = 
                gauss_mixture_create_image(width, height,3);
        memcpy(gauss_img->gi_ucdata, cv_img.data, sizeof(unsigned char) * 
                width * height * cv_img.channels());

        
        gauss_mixture_release_image(&gauss_img);
        gauss_mixture_final(&bg_model, &bg_model_used);
        
        return 0;
}