/**
 * gauss_mixture.h
 *
 * Author:
 * vanilla
 *
 * Abstract:
 * @1. implementation of the Gaussian mixture model background subtraction
 * @2. The code is based on the following papers:
 *     a. Improved adaptive Gausian mixture model for background subtraction
 *     b. Efficient Adaptive Density Estimapion per Image Pixel for the Task of 
 *        Background Subtraction
 *     c. Recursive unsupervised learning of finite mixture models
 * @3. This c implementation is based on the c++ implementation of OpenCV2.3.1
 *     modules\video\src\bgfg_gaussmix2.cpp
 */

#ifndef _GAUSS_MIXTURE_H_
#define _GAUSS_MIXTURE_H_

#include "gauss_mixture_param.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct gaussmix_single_gaussian_s
{
        float gsg_fweight;
        float gsg_fmean[GAUSSMIX_DEFAULT_CHANNELS]; 
        float gsg_fvariance;
} gaussmix_single_gaussian_t;

typedef struct gaussmix_image_s
{
        int gi_iwidth;
        int gi_iheight;
        unsigned char* gi_ucdata;
} gaussmix_image_t;


gaussmix_image_t* gauss_mixture_create_image(int width, int height, 
                                             int channels);

void gauss_mixture_release_image(gaussmix_image_t **image);

/**
 * gauss_mixture_initialize --- initialize the Gaussian background model
 *
 * Parameters:
 * @img_width -------- the processed image width
 * @img_height ------- the processed image height
 * @bg_model --------- double pointer to the background model
 * @bg_model_used ---- double pointer to the array which indicates the number 
 *                     of modes used by each pixel
 *
 * Return:
 * return 1 if initialization succeed, -1 if failed.
 *
 * Note:
 * If initialization succeed, *bg_model and *bg_model_used will be set to point
 * to the corresponding memory.
 * If initialization failed, *bg_model and *bg_model_used will be set to NULL.
 */
int gauss_mixture_initialize(int img_width, int img_height, 
                             float **bg_model,
                             unsigned char **bg_model_used);

void gauss_mixture_final(float **bg_model, 
                         unsigned char **bg_model_used);

/**
 * gauss_mixture_update ------- update the Gaussian background model based on 
 *                              the input image, and output the foreground mask.
 *
 * Parameters:
 * @image ----------- the input 8-bit RGB image
 * @fg_mask --------- the output 8-bit binary image, foreground mask
 * @bg_model -------- point to the background model
 * @bg_model_used --- pointer to the array which indicates the number 
 *                    of modes used by each pixel
 */
void gauss_mixture_update(gaussmix_image_t *image, 
                          gaussmix_image_t *fg_mask, 
                          float *bg_model, 
                          unsigned char *bg_model_used);

#ifdef __cplusplus
}
#endif

#endif /*_GAUSS_MIXTURE_H_*/
