#include "gauss_mixture.h"
#include <stdio.h>
#include <stdlib.h>

gaussmix_image_t* gauss_mixture_create_image(int width, int height)
{
        gaussmix_image_t *image = (gaussmix_image_t*)malloc(
                sizeof(gaussmix_image_t));
        if (!image) return NULL;

        image->gi_ucdata = (unsigned char*)calloc(width * height, 
                sizeof(unsigned char));
        if (!image->gi_ucdata) {
                free(image);
                return NULL;
        }

        image->gi_iwidth = width;
        image->gi_iheight = height;

        return image;
}

void gauss_mixture_release_image(gaussmix_image_t **image)
{
        if (*image) {
                if ( (*image)->gi_ucdata ) {
                        free( (*image)->gi_ucdata );
                }
                free(*image);
                (*image) = NULL;
        }
}

int gauss_mixture_initialize(int img_width, int img_height, 
                             int nmax_gaussians,
                             float **bg_model,
                             unsigned char **bg_model_used)
{
        if ( !(*bg_model) || !(*bg_model_used) ) goto __failed;

        /* the form 2 + GAUSSMIX_DEFAULT_CHANNELS refers to the struct
           gaussmix_single_gaussian_t */
        *bg_model = (float*)calloc(img_width * img_height * nmax_gaussians * 
                (2 + GAUSSMIX_DEFAULT_CHANNELS), sizeof(float) );
        if ( !(*bg_model) ) goto __failed;

        *bg_model_used = (unsigned char *)calloc(img_width * img_height, 
                sizeof(unsigned char));
        if ( !(*bg_model_used) ) goto __failed;

        return 1;

__failed:
        if (*bg_model) free(*bg_model);
        if (*bg_model_used) free(*bg_model_used);
        *bg_model = NULL;
        *bg_model_used = NULL;
        return -1;
}