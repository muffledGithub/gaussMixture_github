#include "gauss_mixture.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

static gaussmix_model_param_t g_model_param;

void _set_model_param(int img_width, int img_height)
{
        memset(&g_model_param, 0, sizeof(gaussmix_model_param_t));

        g_model_param.gmp_iwidth = img_width;
        g_model_param.gmp_iheight = img_height;
        g_model_param.gmp_ichannels = GAUSSMIX_DEFAULT_CHANNELS;

        g_model_param.gmp_bpost_filtering = 1;
        g_model_param.gmp_dmin_area = GAUSSMIX_MINAREA;
        g_model_param.gmp_binit = 1;

        g_model_param.gmp_falpha = 1.0f / GAUSSMIX_WINDOW_SIZE;
        g_model_param.gmp_fcthr = GAUSSMIX_STD_THRESHOLD * 
                GAUSSMIX_STD_THRESHOLD;
        g_model_param.gmp_fthres_smd = GAUSSMIX_STD_THRESHOLD_GENERATE * 
                GAUSSMIX_STD_THRESHOLD_GENERATE;
        g_model_param.gmp_fone_minus_cf = GAUSSMIX_BACKGROUND_THRESHOLD;

        g_model_param.gmp_fvar_init = GAUSSMIX_VAR_INIT;
        g_model_param.gmp_fvar_max = GAUSSMIX_VAR_MAX;
        g_model_param.gmp_fvar_min = GAUSSMIX_VAR_MIN;

        g_model_param.gmp_fct = GAUSSMIX_CT;
        g_model_param.gmp_inmodels = GAUSSMIX_NGAUSSIANS;

        g_model_param.gmp_bshadow_detection = 1;
        g_model_param.gmp_ucshadow_value = GAUSSMIX_SHADOW_VALUE;
        g_model_param.gmp_ftau = GAUSSMIX_SHADOW_TAU;
}

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

        _set_model_param(img_width, img_height);

        return 1;

__failed:
        if (*bg_model) free(*bg_model);
        if (*bg_model_used) free(*bg_model_used);
        *bg_model = NULL;
        *bg_model_used = NULL;
        return -1;
}