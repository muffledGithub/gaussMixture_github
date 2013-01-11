#include "gauss_mixture.h"
#include <stdio.h>
#include <stdlib.h>
#include <memory.h>

static gaussmix_model_param_t g_model_param;

static void _set_model_param(int img_width, int img_height)
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

static __inline float _maha_distance(float r, float g, float b,
                                     float *mean)
{
        float maha_dis = 0.f;

        maha_dis += ( (mean[0] - r) * (mean[0] - r) );
        maha_dis += ( (mean[1] - g) * (mean[1] - g) );
        maha_dis += ( (mean[2] - b) * (mean[2] - b) );

        return maha_dis;
}

static __inline void _gaussian_update(gaussmix_single_gaussian_t *pgauss,
                                      float r, float g, float b,
                                      float maha_dis, 
                                      float weight, 
                                      float alpha)
{
        weight += alpha;
        pgauss->gsg_fmean[0] += alpha /
                weight * (r - pgauss->gsg_fmean[0]);
        pgauss->gsg_fmean[1] += alpha /
                weight * (g - pgauss->gsg_fmean[1]);
        pgauss->gsg_fmean[2] += alpha /
                weight * (b - pgauss->gsg_fmean[2]);
        pgauss->gsg_fvariance += alpha / 
                weight * (maha_dis - pgauss->gsg_fvariance);
        if (pgauss->gsg_fvariance < g_model_param.gmp_fvar_min) {
                pgauss->gsg_fvariance = 
                        g_model_param.gmp_fvar_min;
        }
        if (pgauss->gsg_fvariance > g_model_param.gmp_fvar_max) {
                pgauss->gsg_fvariance = 
                        g_model_param.gmp_fvar_max;
        }
        pgauss->gsg_fweight = weight;
}

static __inline _gaussian_sort(gaussmix_single_gaussian_t* pgauss,
                               int index)
{
        gaussmix_single_gaussian_t swap_temp;
        
        for (; index > 0; index--)
        {
                if (pgauss[index].gsg_fweight < pgauss[index - 1].gsg_fweight)
                        break;
                else {
                        swap_temp = pgauss[index];
                        pgauss[index] = pgauss[index - 1];
                        pgauss[index - 1] = swap_temp;
                }
        }
}

static __inline _generate_new_gaussian(gaussmix_single_gaussian_t* pgauss, 
                                       unsigned char *pngaussians_used, 
                                       float r, float g, float b,
                                       float alpha)
{
        int postion = 0;

        /* maximum positions have been used */
        if (*pngaussians_used == g_model_param.gmp_inmodels) {
                postion = (*pngaussians_used) - 1;
        }
        else {
                postion = (*pngaussians_used);
                (*pngaussians_used)++;
        }

        pgauss[postion].gsg_fweight = alpha;
        pgauss[postion].gsg_fmean[0] = r;
        pgauss[postion].gsg_fmean[1] = g;
        pgauss[postion].gsg_fmean[2] = b;
        pgauss[postion].gsg_fvariance = g_model_param.gmp_fvar_init;

        /* give the new gaussian correct position */
        _gaussian_sort(pgauss, postion);
}

/* calculate distances to the modes (+ sort)
   here we need to go in descending order!!! */
static unsigned char _update(float r, float g, float b, 
                             float falpha, 
                             gaussmix_single_gaussian_t *psg, 
                             unsigned char *pngaussians_used)
{
        /* return value, 1 => the pixel is classified as background */
        unsigned char bbackground = 0; 

        /* internal: */
        /* if remains 0, a new GMM model will be added */
        unsigned char bfits = 0; 
        float ftotal_weight = 0.f; /* used for background portion decision */
        float fsingle_weight = 0.f;
        float maha_dis = 0.f; /* Mahalanobis distance */
        int imodes = 0;
        for (; imodes < (*pngaussians_used); imodes++, psg++) {
                fsingle_weight = psg->gsg_fweight;
                fsingle_weight = (1 - falpha) * fsingle_weight
                        - falpha * g_model_param.gmp_fct;

                /* fit not found yet */
                if (!bfits) {
                        maha_dis = _maha_distance(r, g, b, psg->gsg_fmean);

                        /* background checking */
                        if ((ftotal_weight < g_model_param.gmp_fone_minus_cf) && 
                             (maha_dis < g_model_param.gmp_fcthr * 
                             psg->gsg_fvariance)) {
                                        bbackground = 1;
                        }

                        /* check if fits the current gaussian */
                        if (maha_dis < g_model_param.gmp_fthres_smd * 
                                psg->gsg_fvariance) {

                                bfits = 1; /* belongs to the current gaussian */

                                /* update the current gaussian distribution */
                                _gaussian_update(psg, r, g, b, 
                                                 maha_dis, fsingle_weight,
                                                 falpha);

                                /* sort
                                   all other weights are at the same place and 
                                   only the matched (imodes) is higher -> 
                                   just find the new place for it */
                                _gaussian_sort(psg, imodes);
                        }
                        else {
                                if (fsingle_weight < -falpha * 
                                        g_model_param.gmp_fct) {
                                        fsingle_weight = 0.f;
                                        (*pngaussians_used)--;
                                }
                                psg->gsg_fweight = fsingle_weight;
                        }
                }
                ftotal_weight += psg->gsg_fweight;
        }

        if (!bfits) {
                _generate_new_gaussian(psg, pngaussians_used, r, g, b, falpha);
        }

        // renormalize the weights
        ftotal_weight = 0.f;
        for (imodes = 0; imodes < (*pngaussians_used); ++imodes)
                ftotal_weight += psg[imodes].gsg_fweight;
        for (imodes = 0; imodes < (*pngaussians_used); ++imodes)
                psg[imodes].gsg_fweight /= ftotal_weight;

        return bbackground;                
}

gaussmix_image_t* gauss_mixture_create_image(int width, int height, 
                                             int channels)
{
        gaussmix_image_t *image = (gaussmix_image_t*)malloc(
                sizeof(gaussmix_image_t));
        if (!image) return NULL;

        image->gi_ucdata = (unsigned char*)calloc(width * height, 
                sizeof(unsigned char) * channels);
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
                             float **bg_model,
                             unsigned char **bg_model_used)
{
        int nmax_gaussians;

        _set_model_param(img_width, img_height);

        nmax_gaussians = g_model_param.gmp_inmodels;

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

void gauss_mixture_final(float **bg_model, 
                         unsigned char **bg_model_used)
{
        if (*bg_model) {
                free(*bg_model);
                *bg_model = NULL;
        }

        if (*bg_model_used) {
                free(*bg_model_used);
                *bg_model_used = NULL;
        }
}

void gauss_mixture_update(gaussmix_image_t *image, 
                          gaussmix_image_t *fg_mask, 
                          float *bg_model, 
                          unsigned char *bg_model_used)
{
        static long long nframe = 1; /* how many frames have been processed */

        int width = image->gi_iwidth;
        int height = image->gi_iheight;
        unsigned char *img_data = image->gi_ucdata;
        unsigned char *mask_data = fg_mask->gi_ucdata;
        gaussmix_single_gaussian_t *psg = (gaussmix_single_gaussian_t*)bg_model;
        int i = 0;
        unsigned char bbackground = 0;
        float falpha = 0.f;

        /* at the start, use faster learning speed */
        if (nframe++ < GAUSSMIX_WINDOW_SIZE / 2) {
                falpha = 1.0f / (2 * nframe);
        }
        else {
                falpha = g_model_param.gmp_falpha;
        }

        printf("nfram is: %d.\n", nframe);
        memset(fg_mask->gi_ucdata, 0, sizeof(unsigned char) * width * height);
        while (i++ < width * height) {
                bbackground = _update( (float)img_data[0], 
                                       (float)img_data[1], 
                                       (float)img_data[2],
                                       falpha,
                                       psg, 
                                       bg_model_used );
                if (bbackground) *mask_data = 255;

                img_data += 3; /* RGB 3 channels */
                mask_data++;
                psg += g_model_param.gmp_inmodels;
                bg_model_used++;
        }
}