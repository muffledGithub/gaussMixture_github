/**
 * guass_mixture_param.h
 *
 * Author:
 * vanilla
 *
 * Abstract:
 * @1. predefined parameters for the Gaussian mixture model 
 *     background subtraction
 */

#ifndef _GAUSS_MIXTURE_PARAM_H_
#define _GAUSS_MIXTURE_PARAM_H_

#ifdef __cplusplus
extern "C" {
#endif

#define GAUSSMIX_DEFAULT_CHANNELS         3                 

/* default parameters of gaussian background detection algorithm */

/* lambda = 2.5 is 99% */
#define GAUSSMIX_STD_THRESHOLD            4.0f 

/* Learning rate; alpha = 1/CV_GBG_WINDOW_SIZE */
#define GAUSSMIX_WINDOW_SIZE              500 

/* threshold sum of weights for background test */
#define GAUSSMIX_BACKGROUND_THRESHOLD     0.9f   

/* lambda = 2.5 is 99% */
#define GAUSSMIX_STD_THRESHOLD_GENERATE   3.0f  

/* = K = number of Gaussians in mixture */
#define GAUSSMIX_NGAUSSIANS               4 

/* initial variance for new components*/
#define GAUSSMIX_VAR_INIT                 10.0f   

#define GAUSSMIX_VAR_MIN                  4.0f

#define GAUSSMIX_VAR_MAX                  5 * GAUSSMIX_VAR_INIT

/* for postfiltering */
#define GAUSSMIX_MINAREA                  15.0f    


/* additional parameters */

/* complexity reduction prior constant 0 - no reduction 
   of number of components*/
#define GAUSSMIX_CT                       0.01f   

/* value to use in the segmentation mask for shadows, 
   sot 0 not to do shadow detection*/
#define GAUSSMIX_SHADOW_VALUE             127   

/* Tau - shadow threshold, see the paper for explanation*/
#define GAUSSMIX_SHADOW_TAU               0.5f   

typedef struct gaussmix_model_param_s
{
        /* image info */
        int gmp_iwidth;
        int gmp_iheight;
        int gmp_ichannels;

        /* defult 1 - do postfiltering 
           which will make shadow detection results also give value 255 */
        unsigned char gmp_bpost_filtering;
        double gmp_dmin_area;  /* for postfiltering */

        unsigned char gmp_binit; /* default 1, faster updates at start */

        /* very important parameters - things you will change */

        /* alpha - speed of update
           If the time interval you want to average over is T set alpha = 1/T. 
           It is also usefull at start to make T slowly increase from 1 until 
           the desired T. */
        float gmp_falpha;

        /* cthr - threshold on the squared Mahalan. dist. to decide 
                  if it is well described by the background model or not. 
                  Related to Cthr from the paper.
           This does not influence the update of the background. 
           A typical value could be 4 sigma and that is cthr = 4 * 4 = 16 */
        float gmp_fcthr;

        /* less important parameters 
           - things you might change but be carefull */

        /* thres_smd - threshold on the squared Mahalan. dist. to decide
                       whether a sample is close to the existing components.
                       If it is not closeto any a new component will be 
                       generated.
           I use 3 sigma => thres_smd = 3 * 3 = 9.
           Smaller thres_smd leads to more generated components and 
           higher thres_smd might lead to small number of components 
           but they can grow too large. */
        float gmp_fthres_smd;

        /* 1 - cf from the paper
           one_minus_cf - threshold when the component becomes significant 
                          enough to be included into the background model. 
            So I use cf = 0.1 => one_minus_cf = 0.9
            For alpha = 0.001 it means that the mode should exist for 
            approximately 105 frames before it is considered foreground */
        float gmp_fone_minus_cf;

        /* Initial standard deviation for the newly generated components.
           It will influence the speed of adaptation, and a good guess should 
              be made.
           A simple way is to estimate the typical standard deviation from 
              the images.
           I used here 10 as a reasonable value. */
        float gmp_fvar_init;

        /* the variance of every Gaussian model should be limited in
           [var_min, var_max] */
        float gmp_fvar_max, gmp_fvar_min;

        /* ct - complexity reduction prior
           This is related to the number of samples needed to accept that a 
           component actually exists.We use ct = 0.05 of all the samples. 
           By setting ct = 0 you getthe standard Stauffer&Grimson algorithm 
           (maybe not exact but very similar). */
        float gmp_fct;

        /* even less important parameters */

        /* max number of modes - const - 4 is usually enough */
        int gmp_inmodels;

        /* shadow detection parameters */

        /* default 1 - do shadow detection */
        unsigned char gmp_bshadow_detection;

        /* do shadow detection - insert this value as the detection result */
        unsigned char gmp_ucshadow_value;

        /* tau - shadow threshold. 
           The shadow is detected if the pixel is darker version of the 
           background. 
           tau is a threshold on how much darker the shadow can be.
           tau = 0.5 means that if pixel is more than 2 times darker then it is 
           not shadow.
           See: Prati,Mikic, Trivedi, Cucchiarra, "Detecting Moving Shadows...",
           IEEE PAMI, 2003. */
        float gmp_ftau;     
} gaussmix_model_param_t;

#ifdef __cplusplus
}
#endif

#endif /* _GAUSS_MIXTURE_PARAM_H_ */
