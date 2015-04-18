//
//  main.cpp
//  sbademo
//
//  Created by  刘骥 on 15/4/15.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#include <iostream>
#include "readparams.h"
#include "imgproj.h"
#include <math.h>
#include "sba.h"

#define MAXITER         100
#define MAXITER2        150
template <class T>
void print_array(T first,T last) {
    int i=0;
    while(first!=last)
    {
        std::cout<<i<<":"<<*first<<std::endl;
        ++first;
        ++i;
    }
}
struct globs_{
    double *rot0params; /* initial rotation parameters, combined with a local rotation parameterization */
    double *intrcalib; /* the 5 intrinsic calibration parameters in the order [fu, u0, v0, ar, skew],
                        * where ar is the aspect ratio fv/fu.
                        * Used only when calibration is fixed for all cameras;
                        * otherwise, it is null and the intrinsic parameters are
                        * included in the set of motion parameters for each camera
                        */
    int nccalib; /* number of calibration parameters that must be kept constant.
                  * 0: all parameters are free
                  * 1: skew is fixed to its initial value, all other parameters vary (i.e. fu, u0, v0, ar)
                  * 2: skew and aspect ratio are fixed to their initial values, all other parameters vary (i.e. fu, u0, v0)
                  * 3: meaningless
                  * 4: skew, aspect ratio and principal point are fixed to their initial values, only the focal length varies (i.e. fu)
                  * 5: all intrinsics are kept fixed to their initial values
                  * >5: meaningless
                  * Used only when calibration varies among cameras
                  */
    
    int ncdist; /* number of distortion parameters in Bouguet's model that must be kept constant.
                 * 0: all parameters are free
                 * 1: 6th order radial distortion term (kc[4]) is fixed
                 * 2: 6th order radial distortion and one of the tangential distortion terms (kc[3]) are fixed
                 * 3: 6th order radial distortion and both tangential distortion terms (kc[3], kc[2]) are fixed [i.e., only 2nd & 4th order radial dist.]
                 * 4: 4th & 6th order radial distortion terms and both tangential distortion ones are fixed [i.e., only 2nd order radial dist.]
                 * 5: all distortion parameters are kept fixed to their initial values
                 * >5: meaningless
                 * Used only when calibration varies among cameras and distortion is to be estimated
                 */
    int cnp, pnp, mnp; /* dimensions */
    
    double *ptparams; /* needed only when bundle adjusting for camera parameters only */
    double *camparams; /* needed only when bundle adjusting for structure parameters only */
} ;
/* unit quaternion from vector part */
#define _MK_QUAT_FRM_VEC(q, v){                                     \
(q)[1]=(v)[0]; (q)[2]=(v)[1]; (q)[3]=(v)[2];                      \
(q)[0]=sqrt(1.0 - (q)[1]*(q)[1] - (q)[2]*(q)[2]- (q)[3]*(q)[3]);  \
}
inline static void quatMultFast(double q1[FULLQUATSZ], double q2[FULLQUATSZ], double p[FULLQUATSZ])
{
    double t1, t2, t3, t4, t5, t6, t7, t8, t9;
    //double t10, t11, t12;
    
    t1=(q1[0]+q1[1])*(q2[0]+q2[1]);
    t2=(q1[3]-q1[2])*(q2[2]-q2[3]);
    t3=(q1[1]-q1[0])*(q2[2]+q2[3]);
    t4=(q1[2]+q1[3])*(q2[1]-q2[0]);
    t5=(q1[1]+q1[3])*(q2[1]+q2[2]);
    t6=(q1[1]-q1[3])*(q2[1]-q2[2]);
    t7=(q1[0]+q1[2])*(q2[0]-q2[3]);
    t8=(q1[0]-q1[2])*(q2[0]+q2[3]);
    
#if 0
    t9 =t5+t6;
    t10=t7+t8;
    t11=t5-t6;
    t12=t7-t8;
    
    p[0]= t2 + 0.5*(-t9+t10);
    p[1]= t1 - 0.5*(t9+t10);
    p[2]=-t3 + 0.5*(t11+t12);
    p[3]=-t4 + 0.5*(t11-t12);
#endif
    
    /* following fragment it equivalent to the one above */
    t9=0.5*(t5-t6+t7+t8);
    p[0]= t2 + t9-t5;
    p[1]= t1 - t9-t6;
    p[2]=-t3 + t9-t8;
    p[3]=-t4 + t9-t7;
}

static void img_projsRTS_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *hx, void *adata)
{
    int i, j;
    int cnp, pnp, mnp;
    double *pa, *pb, *pqr, *pt, *ppt, *pmeas, *Kparms, *pr0, lrot[FULLQUATSZ], trot[FULLQUATSZ];
    //int n;
    int m, nnz;
    struct globs_ *gl;
    
    gl=(struct globs_ *)adata;
    cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
    Kparms=gl->intrcalib;
    
    //n=idxij->nr;
    m=idxij->nc;
    pa=p; pb=p+m*cnp;
    
    for(j=0; j<m; ++j){
        /* j-th camera parameters */
        pqr=pa+j*cnp;
        pt=pqr+3; // quaternion vector part has 3 elements
        pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
        _MK_QUAT_FRM_VEC(lrot, pqr);
        quatMultFast(lrot, pr0, trot); // trot=lrot*pr0
        
        nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */
        
        for(i=0; i<nnz; ++i){
            ppt=pb + rcsubs[i]*pnp;
            pmeas=hx + idxij->val[rcidxs[i]]*mnp; // set pmeas to point to hx_ij
            
            calcImgProjFullR(Kparms, trot, pt, ppt, pmeas); // evaluate Q in pmeas
            std::cout<<trot[0]<<" "<<trot[1]<<" "<<trot[2]<<" "<<trot[3]<<std::endl;
            std::cout<<Kparms[0]<<" "<<Kparms[1]<<" "<<Kparms[2]<<" "<<Kparms[3]<<" "<<Kparms[4]<<std::endl;
            std::cout<<pmeas[0]<<"\t"<<pmeas[1]<<std::endl;
        }
    }
}
static void img_projsRTS_jac_x(double *p, struct sba_crsm *idxij, int *rcidxs, int *rcsubs, double *jac, void *adata)
{
    int i, j;
    int cnp, pnp, mnp;
    double *pa, *pb, *pqr, *pt, *ppt, *pA, *pB, *Kparms, *pr0;
    //int n;
    int m, nnz, Asz, Bsz, ABsz;
    struct globs_ *gl;
    
    gl=(struct globs_ *)adata;
    cnp=gl->cnp; pnp=gl->pnp; mnp=gl->mnp;
    Kparms=gl->intrcalib;
    
    //n=idxij->nr;
    m=idxij->nc;
    pa=p; pb=p+m*cnp;
    Asz=mnp*cnp; Bsz=mnp*pnp; ABsz=Asz+Bsz;
    
    for(j=0; j<m; ++j){
        /* j-th camera parameters */
        pqr=pa+j*cnp;
        pt=pqr+3; // quaternion vector part has 3 elements
        pr0=gl->rot0params+j*FULLQUATSZ; // full quat for initial rotation estimate
        
        nnz=sba_crsm_col_elmidxs(idxij, j, rcidxs, rcsubs); /* find nonzero hx_ij, i=0...n-1 */
        
        for(i=0; i<nnz; ++i){
            ppt=pb + rcsubs[i]*pnp;
            pA=jac + idxij->val[rcidxs[i]]*ABsz; // set pA to point to A_ij
            pB=pA  + Asz; // set pB to point to B_ij
            
            calcImgProjJacRTS(Kparms, pr0, pqr, pt, ppt, (double (*)[6])pA, (double (*)[3])pB); // evaluate dQ/da, dQ/db in pA, pB
        }
    }
}
/* convert a vector of camera parameters so that rotation is represented by
 * the vector part of the input quaternion. The function converts the
 * input quaternion into a unit one with a non-negative scalar part. Remaining
 * parameters are left unchanged.
 *
 * Input parameter layout: intrinsics (5, optional), distortion (5, optional), rot. quaternion (4), translation (3)
 * Output parameter layout: intrinsics (5, optional), distortion (5, optional), rot. quaternion vector part (3), translation (3)
 */
void quat2vec(double *inp, int nin, double *outp, int nout)
{
    double mag, sg;
    int i;
    
    /* intrinsics & distortion */
    if(nin>7) // are they present?
        for(i=0; i<nin-7; ++i)
            outp[i]=inp[i];
    else
        i=0;
    
    /* rotation */
    /* normalize and ensure that the quaternion's scalar component is non-negative;
     * if not, negate the quaternion since two quaternions q and -q represent the
     * same rotation
     */
    mag=sqrt(inp[i]*inp[i] + inp[i+1]*inp[i+1] + inp[i+2]*inp[i+2] + inp[i+3]*inp[i+3]);
    sg=(inp[i]>=0.0)? 1.0 : -1.0;
    mag=sg/mag;
    outp[i]  =inp[i+1]*mag;
    outp[i+1]=inp[i+2]*mag;
    outp[i+2]=inp[i+3]*mag;
    i+=3;
    
    /* translation*/
    for( ; i<nout; ++i)
        outp[i]=inp[i+1];
}

/* convert a vector of camera parameters so that rotation is represented by
 * a full unit quaternion instead of its input 3-vector part. Remaining
 * parameters are left unchanged.
 *
 * Input parameter layout: intrinsics (5, optional), distortion (5, optional), rot. quaternion vector part (3), translation (3)
 * Output parameter layout: intrinsics (5, optional), distortion (5, optional), rot. quaternion (4), translation (3)
 */
void vec2quat(double *inp, int nin, double *outp, int nout)
{
    double *v, q[FULLQUATSZ];
    int i;
    
    /* intrinsics & distortion */
    if(nin>7-1) // are they present?
        for(i=0; i<nin-(7-1); ++i)
            outp[i]=inp[i];
    else
        i=0;
    
    /* rotation */
    /* recover the quaternion from the vector */
    v=inp+i;
    _MK_QUAT_FRM_VEC(q, v);
    outp[i]  =q[0];
    outp[i+1]=q[1];
    outp[i+2]=q[2];
    outp[i+3]=q[3];
    i+=FULLQUATSZ;
    
    /* translation */
    for( ; i<nout; ++i)
        outp[i]=inp[i-1];
}

int main(int argc, const char * argv[]) {

    struct globs_ globs;
    int verbose=10;
    double K[9], ical[5];
    double opts[SBA_OPTSSZ], info[SBA_INFOSZ];
    int cnp=6;
    int pnp=3;
    int mnp=2;
    int numpts3D=15;
    int nconstpts3D=0;
    int nframes=2;
    int nconstframes=0;
    char *vmask=new char[numpts3D*nframes];
    std::fill_n(vmask, numpts3D*nframes, 1);
    double *motstruct=new double[nframes*cnp+numpts3D*pnp]{0.000000,0.000000,0.000000,0.000000,0.000000,0.000000,0.040564,-0.042632,0.032459,-0.401172,-0.109200,0.909470,1.034947 -0.300480,4.761752};
    int numprojs=numpts3D*nframes;
    double *imgpts=new double[numpts3D*nframes*mnp]{686.495911,309.032562,538.390503,271.069916};
    double *covimgpts=NULL;
    /* set up globs structure */
    globs.cnp=cnp;
    globs.pnp=pnp;
    globs.mnp=mnp;
    globs.rot0params=new double[FULLQUATSZ*nframes]{1.000000,0.000000,0.000000,0.000000,0.997739,0.040564,-0.042632,0.032459};
    ical[0]=864.5; // fu
    ical[1]=499.5; // u0
    ical[2]=374.5; // v0
    ical[3]=1.0; // ar
    ical[4]=0; // s
    globs.intrcalib=ical;
    globs.ptparams=NULL;
    globs.camparams=NULL;
    
    /* call sparse LM routine */
    opts[0]=SBA_INIT_MU; opts[1]=SBA_STOP_THRESH; opts[2]=SBA_STOP_THRESH;
    opts[3]=SBA_STOP_THRESH;
    //opts[3]=0.05*numprojs; // uncomment to force termination if the average reprojection error drops below 0.05
    opts[4]=0.0;
    //opts[4]=1E-05; // uncomment to force termination if the relative reduction in the RMS reprojection error drops below 1E-05
    int n=sba_motstr_levmar_x(numpts3D, nconstpts3D, nframes, nconstframes, vmask, motstruct, cnp, pnp, imgpts, covimgpts, mnp,img_projsRTS_x,img_projsRTS_jac_x,(void*)(&globs),MAXITER2, verbose, opts, info);
    printf("%d\n",n);
 //   printSBAData(stdout, motstruct, cnp, pnp, mnp, vec2quat, cnp+1, nframes, numpts3D, imgpts, numprojs, vmask);
    return 0;
}
