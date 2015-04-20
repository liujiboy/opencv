//
//  opencv_sba.cpp
//  sfm
//
//  Created by  刘骥 on 15/4/16.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#include "opencv_sba.h"
#include "readparams.h"
#include "imgproj.h"
#include "sba.h"
#include <math.h>
#define MAXITER         100
#define MAXITER2        150
/* unit quaternion from vector part */
#define _MK_QUAT_FRM_VEC(q, v){                                     \
(q)[1]=(v)[0]; (q)[2]=(v)[1]; (q)[3]=(v)[2];                      \
(q)[0]=sqrt(1.0 - (q)[1]*(q)[1] - (q)[2]*(q)[2]- (q)[3]*(q)[3]);  \
}
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
/*** MEASUREMENT VECTOR AND JACOBIAN COMPUTATION FOR THE EXPERT DRIVERS ***/

/* FULL BUNDLE ADJUSTMENT, I.E. SIMULTANEOUS ESTIMATION OF CAMERA AND STRUCTURE PARAMETERS */

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * hx the prediction of the measurements, i.e. the projections of 3D points in the m images. The measurements
 * are returned in the order (hx_11^T, .. hx_1m^T, ..., hx_n1^T, .. hx_nm^T)^T, where hx_ij is the predicted
 * projection of the i-th point on the j-th camera.
 * Notice that depending on idxij, some of the hx_ij might be missing
 *
 */
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
            //calcImgProj(Kparms, pr0, pqr, pt, ppt, pmeas); // evaluate Q in pmeas
        }
    }
}

/* Given a parameter vector p made up of the 3D coordinates of n points and the parameters of m cameras, compute in
 * jac the jacobian of the predicted measurements, i.e. the jacobian of the projections of 3D points in the m images.
 * The jacobian is returned in the order (A_11, ..., A_1m, ..., A_n1, ..., A_nm, B_11, ..., B_1m, ..., B_n1, ..., B_nm),
 * where A_ij=dx_ij/db_j and B_ij=dx_ij/db_i (see HZ).
 * Notice that depending on idxij, some of the A_ij, B_ij might be missing
 *
 */
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
//r表示四元数的后三位
void rotmat2quat(double *R,double *r){
    /*
    double q[4];
    double tmp[4];
    int maxpos=0;
    tmp[0]=1.0 + R[0] + R[4] + R[8];
    tmp[1]=1.0 + R[0] - R[4] - R[8];
    tmp[2]=1.0 - R[0] + R[4] - R[8];
    tmp[3]=1.0 - R[0] - R[4] + R[8];
    double mag=-1.0;
    for(int i=0; i<4; i++){
        if(tmp[i]>mag){
            mag=tmp[i];
            maxpos=i;
        }
    }
    
    if(maxpos==0){
        q[0]=sqrt(tmp[0])*0.5;
        q[1]=(R[7] - R[5])/(4.0*q[0]);
        q[2]=(R[2] - R[6])/(4.0*q[0]);
        q[3]=(R[3] - R[1])/(4.0*q[0]);
    }
    else if(maxpos==1){
        q[1]=sqrt(tmp[1])*0.5;
        q[0]=(R[7] - R[5])/(4.0*q[1]);
        q[2]=(R[3] + R[1])/(4.0*q[1]);
        q[3]=(R[2] + R[6])/(4.0*q[1]);
    }
    else if(maxpos==2){
        q[2]=sqrt(tmp[2])*0.5;
        q[0]=(R[2] - R[6])/(4.0*q[2]);
        q[1]=(R[3] + R[1])/(4.0*q[2]);
        q[3]=(R[7] + R[5])/(4.0*q[2]);
    }
    else if(maxpos==3){
        q[3]=sqrt(tmp[3])*0.5;
        q[0]=(R[3] - R[1])/(4.0*q[3]);
        q[1]=(R[2] + R[6])/(4.0*q[3]);
        q[2]=(R[7] + R[5])/(4.0*q[3]);
    }
    else{
    }
    
    // enforce unit length
    mag=q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3];
    if(mag!=1.0)
    {
        mag=1.0/sqrt(mag);
        q[0]*=mag; q[1]*=mag; q[2]*=mag; q[3]*=mag;
     
    }*/
    double q[4];
    q[0]=sqrt(1.0 + R[0] + R[4] + R[8])*0.5;
    q[1]=(R[7] - R[5])/(4.0*q[0]);
    q[2]=(R[2] - R[6])/(4.0*q[0]);
    q[3]=(R[3] - R[1])/(4.0*q[0]);

    r[0]=q[0];
    r[1]=q[1];
    r[2]=q[2];
    r[3]=q[3];
}
//r表示四元数的后三位
void quat2rotmat(double *r,double *R)
{   double q1=r[0];
    double q2=r[1];
    double q3=r[2];
    double q0=sqrt(1-q1*q1-q2*q2-q3*q3);
    R[0]=q0*q0+q1*q1-q2*q2-q3*q3;
    R[1]=2*(q1*q2-q0*q3);
    R[2]=2*(q1*q3+q0*q2);
    
    R[3]=2*(q1*q2+q0*q3);
    R[4]=q0*q0+q2*q2-q1*q1-q3*q3;
    R[5]=2*(q2*q3-q0*q1);
    
    R[6]=2*(q1*q3-q0*q2);
    R[7]=2*(q2*q3+q0*q1);
    R[8]=q0*q0+q3*q3-q1*q1-q2*q2;
}
void cammat2quat(const Mat_<double>&p,double*r,double*t)
{
    double R[9]={p(0,0),p(0,1),p(0,2),p(1,0),p(1,1),p(1,2),p(2,0),p(2,1),p(2,2)};
    rotmat2quat(R, r);
    t[0]=p(0,3);
    t[1]=p(1,3);
    t[2]=p(2,3);
}
void quat2cammat(double*r,double*t,Mat_<double>&p)
{
    double R[9];
    quat2rotmat(r,R);
    p(0,0)=R[0];
    p(0,1)=R[1];
    p(0,2)=R[2];
    p(1,0)=R[3];
    p(1,1)=R[4];
    p(1,2)=R[5];
    p(2,0)=R[6];
    p(2,1)=R[7];
    p(2,2)=R[8];
    p(0,3)=t[0];
    p(1,3)=t[1];
    p(2,3)=t[2];
    
}

void opencv_twoview_sba(const Mat_<double>&k,const Mat_<double>&p1,const Mat_<double>&p2,const Mat_<double>&points,const Mat_<double>&points1,const Mat_<double>&points2,Mat_<double>&new_p1,Mat_<double>&new_p2,Mat_<double>&new_points,int nconstpts3D,int nconstframes,int maxiter,int verbose)
{
    struct globs_ globs;
    double ical[5]={k(0,0),k(0,2),k(1,2),k(1,1)/k(0,0),0};//f cx cy ar s;
    double opts[SBA_OPTSSZ], info[SBA_INFOSZ];
    int cnp=6;  //相机矩阵用3位表示旋转，3位表示位移
    int pnp=3;  //三维点的坐标数
    int mnp=2;  //二维点的坐标数
    int numpts3D=points.cols;   //三维点的数量
    int nframes=2;  //相机的数量
    int numprojs=numpts3D*nframes; //在所有相机下，三维点共计有多少个二维投影
    //vmask[i,j]表示第i个点在第j个镜头下是否可见，此处填充为全1，因为点在两个镜头下均可见
    char *vmask=new char[numpts3D*nframes];
    std::fill_n(vmask, numpts3D*nframes, 1);
    //将相机矩阵p1和p2转换为四元数表示
    //r1和r2是四元数，表示旋转矩阵
    //t1和t2是位移向量
    double r1[4],t1[3],r2[4],t2[3];
    cammat2quat(p1,r1,t1);
    cammat2quat(p2,r2,t2);
    //motstruct是待优化的相机矩阵和三维点，其结构为(r1,t1,r2,t2,X[1],X[2]...X[n])
    int motstruct_size=nframes*cnp+numpts3D*pnp;
    double *motstruct=new double[motstruct_size]();
    //拷贝相机矩阵
    std::copy(r1+1, r1+4, motstruct);
    std::copy(t1, t1+3, motstruct+3);
    std::copy(r2+1, r2+4, motstruct+6);
    std::copy(t2, t2+3, motstruct+9);
    //拷贝三维点
    int pstart=nframes*cnp; //三维点的开始位置
    for(int i=0;i<numpts3D;i++)
    {
        motstruct[pstart+i*pnp]=points(0,i);
        motstruct[pstart+i*pnp+1]=points(1,i);
        motstruct[pstart+i*pnp+2]=points(2,i);
    }
    //如果要对相机旋转矩阵和三维点的位置同时优化，必须将相机矩阵的旋转初始化为0，即四元数表示的(1,0,0,0)
    //并用globs.rot0params保存了相机旋转矩阵
    //若只对三维点的位置进行优化，此步不做
    for(int i=0; i<nframes; ++i){
        int j=(i+1)*cnp; // 跳过位移向量
        motstruct[j-4]=motstruct[j-5]=motstruct[j-6]=0.0; // 设置为(1,0,0,0)
    }
    //imgpts保存三维点在每个相机下的投影，即二维点
    //imgpts[i,j]是第i个三维点在第j个相机下的投影
    double *imgpts=new double[numprojs*mnp]();
    for(int i=0;i<numpts3D;i++)
    {
        imgpts[i*nframes*mnp]=points1(0,i);
        imgpts[i*nframes*mnp+1]=points1(1,i);
        imgpts[i*nframes*mnp+2]=points2(0,i);
        imgpts[i*nframes*mnp+3]=points2(1,i);
    }
    double *covimgpts=NULL;
    //设置globs
    globs.cnp=cnp;
    globs.pnp=pnp;
    globs.mnp=mnp;
    globs.rot0params=new double[FULLQUATSZ*nframes]();
    std::copy(r1,r1+4, globs.rot0params);
    std::copy(r2,r2+4, globs.rot0params+4);
    globs.intrcalib=ical;
    globs.ptparams=NULL;
    globs.camparams=NULL;
    
    //设置优化选项
    opts[0]=SBA_INIT_MU; opts[1]=SBA_STOP_THRESH; opts[2]=SBA_STOP_THRESH;
    opts[3]=SBA_STOP_THRESH;
    //opts[3]=0.05*numprojs; // uncomment to force termination if the average reprojection error drops below 0.05
    opts[4]=0.0;
    //opts[4]=1E-05; // uncomment to force termination if the relative reduction in the RMS reprojection error drops below 1E-05
     
    /*opts[0]=1E-8; opts[1]=1E-15; opts[2]=1E-15;
    opts[3]=1E-15;
    //opts[3]=0.05*numprojs; // uncomment to force termination if the average reprojection error drops below 0.05
    opts[4]=0.0;*/
    //优化
    int n=sba_motstr_levmar_x(numpts3D, nconstpts3D, nframes, nconstframes, vmask, motstruct, cnp, pnp, imgpts, covimgpts, mnp,img_projsRTS_x,img_projsRTS_jac_x,(void*)(&globs),maxiter, verbose, opts, info);
    if(n!=SBA_ERROR)
    {
        
        /* combine the local rotation estimates with the initial ones */
        for(int i=0; i<nframes; ++i){
            double *v, qs[FULLQUATSZ], *q0, prd[FULLQUATSZ];
            
            /* retrieve the vector part */
            v=motstruct + (i+1)*cnp - 6; // note the +1, we access the motion parameters from the right, assuming 3 for translation!
            _MK_QUAT_FRM_VEC(qs, v);
            
            q0=globs.rot0params+i*FULLQUATSZ;
            quatMultFast(qs, q0, prd); // prd=qs*q0
            
            /* copy back vector part making sure that the scalar part is non-negative */
            if(prd[0]>=0.0){
                v[0]=prd[1];
                v[1]=prd[2];
                v[2]=prd[3];
            }
            else{ // negate since two quaternions q and -q represent the same rotation
                v[0]=-prd[1];
                v[1]=-prd[2];
                v[2]=-prd[3];
            }
        }
        //printSBAData(stdout, motstruct, cnp, pnp, mnp, vec2quat, cnp+1, nframes, numpts3D, imgpts, numprojs, vmask);

        new_p1.create(3, 4);
        new_p2.create(3,4);
        new_points.create(4,numpts3D);
        std::copy(motstruct, motstruct+3, r1);
        std::copy(motstruct+3, motstruct+6, t1);
        std::copy(motstruct+6, motstruct+9, r2);
        std::copy(motstruct+9, motstruct+12, t2);
        quat2cammat(r1, t1, new_p1);
        quat2cammat(r2, t2, new_p2);
        int i=0;
        for (double *p=motstruct+12; p!=motstruct+motstruct_size; p+=3) {
            new_points(0,i)=*p;
            new_points(1,i)=*(p+1);
            new_points(2,i)=*(p+2);
            new_points(3,i)=1;
            i++;
        }
               //std::cout<<new_points<<std::endl;
    }

}
void nviewSba(const Mat_<double>&cameraMatrix, vector<Mat_<double> >& pmatrices,vector<CloudPoint> &cloudPoints,int nframes,int nconstpts3D,int nconstframes,int maxiter,int verbose)
{
    int cnp=6;  //相机矩阵用3位表示旋转，3位表示位移
    int pnp=3;  //三维点的坐标数
    int mnp=2;  //二维点的坐标数
    double ical[5]={cameraMatrix(0,0),cameraMatrix(0,2),cameraMatrix(1,2),cameraMatrix(1,1)/cameraMatrix(0,0),0};//f cx cy ar s;
    double opts[SBA_OPTSSZ], info[SBA_INFOSZ];
    struct globs_ globs;
    //设置globs
    globs.cnp=cnp;
    globs.pnp=pnp;
    globs.mnp=mnp;
    globs.rot0params=new double[FULLQUATSZ*nframes]();
    globs.intrcalib=ical;
    globs.ptparams=NULL;
    globs.camparams=NULL;
    
    //设置优化选项
    opts[0]=SBA_INIT_MU;
    opts[1]=SBA_STOP_THRESH;
    opts[2]=SBA_STOP_THRESH;
    opts[3]=SBA_STOP_THRESH;
    opts[4]=0.0;
    
    
    int numpts3D=(int)cloudPoints.size();   //三维点的数量
    int numprojs=0; //在所有相机下，三维点共计有多少个二维投影
    //vmask[i,j]表示第i个点在第j个镜头下是否可见，此处填充为全1，因为点在两个镜头下均可见
    char *vmask=new char[numpts3D*nframes]();
    for(int i=0;i<numpts3D;i++)
    {
        CloudPoint cp=cloudPoints[i];
        for (int j=0; j<nframes; j++) {
            int index=i*nframes+j;
            if(cp.getPointIndex(j)!=-1)
            {
                vmask[index]=1;
                numprojs++;
            }
        }
    }
    //motstruct是待优化的相机矩阵和三维点，其结构为(r1,t1,r2,t2,X[1],X[2]...X[n])
    int motstruct_size=nframes*cnp+numpts3D*pnp;
    double *motstruct=new double[motstruct_size]();
    for(int i=0;i<nframes;i++)
    {
        Mat_<double> p=pmatrices[i];
        double r[4],t[3];
        cammat2quat(p, r, t);
        copy(r+1, r+4, motstruct+i*cnp);
        copy(t,t+3,motstruct+i*cnp+3);
        copy(r,r+4, globs.rot0params+i*FULLQUATSZ);
    }
    //拷贝三维点
    int pstart=nframes*cnp; //三维点的开始位置
    for(int i=0;i<numpts3D;i++)
    {
        CloudPoint cp=cloudPoints[i];
        motstruct[pstart+i*pnp]=cp.x;
        motstruct[pstart+i*pnp+1]=cp.y;
        motstruct[pstart+i*pnp+2]=cp.z;
    }
    //如果要对相机旋转矩阵和三维点的位置同时优化，必须将相机矩阵的旋转初始化为0，即四元数表示的(1,0,0,0)
    //并用globs.rot0params保存了相机旋转矩阵
    //若只对三维点的位置进行优化，此步不做
    for(int i=0; i<nframes; ++i){
        int j=(i+1)*cnp; // 跳过位移向量
        motstruct[j-4]=motstruct[j-5]=motstruct[j-6]=0.0; // 设置为(1,0,0,0)
    }
    //imgpts保存三维点在每个相机下的投影，即二维点
    double *imgpts=new double[numprojs*mnp]();
    for(int i=0,n=0;i<numpts3D;i++)
    {
        CloudPoint cp=cloudPoints[i];
        for (int j=0; j<nframes; j++) {
            Point2d point;
            if(cp.getPointInFrame(j, point))
            {
                imgpts[n*mnp]=point.x;
                imgpts[n*mnp+1]=point.y;
                n++;
            }
        }
    }

    double *covimgpts=NULL;
    
   
    //优化
    int n=sba_motstr_levmar_x(numpts3D, nconstpts3D, nframes, nconstframes, vmask, motstruct, cnp, pnp, imgpts, covimgpts, mnp,img_projsRTS_x,img_projsRTS_jac_x,(void*)(&globs),maxiter, verbose, opts, info);
    if(n!=SBA_ERROR)
    {
        
        /* combine the local rotation estimates with the initial ones */
        for(int i=0; i<nframes; ++i){
            double *v, qs[FULLQUATSZ], *q0, prd[FULLQUATSZ];
            
            /* retrieve the vector part */
            v=motstruct + (i+1)*cnp - 6; // note the +1, we access the motion parameters from the right, assuming 3 for translation!
            _MK_QUAT_FRM_VEC(qs, v);
            
            q0=globs.rot0params+i*FULLQUATSZ;
            quatMultFast(qs, q0, prd); // prd=qs*q0
            
            /* copy back vector part making sure that the scalar part is non-negative */
            if(prd[0]>=0.0){
                v[0]=prd[1];
                v[1]=prd[2];
                v[2]=prd[3];
            }
            else{ // negate since two quaternions q and -q represent the same rotation
                v[0]=-prd[1];
                v[1]=-prd[2];
                v[2]=-prd[3];
            }
        }
        //printSBAData(stdout, motstruct, cnp, pnp, mnp, vec2quat, cnp+1, nframes, numpts3D, imgpts, numprojs, vmask);
        for(int i=0;i<nframes;i++)
        {
            Mat_<double> p(3,4);
            double r[3],t[3];
            copy(motstruct+i*cnp, motstruct+3+i*cnp, r);
            copy(motstruct+3+i*cnp, motstruct+6+i*cnp, t);
            quat2cammat(r,t,p);
            pmatrices[i]=p;
        }
        for (int i=0; i<numpts3D; i++) {
            CloudPoint &cp=cloudPoints[i];
            cp.x=motstruct[nframes*cnp+i*pnp];
            cp.y=motstruct[nframes*cnp+i*pnp+1];
            cp.z=motstruct[nframes*cnp+i*pnp+2];
        }
        /*
        new_p1.create(3, 4);
        new_p2.create(3,4);
        new_points.create(4,numpts3D);
        std::copy(motstruct, motstruct+3, r1);
        std::copy(motstruct+3, motstruct+6, t1);
        std::copy(motstruct+6, motstruct+9, r2);
        std::copy(motstruct+9, motstruct+12, t2);
        quat2cammat(r1, t1, new_p1);
        quat2cammat(r2, t2, new_p2);
        int i=0;
        for (double *p=motstruct+12; p!=motstruct+motstruct_size; p+=3) {
            new_points(0,i)=*p;
            new_points(1,i)=*(p+1);
            new_points(2,i)=*(p+2);
            new_points(3,i)=1;
            i++;
        }*/
        //std::cout<<new_points<<std::endl;
    }

}
