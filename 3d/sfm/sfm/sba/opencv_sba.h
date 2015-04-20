//
//  opencv_sba.h
//  sfm
//
//  Created by  刘骥 on 15/4/16.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#ifndef __sfm__opencv_sba__
#define __sfm__opencv_sba__

#include <stdio.h>
#include <iostream>
#include <vector>
#include "opencv2/opencv.hpp"
#include "CloudPoint.h"
using namespace cv;
using namespace std;
void opencv_twoview_sba(const Mat_<double>&k,const Mat_<double>&p1,const Mat_<double>&p2,const Mat_<double>&points,const Mat_<double>&points1,const Mat_<double>&points2,Mat_<double>&new_p1,Mat_<double>&new_p2,Mat_<double>&new_points,int nconstpts3D=0,int nconstframes=1,int maxiter=150,int verbose=0);
void nviewSba(const Mat_<double>&cameraMatrix, vector<Mat_<double> >& pmatrices,vector<CloudPoint> &cloudPoints,int nframes,int nconstpts3D=0,int nconstframes=0,int maxiter=150,int verbose=0);
#endif /* defined(__sfm__opencv_sba__) */
