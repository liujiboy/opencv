//
//  SFM.h
//  sfm
//
//  Created by  刘骥 on 15/4/21.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#ifndef __sfm__SFM__
#define __sfm__SFM__

#include <iostream>
#include <opencv2/opencv.hpp>
#include "Cloud.h"
#include "readparams.h"
#include "imgproj.h"
#include "sba.h"
#include <math.h>
using namespace cv;
using namespace std;
class SFM{
private:
    const Mat &cameraMatrix;
    vector<Mat_<double> >pmatrices;
    const vector<vector<Point2d> > &keypoints; //keypoints[i][j]表示第i张第j个图片的特征点
    const vector<vector<Mat> > &fmatrices; //fmatrices[i][j]表示第i张图片和第j张图片间的基础矩阵
    const vector<vector<vector<DMatch>> > &matches; //matches[i][j]表示第i张图片和第j张图片间特征点的匹配
    int nframes;
    Cloud cloud;
    void getMatchPoints(const vector<Point2d>& keypoints1,const vector<Point2d>& keypoints2,const vector<DMatch>&matches,Mat_<double>&points1,Mat_<double>&points2);
    void computeP(const Mat&E,vector<Mat>&pVector);
    Mat_<double> trianglulateAndFindCameraMatrix(const Mat_<double>&points1,const Mat_<double>&points2,const Mat&P1,const Mat&F,Mat&P2);
    void triangulate(const Mat_<double>&points1,const Mat_<double>&points2,const Mat&P1,const Mat&P2,Mat_<double>&points);
    void triangulatePoint(const Mat_<double>&point1,const Mat_<double>&point2,const Mat_<double>&P1,const Mat_<double>&P2,Mat_<double>& point);
    Mat_<double> findPmatrixByKnownPoints(const vector<DMatch>&match,bool* known,const Mat_<double>&points,int frame);
    void reconstructByKnownPoints(int frame,int prevFrame,const vector<DMatch>&match,bool *known,Mat_<double> &points);
    void nviewSba(int frameNum,int nconstframes=1,int nconstpts3D=0,int maxiter=150,int verbose=0);
public:
    SFM(const Mat &cameraMatrix,const vector<vector<Point2d> > &keypoints,const vector<vector<Mat> > &fmatrices,const vector<vector<vector<DMatch>> > &matches,int nframes);
    void initialReconstruct();
    void addView(int frame);
    Cloud& getCloud();
};
#endif /* defined(__sfm__SFM__) */
