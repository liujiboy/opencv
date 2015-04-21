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
//Structure From Motion
class SFM{
private:
    //相机矩阵
    const Mat &cameraMatrix;
    //相机投影矩阵，pmatrices[i]表示视图i的投影矩阵
    vector<Mat_<double> >pmatrices;
    //keypoints[i][j]表示视图i的第j个特征点
    const vector<vector<Point2d> > &keypoints;
    //fmatrices[i][j]表示视图i和视图j间的基础矩阵
    const vector<vector<Mat> > &fmatrices;
    //matches[i][j]表示视图i和视图j的特征点匹配
    const vector<vector<vector<DMatch>> > &matches;
    //视图数量
    int nframes;
    //重建的点云
    Cloud cloud;
    //根据匹配取特征点
    //keypoints1：视图1的特征点
    //keypoints2：视图2的特征点
    //matches：视图1和视图2的匹配
    //points1，points2：返回的视图1和视图2匹配的特征点
    void getMatchPoints(const vector<Point2d>& keypoints1,const vector<Point2d>& keypoints2,const vector<DMatch>&matches,Mat_<double>&points1,Mat_<double>&points2);
    //根据本质矩阵计算出4个投影矩阵
    //E：本质矩阵
    //pVector：返回的4个投影矩阵
    void computeP(const Mat&E,vector<Mat>&pVector);
    //对视图1和视图2进行三角化，并返回视图2的投影矩阵
    //在初始重建中假设P1为已知，P2的可能有4种
    //points1,points2：视图1和视图2的特征点
    //P1：视图1的投影矩阵
    //F：基础矩阵
    //P2：返回的视图2的投影矩阵
    Mat_<double> trianglulateAndFindCameraMatrix(const Mat_<double>&points1,const Mat_<double>&points2,const Mat&P1,const Mat&F,Mat&P2);
    //视图1和视图2的特征点进行三角化
    //points1,points2：视图1和视图2的特征点
    //P1,P2：视图1和视图2的投影矩阵
    //points：重建的三维点
    void triangulate(const Mat_<double>&points1,const Mat_<double>&points2,const Mat&P1,const Mat&P2,Mat_<double>&points);
    //一对特征点进行三角化
    //point1，point2：一对特征点
    //P1,P2：视图1和视图2的投影矩阵
    //point：重建的三维点
    void triangulatePoint(const Mat_<double>&point1,const Mat_<double>&point2,const Mat_<double>&P1,const Mat_<double>&P2,Mat_<double>& point);
    //根据视图frame中已知的三维坐标的特征点，计算投影矩阵
    //match：frame和prevFrame对应的匹配
    //known：已知三维坐标的特征点标记，known[i]表示match[i]中的特征点已知三维坐标,即keypoints[frame][match[i].trainIdx]对应的三维点已知
    //points：三维点坐标
    Mat_<double> findPmatrixByKnownPoints(const vector<DMatch>&match,const bool* known,const Mat_<double>&points,int frame);
    //根据视图frame和视图prevFrame已知三维坐标的点，重建其他点
    //match：frame和prevFrame对应的匹配
    //known：已知三维坐标的特征点标记，known[i]表示match[i]中的特征点已知三维坐标,即keypoints[frame][match[i].trainIdx]对应的三维点已知
    //points：三维点坐标
    void reconstructByKnownPoints(int frame,int prevFrame,const vector<DMatch>&match,bool *known,Mat_<double> &points);
    //N视图bundle ajdustment
    //frameNum：需要做bundle ajdustment的视图数量
    //nconstframes：几个视图的投影矩阵固定不变，默认第1个视图投影矩阵不变
    //nconstpts3D：几个三维点的坐标保持不变
    //maxiter：最大迭代次数
    //verbose：调试信息
    void nviewSba(int frameNum,int nconstframes=1,int nconstpts3D=0,int maxiter=250,int verbose=0);
public:
    SFM(const Mat &cameraMatrix,const vector<vector<Point2d> > &keypoints,const vector<vector<Mat> > &fmatrices,const vector<vector<vector<DMatch>> > &matches,int nframes);
    //初始重建
    void initialReconstruct();
    //增加视图
    //frame：增加的视图编号
    void addView(int frame);
    //返回点云
    Cloud& getCloud();
};
#endif /* defined(__sfm__SFM__) */
