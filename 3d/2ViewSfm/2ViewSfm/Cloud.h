//
//  Cloud.h
//  sfm
//
//  Created by  刘骥 on 15/4/20.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#ifndef __sfm__Cloud__
#define __sfm__Cloud__

#include <iostream>
#include <fstream>
#include <vector>
#include "CloudPoint.h"
#include <opencv2/opencv.hpp>
#include "OpencvTools.h"
using namespace cv;
using namespace std;
using namespace opencvTool;
//此类描述所有的三维点
class Cloud{
private:
    //重建的所有三维点
    vector<CloudPoint> cloudPoints;
    //相机矩阵
    const Mat &cameraMatrix;
    //相机投影矩阵，pmatrices[i]表示视图i的投影矩阵
    const vector<Mat_<double> > &pmatrices;
    //视图数量
    int nframes;
public:
    Cloud(const Mat&cameraMatrix,const vector<Mat_<double> >&pmatrices,int nframes);
    //增加三维点
    void addPoint(const CloudPoint&cp);
    //计算从视图0到视图frameNum-1的投影误差
    void reprojectError(int frameNum);
    //在视图中绘制投影后的三维点
    void plotCloudPoints(const vector<Mat>&images);
    //保存三维点到图像
    void saveCloudPoints(const char* fileName);
    //保存三维点到ply文件
    void saveCloudPointsToPly(const vector<Mat>&images,const char* fileName);
    //获取当前视图中已知三维位置的点
    //frame：当前视图
    //prevFrame：前一视图
    //known：known[i]==true表示match[i]的三维位置已知，即keypoints[frame][match[i].trainIdx]对应的三维点已知
    //points：返回已知位置的三维点
    void findKnownPoints(int frame,int prevFrame,const vector<DMatch>&match,bool *known,Mat_<double>&points);
    //获取三维点数量
    int getPointSize();
    //获取三维点
    CloudPoint&getPoint(int i);
};
#endif /* defined(__sfm__Cloud__) */
