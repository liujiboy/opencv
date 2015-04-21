//
//  CloudPoint.h
//  sfm
//
//  Created by  刘骥 on 15/4/18.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#ifndef __sfm__CloudPoint__
#define __sfm__CloudPoint__

#include <iostream>
#include <vector>
#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;
//用于表示三维点云点，其结构如下：
//x y z i0  i1  i2  i3....
//x、y、z表示点的三维坐标
//i0、i1、i2、i3表示该点在视图0、视图1、视图2和视图3等中对应特征点的索引
class CloudPoint {
private:
    //各个视图的特征点索引，pointIndices[i]表示在视图i下的特征点索引
    vector<int> pointIndices;
    //特征点引用，keypoints[i][j]返回视图i中第j个特征点
    //将keypoints作为成员变量，只是为了便于取特征点
    //由于此处是引用，故此只会占用很少的内存（在64位系统下是64位）
    const vector<vector<Point2d> >& keypoints;
public:
    double x;
    double y;
    double z;
    //nframes：视图数量
    //keypoints：特征点
    CloudPoint(double x,double y,double z,int nframes,const vector<vector<Point2d> >& keypoints);
    //设置特征点索引
    //frame：视图编号
    //index：特征点索引
    void setPointIndex(int frame,int index);
    //获取特征点索引
    //frame：视图编号
    int  getPointIndex(int frame) const;
    //获取特征点，返回true表示在视图frame下特征点存在，false表示特征点不存在
    //frame：视图编号
    //point：返回的特征点
    bool getPointInFrame(int frame,Point2d&point) const;
    //将三维点投影到二维
    //cameraMatrix：相机矩阵
    //pMatrix：投影矩阵
    Point2d project(const Mat &cameraMatrix,const Mat_<double>& pMatrix);
    //重载<<运算符
    friend ostream&operator<<(ostream&out,const CloudPoint&cp);
};

#endif /* defined(__sfm__CloudPoint__) */
