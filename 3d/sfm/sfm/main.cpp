//
//  main.cpp
//  sfm
//
//  Created by  刘骥 on 15/3/31.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//
#include <opencv2/opencv.hpp>
#include <opencv2/nonfree/nonfree.hpp>
#include <opencv2/features2d/features2d.hpp>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include "SFM.h"
using namespace std;
using namespace cv;
void loadImages(vector<Mat>&images,const char* fileName)
{
    ifstream in(fileName);
    string imageName;
    in>>imageName;
    while(in.good())
    {
        Mat im=imread(imageName);
        images.push_back(im);
        in>>imageName;
    }
}
void computeFeatures(const vector<Mat>&images,vector<vector<Point2d> >&keypoints,vector<Mat>& descriptors)
{
    SIFT sift;
    for(vector<Mat>::const_iterator it=images.begin();it!=images.end();it++)
    {
        vector<KeyPoint> kp;
        Mat des;
        sift(*it,Mat(),kp,des);
        vector<Point2d> points;
        for (vector<KeyPoint>::iterator it=kp.begin(); it!=kp.end(); it++) {
            points.push_back(it->pt);
        }
        keypoints.push_back(points);
        descriptors.push_back(des);
    }
}

//计算基础矩阵，用基础矩阵获取更好的匹配点
Mat findMatchesByFundamentalMat(const vector<Point2d>& keypoints1,const vector<Point2d>& keypoints2,const Mat& descriptor1,const Mat& descriptor2,vector<DMatch>&matches)
{
    BFMatcher bfMatcher(cv::NORM_L2,true);
    vector<DMatch> bfMatches;
    bfMatcher.match(descriptor1, descriptor2, bfMatches);
    vector<Point2d> points1;
    vector<Point2d> points2;
    for(vector<DMatch>::iterator it=bfMatches.begin();it!=bfMatches.end();it++)
    {
        points1.push_back(keypoints1[it->queryIdx]);
        points2.push_back(keypoints2[it->trainIdx]);
    }
    vector<uchar> mask;
    Mat F=cv::findFundamentalMat(points1,points2,FM_RANSAC,1,0.99,mask);
    for(vector<DMatch>::size_type i=0;i<bfMatches.size();i++)
    {
        if(mask[i])
            matches.push_back(bfMatches[i]);
    }
    return F;
}
void saveMat(const Mat&m,const char* fileName)
{
    ofstream out(fileName);
    out<<cv::format(m, "csv");
    out.close();
}
int main(int argc, const char * argv[]) {
    
    //加载相机矩阵
    FileStorage fs(argv[1],FileStorage::READ);
    Mat cameraMatrix;
    Mat distortion;
    fs["Camera_Matrix"]>>cameraMatrix;
    fs["Distortion_Coefficients"]>>distortion;
    SIFT sift(0, 3, 0.04, 10, 1.6);
    vector<Mat> images; //图片
    cout<<"加载图像"<<endl;
    loadImages(images,argv[2]);
    int nframes=(int)images.size(); //图片数量
    vector<vector<Point2d> > keypoints; //keypoints[i][j]表示第i张第j个图片的特征点
    vector<Mat> descriptors; //descriptors[i]表示第i张图片特征点的描述
    vector<vector<Mat> > fmatrices(nframes,vector<Mat>(nframes)); //fmatrices[i][j]表示第i张图片和第j张图片间的基础矩阵
    vector<vector<vector<DMatch>> > matches(nframes,vector<vector<DMatch>>(nframes)); //matches[i][j]表示第i张图片和第j张图片间特征点的匹配
    cout<<"计算图像特征"<<endl;
    computeFeatures(images,keypoints,descriptors);
    for (int i=0; i<nframes; i++) {
        for (int j=i+1; j<nframes; j++) {
            cout<<"求图像"<<i<<"和图像"<<j<<"之间的匹配以及基础矩阵"<<endl;
            vector<DMatch> match;
            Mat f=findMatchesByFundamentalMat(keypoints[i],keypoints[j],descriptors[i],descriptors[j],match);
            fmatrices[i][j]=f;
            matches[i][j]=match;
        }
    }
    SFM sfm(cameraMatrix,keypoints,fmatrices,matches,nframes);
    cout<<"初始重建(图像0和图像1)"<<endl;
    //用视图1和视图2进行重建
    sfm.initialReconstruct();
    //增加视图
    for (int frame=2; frame<nframes; frame++) {
       
        cout<<"增加视图"<<frame<<endl;
        sfm.addView(frame);
    }
    Cloud& cloud=sfm.getCloud();
    cloud.plotCloudPoints(images);
    //printCloudPoints(cloudPoints);
    cloud.saveCloudPoints("/Users/liuji/clouds.txt");
    cloud.saveCloudPointsToPly(images,"/Users/liuji/clouds.ply");
    //testTwoView(images,K,nframes,keypoints,descriptors,matches,fmatrices);
    return 0;
}
