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

//加载视图
//fileName：视图列表文件
//images：加载的视图
//images实际上是函数的返回值，为什么直接用return一个vector<Mat>类型的返回值呢？
//因为这样会导致内存的拷贝，用引用就没有这个问题。
//此处非常怀念Java，因为Java所有的对象都是引用，直接用return也不会有内存拷贝问题。
//今后在程序中，我们约定：凡是函数参数中没有加const的“引用”或者“指针”均是作为函数
//的返回值。加const的“引用”或者“指针”是为了在参数传递时不会出现内存的拷贝。
void loadImages(const char* fileName,vector<Mat>&images)
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
//计算视图特征
//sift：视图特征的计算方法
//images：视图列表
//keyPoints：返回的特征点
//descriptor：返回的特征描述
//如前所述，加const的引用表示只读，不加const的引用表示返回值
void computeFeatures(const SIFT sift,const vector<Mat>&images,vector<vector<Point2d> >&keypoints,vector<Mat>& descriptors)
{
    
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

//用两个视图的特征点与特征描述计算基础矩阵，并用ransac算法产生的mask获取更好的特征匹配
//keypoints1：视图1的特征点
//keypoints2：视图2的特征点
//descriptor1：视图1的特征描述
//descriptor2：视图2的特征描述
//betterMatch：返回更好的匹配
Mat computeFundamentalMatAndMatch(const vector<Point2d>& keypoints1,const vector<Point2d>& keypoints2,const Mat& descriptor1,const Mat& descriptor2,vector<DMatch>&betterMatch)
{
    //求得初始匹配（蛮力法）
    BFMatcher bfMatcher(cv::NORM_L2,true);
    vector<DMatch> initialMatch;
    bfMatcher.match(descriptor1, descriptor2, initialMatch);
    //计算基础矩阵
    vector<Point2d> points1;
    vector<Point2d> points2;
    for(vector<DMatch>::iterator it=initialMatch.begin();it!=initialMatch.end();it++)
    {
        points1.push_back(keypoints1[it->queryIdx]);
        points2.push_back(keypoints2[it->trainIdx]);
    }
    vector<uchar> mask;
    Mat F=cv::findFundamentalMat(points1,points2,FM_RANSAC,1,0.99,mask);
    //用ransac算法产生的mask获取更好的特征匹配
    for(vector<DMatch>::size_type i=0;i<initialMatch.size();i++)
    {
        if(mask[i])
            betterMatch.push_back(initialMatch[i]);
    }
    return F;
}
//程序的基本流程如下：
//(1)加载相机矩阵
//(2)加载视图（多个）
//(3)计算视图特征
//(4)计算基础矩阵和匹配
//(5)用视图0和视图1进行初始重建
//(6)增加视图，直到全部视图都计算完毕，重建结束
//(7)显示重建结果，保存输出
int main(int argc, const char * argv[]) {
    if(argc!=5)
    {
        cout<<"使用方法：程序名 相机参数文件路径 视图列表文件路径 输出点云文件路径 输出ply文件路径"<<endl;
        exit(1);
    }
    //(1)加载相机矩阵
    cout<<"(1)加载相机矩阵"<<endl;
    FileStorage fs(argv[1],FileStorage::READ);
    Mat cameraMatrix; //相机矩阵
    fs["Camera_Matrix"]>>cameraMatrix;
    //(2)加载视图
    cout<<"(2)加载图像"<<endl;
    vector<Mat> images; //视图
    loadImages(argv[2],images);
    int nframes=(int)images.size(); //视图数量
    //(3)计算视图特征
    cout<<"(3)计算视图特征"<<endl;
    SIFT sift(0, 3, 0.04, 10, 1.6); //使用sift计算特征
    //vector<vector<XXX> >定义一个动态二维数组结构，元素类型为XXX
    //与普通二维数组不同，此数组的行列维数可以变化
    vector<vector<Point2d> > keypoints;//keypoints[i][j]表示第i视图第j个特征点
    vector<Mat> descriptors; //descriptors[i]表示第i视图特征点的描述
    computeFeatures(sift,images,keypoints,descriptors);
    
    //(4)计算基础矩阵和匹配
    cout<<"(4)计算基础矩阵和匹配"<<endl;
    //fmatrices的定义看上去比较蛋疼
    //fmatrices是一个二维数组，它有nframes行，nframes列
    //定义成Mat fmatrices[nframes][nframes]不行吗？不行，不信你试试，编译都过不了。
    //在C和C++中部支持二维数组维数的动态定义，如要动态，只能用new来定义，结构如下：
    //Mat** fmatrices=new Mat*[nframes];
    //for(int i=0;i<nframes;i++)
    //      fmatrices[i]=new Mat[nframes];
    //怎么样更蛋疼吧？算了还是下面这个定义稍微简单一些，而且不用delete
    vector<vector<Mat> > fmatrices(nframes,vector<Mat>(nframes)); //fmatrices[i][j]表示第i视图和第j视图间的基础矩阵
    //matches也是一个二维数组，它其实应该可以优雅的定义为：
    //vector<DMatch> matches[nframes][nframes];
    //抱歉语法不支持，所以只能看到如下蛋疼的定义
    vector<vector<vector<DMatch>> > matches(nframes,vector<vector<DMatch>>(nframes)); //matches[i][j]表示第i视图和第j视图间特征点的匹配
    for (int i=0; i<nframes; i++) {
        for (int j=i+1; j<nframes; j++) {
            cout<<"计算图像"<<i<<"和图像"<<j<<"之间的匹配以及基础矩阵"<<endl;
            vector<DMatch> match;
            Mat f=computeFundamentalMatAndMatch(keypoints[i],keypoints[j],descriptors[i],descriptors[j],match);
            fmatrices[i][j]=f;
            matches[i][j]=match;
        }
    }
    
    SFM sfm(cameraMatrix,keypoints,fmatrices,matches,nframes);
    //(5)用视图0和视图1进行初始重建
    cout<<"(5)用视图0和视图1进行初始重建"<<endl;
    sfm.initialReconstruct();
    //(6)增加视图，直到全部视图都计算完毕，重建结束
    cout<<"(6)增加视图"<<endl;
    for (int frame=2; frame<nframes; frame++) {
       
        cout<<"增加视图"<<frame<<endl;
        sfm.addView(frame);
    }
    //(7)显示重建结果，保存输出
    cout<<"(7)显示重建结果，保存输出"<<endl;
    Cloud& cloud=sfm.getCloud();
    cloud.plotCloudPoints(images);
    cloud.saveCloudPoints(argv[3]);
    cloud.saveCloudPointsToPly(images,argv[4]);
    return 0;
}
