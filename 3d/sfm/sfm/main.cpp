//
//  main.cpp
//  sfm
//
//  Created by  刘骥 on 15/3/31.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//
#include "opencv2/opencv.hpp"
#include "opencv2/nonfree/nonfree.hpp"
#include "opencv2/features2d/features2d.hpp"
#include <iostream>
#include <vector>
#include <fstream>
#include "opencv_sba.h"
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
void computeFeatures(const vector<Mat>&images,vector<vector<KeyPoint> >&keypointsVector,vector<Mat>& descriptorVector)
{
    SIFT sift;
    for(vector<Mat>::const_iterator it=images.begin();it!=images.end();it++)
    {
        vector<KeyPoint> keypoints;
        Mat descriptor;
        sift(*it,Mat(),keypoints,descriptor);
        keypointsVector.push_back(keypoints);
        descriptorVector.push_back(descriptor);
    }
}
void computeP(const Mat&E,vector<Mat>&pVector)
{
    SVD svd(E);
    Mat u=svd.u;
    Mat vt=svd.vt;
    if(determinant(u*vt)<0)
    {
        vt=-vt;
    }
    //Matx33d s(1,0,0,0,1,0,0,0,0);
    //Mat e=u*Mat(s)*vt;
    Matx33d w(0,-1,0,1,0,0,0,0,1);
    Mat_<double> R1=u*Mat(w)*vt;
    Mat_<double> R2=u*Mat(w).t()*vt;
    Mat_<double> t1=u.col(2);
    Mat_<double> t2=-t1;
    Matx34d P1(R1(0,0),R1(0,1),R1(0,2),t1(0),
              R1(1,0),R1(1,1),R1(1,2),t1(1),
              R1(2,0),R1(2,1),R1(2,2),t1(2));
    Matx34d P2(R2(0,0),R2(0,1),R2(0,2),t1(0),
               R2(1,0),R2(1,1),R2(1,2),t1(1),
               R2(2,0),R2(2,1),R2(2,2),t1(2));
    Matx34d P3(R1(0,0),R1(0,1),R1(0,2),t2(0),
               R1(1,0),R1(1,1),R1(1,2),t2(1),
               R1(2,0),R1(2,1),R1(2,2),t2(2));
    Matx34d P4(R2(0,0),R2(0,1),R2(0,2),t2(0),
               R2(1,0),R2(1,1),R2(1,2),t2(1),
               R2(2,0),R2(2,1),R2(2,2),t2(2));
    pVector.push_back(Mat(P1));
    pVector.push_back(Mat(P2));
    pVector.push_back(Mat(P3));
    pVector.push_back(Mat(P4));
    
   /* cout<<P1<<endl;
    cout<<P2<<endl;
    cout<<P3<<endl;
    cout<<P4<<endl;*/
    
}

//计算基础矩阵，用基础矩阵获取更好的匹配点
Mat findMatchesByFundamentalMat(const vector<KeyPoint>& keypoints1,const vector<KeyPoint>& keypoints2,const Mat& descriptor1,const Mat& descriptor2,vector<DMatch>&matches)
{
    BFMatcher bfMatcher(cv::NORM_L2,true);
    vector<DMatch> bfMatches;
    bfMatcher.match(descriptor1, descriptor2, bfMatches);
    vector<Point2d> points1;
    vector<Point2d> points2;
    for(vector<DMatch>::iterator it=bfMatches.begin();it!=bfMatches.end();it++)
    {
        points1.push_back(keypoints1[it->queryIdx].pt);
        points2.push_back(keypoints2[it->trainIdx].pt);
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
//返回图像中的匹配点
//points1和points2是3*N的矩阵，每列是一个齐次坐标（x,y,1）
void getMatchPoints(const vector<KeyPoint>& keypoints1,const vector<KeyPoint>& keypoints2,const vector<DMatch>&matches,Mat_<double>&points1,Mat_<double>&points2)
{
    int col=(int)matches.size();
    points1.create(3, col);
    points2.create(3, col);
    for(vector<DMatch>::size_type i=0;i<matches.size();i++)
    {
        DMatch match=matches[i];
        Point2d point1=keypoints1[match.queryIdx].pt;
        Point2d point2=keypoints2[match.trainIdx].pt;
        points1(0,(int)i)=point1.x;
        points1(1,(int)i)=point1.y;
        points1(2,(int)i)=1;
        points2(0,(int)i)=point2.x;
        points2(1,(int)i)=point2.y;
        points2(2,(int)i)=1;
    }
}
//三角化点
//point1和point2是3x1的矩阵
//point是4x1的矩阵
void triangulatePoint(const Mat_<double>&point1,const Mat_<double>&point2,const Mat_<double>&P1,const Mat_<double>&P2,Mat_<double>& point)
{
    Mat_<double> M=Mat::zeros(6, 6, CV_64F);
    P1.copyTo(M(Rect(0,0,4,3)));//Rect(x,y,width,height),矩阵的row对应y和height，col对应x和width
    P2.copyTo(M(Rect(0,3,4,3)));

    M(0,4)=-point1(0);
    M(1,4)=-point1(1);
    M(2,4)=-point1(2);
    M(3,5)=-point2(0);
    M(4,5)=-point2(1);
    M(5,5)=-point2(2);
    Mat_<double> X;
    SVD::solveZ(M, X);
    X=X/X(3);
    point(0,0)=X(0);
    point(1,0)=X(1);
    point(2,0)=X(2);
    point(3,0)=1;
   /* Mat_<double> A=Mat_<double>::zeros(6, 4);
    Mat_<double> B=Mat_<double>::zeros(6, 1);
    Mat_<double> X;
    P1.copyTo(A(Rect(0,0,4,3)));
    P2.copyTo(A(Rect(0,3,4,3)));
    B(0)=point1(0);
    B(1)=point1(1);
    B(2)=point1(2);
    B(3)=point2(0);
    B(4)=point2(1);
    B(5)=point2(2);
    cv::solve(A, B, X,DECOMP_SVD);
    cout<<X<<endl;
    cout<<norm(A*X-B)<<endl;
    cin.get();
    X=X/X(3);
    point(0,0)=X(0);
    point(1,0)=X(1);
    point(2,0)=X(2);
    point(3,0)=1;
    cout<<A<<endl;
    cout<<B<<endl;
    cout<<X<<endl;
    cout<<A*X-B<<endl;
    cin.get();*/
    /*Matx<double,6,4> A;
    A(0) = point1(0)*P1(2)-P1(0);
    A(1) = point1(1)*P1(2)-P1(1);
    A(2) = point1(0)*P1(1)-point1(1)*P1(0);
    A(3) = point2(0)*P2(2)-P2(0);
    A(4) = point2(1)*P2(2)-P2(1);
    A(5) = point2(0)*P2(1)-point2(1)*P2(0);
    Mat_<double> X;
    SVD::solveZ(A, X);
    X=X/X(3);
    point(0,0)=X(0);
    point(1,0)=X(1);
    point(2,0)=X(2);
    point(3,0)=1;*/
}
void saveMat(const Mat&m,const char* fileName)
{
    ofstream out(fileName);
    out<<cv::format(m, "csv");
    out.close();
}
void triangulate(const Mat_<double>&points1,const Mat_<double>&points2,const Mat&P1,const Mat&P2,Mat_<double>&points)
{
    /*saveMat(points1,"/Users/liuji/x1.txt");
    saveMat(points2,"/Users/liuji/x2.txt");
    saveMat(P1,"/Users/liuji/P1.txt");
    saveMat(P2,"/Users/liuji/P2.txt");*/
    points.create(4, points1.cols);
    for(int i=0;i<points1.cols;i++)
    {
      
        Mat_<double> point=points.col(i);
    
        triangulatePoint(points1.col(i), points2.col(i), P1, P2,point);
       // points(0,i)=point(0);
        //points(1,i)=point(0);
        //points(2,i)=point(0);
        //points(3,i)=point(0);
    }
    /*saveMat(points,"/Users/liuji/points.txt");
    cout<<"ok"<<endl;
    cin.get();*/
 }
Mat_<double> trianglulateAndFindCameraMatrix(const Mat_<double>&points1,const Mat_<double>&points2,const Mat&K,const Mat&P1,const Mat&F,Mat&P2)
{

    Mat_<double> E=K.t()*F*K;
   // E=E/E(2,2);
    //cout<<E<<endl;
   // cout<<E<<endl;
    //Mat x1n=K.inv()*points1;
    //Mat x2n=K.inv()*points2;
    //Mat E=findFundamentalMat(x1n.rowRange(0, 2).reshape(2,1),x2n.rowRange(0, 2).reshape(2,1));

    //cout<<cv::format(E,"numpy")<<endl;
    vector<Mat> P2Vector;
    computeP(E,P2Vector);
    // vector<vector<Point3d> > pointsVector;
    int maxCount=0;
    vector<Mat>::size_type bestIndex=0;
    vector<Mat_<double> > points3dVector(4);
    for(vector<Mat>::size_type i=0;i<P2Vector.size();i++)
    {
        Mat p2=P2Vector[i];
        //cout<<"P2------------"<<endl;
        //cout<<p2<<endl;
       // triangulatePoints(P1, p2, points1.rowRange(0, 2), points2.rowRange(0,2), points3dVector[i]);
        //triangulate(points1, points2, K*P1,K*p2, points3dVector[i]);
        triangulate(K.inv()*points1, K.inv()*points2, P1,p2, points3dVector[i]);
        //计算在相机前面的点数量
        //将计算出的3d点用p2进行投影（此处只考虑了p2，因为P1是固定的）
        //统计投影后z坐标大于0的点的数目
        //取得最大数目的p2就是最佳的
        //vector<Point3d> points3d_projected(points3d.size());
        Mat P4x4=Mat::eye(4, 4, CV_64F);
        Mat P3x4=P4x4(Rect(0,0,4,3));
        p2.copyTo(P3x4);
        Mat_<double> points3d_projected=P3x4*points3dVector[i];
        int count=0;
        for(int colIdx=0;colIdx<points3d_projected.cols;colIdx++)
        {
            Mat_<double> point3d_projected=points3d_projected.col(colIdx);
            if(points3d_projected(2)>0)
            {
                count++;
            }
        }
        if(count>maxCount)
        {
            maxCount=count;
            bestIndex=i;
        }
    }
    P2Vector[bestIndex].copyTo(P2);
    return points3dVector[bestIndex];
}
Scalar reprojectError(const Mat_<double>& points2d,const Mat& K,const Mat&P,const Mat_<double> points3d)
{
    vector<double> errors;
    //cout<<P<<endl;
    //将3d点投影到2d
    Mat_<double> points3d_projected=K*P*points3d;
    //将投影后的2d点，归一化为(x/z,y/z,1)
    //将归一化后的点与points2d相减
    for(int i=0;i<points3d_projected.cols;i++)
    {
        Mat_<double> point3d=points3d_projected.col(i);
        point3d=point3d/point3d(2);
        point3d(2)=1;//避免0.99999999的情况
       // cout<<points2d.col(i)<<endl;
       // cout<<point3d<<endl;
        double e=norm(points2d.col(i)-point3d);
       // cout<<e<<endl;
        errors.push_back(e);
    }
    return mean(errors);
}
void plotCircle( Mat&img,const Mat_<double>&points,const Scalar &color,int radius=3,int thickness=-1, int lineType=8, int shift=0)
{

    for(int i=0;i<points.cols;i++)
    {
        Point2d center;
        center.x=points(0,i)/points(2,i);
        center.y=points(1,i)/points(2,i);
        circle(img, center, radius, color,thickness);
    }
}
void showImage(const string&name,const Mat&img,int flags = WINDOW_NORMAL)
{
    namedWindow(name,flags);
    imshow(name, img);
}
int main(int argc, const char * argv[]) {
    Mat a(3,3,CV_64F);
    Mat b(3,3,CV_64F);
    cout<<a*b<<endl;
    SIFT sift;
    vector<Mat> images;
    vector<vector<KeyPoint> > keypointsVector;
    vector<Mat> descriptorVector;
    cout<<"加载图像----------------"<<endl;
    loadImages(images,argv[1]);
    cout<<"计算图像特征----------------"<<endl;
    computeFeatures(images,keypointsVector,descriptorVector);
    cout<<"求图像0和图像1之间的匹配以及基础矩阵----------------"<<endl;
    vector<DMatch> matches01;
    Mat F01;
    F01=findMatchesByFundamentalMat(keypointsVector[0],keypointsVector[1],descriptorVector[0],descriptorVector[1],matches01);
    //cout<<F01<<endl;
    Mat_<double>points0,points1;
    getMatchPoints(keypointsVector[0], keypointsVector[1], matches01, points0, points1);
        cout<<"对图像0和图像1的匹配点进行三角化----------------"<<endl;
   Mat K=Mat(Matx33d(8.6454312993195060e+02,0.0,4.9950000000000000e+02,0.0,8.6454312993195060e+02,3.7450000000000000e+02,0.0,0.0,1.0));
   // Mat K=Mat(Matx33d( 2394,0,932,0,2398,628,0,0,1));
  
    Mat P1=Mat(Matx34d(1,0,0,0,0,1,0,0,0,0,1,0));
    Mat_<double> P2;
    Mat_<double> points=trianglulateAndFindCameraMatrix(points0, points1, K, P1, F01, P2);
    cout<<"计算图像0和图像1三角化的误差----------------"<<endl;
    cout<<reprojectError(points0,K,P1,points)<<endl;
    cout<<reprojectError(points1,K,P2,points)<<endl;
  

    /*
    saveMat(points0, "/Users/liuji/Documents/workspace/SFM/x1.txt");
    saveMat(points1, "/Users/liuji/Documents/workspace/SFM/x2.txt");
    saveMat(P1, "/Users/liuji/Documents/workspace/SFM/P1.txt");
    saveMat(P2, "/Users/liuji/Documents/workspace/SFM/P2.txt");
    saveMat(points, "/Users/liuji/Documents/workspace/SFM/X.txt");*/
    Mat_<double> new_P1;
    Mat_<double> new_P2;
    Mat_<double> new_points;
    opencv_twoview_sba_struct(K,P1,P2,points,points0,points1,new_P1,new_P2,new_points,250,10);
    cout<<reprojectError(points0,K,new_P1,new_points)<<endl;
    cout<<reprojectError(points1,K,new_P2,new_points)<<endl;
    Mat outImg1=images[0].clone();
    Mat outImg2=images[1].clone();
    plotCircle(outImg1, points0, Scalar(255,0,0));
    plotCircle(outImg1, K*new_P1*new_points, Scalar(0,0,255));
    plotCircle(outImg2, points1, Scalar(255,0,0));
    plotCircle(outImg2, K*new_P2*new_points, Scalar(0,0,255));
    showImage("outImg1", outImg1);
    showImage("outImg2", outImg1);
    waitKey();
    /*ofstream out1("/Users/liuji/points1.txt");
    out1<<cv::format(points1, "csv");
    out1.close();
    ofstream out2("/Users/liuji/points1_pro.txt");
    out2<<cv::format(K*P1*points, "csv");
    out2.close();*/
    /*saveMat(points0, "/Users/liuji/Documents/workspace/SFM/ch5/x1.txt");
    saveMat(K*P1*points, "/Users/liuji/Documents/workspace/SFM/ch5/x1p.txt");
    saveMat(points1, "/Users/liuji/Documents/workspace/SFM/ch5/x2.txt");
    saveMat(K*P2*points, "/Users/liuji/Documents/workspace/SFM/ch5/x2p.txt");
    saveMat(points, "/Users/liuji/Documents/workspace/SFM/ch5/X.txt");*/
    return 0;
}
