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
#include <sstream>
#include "opencv_sba.h"
#include "CloudPoint.h"
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
//返回图像中的匹配点
//points1和points2是3*N的矩阵，每列是一个齐次坐标（x,y,1）
void getMatchPoints(const vector<Point2d>& keypoints1,const vector<Point2d>& keypoints2,const vector<DMatch>&matches,Mat_<double>&points1,Mat_<double>&points2)
{
    int col=(int)matches.size();
    points1.create(3, col);
    points2.create(3, col);
    for(vector<DMatch>::size_type i=0;i<matches.size();i++)
    {
        DMatch match=matches[i];
        Point2d point1=keypoints1[match.queryIdx];
        Point2d point2=keypoints2[match.trainIdx];
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
}
/*
void triangulatePoint(const Point2d&point1,Point2d&point2,const Mat_<double>&P1,const Mat_<double>&P2,Point3d& point)
{
    Mat_<double> M=Mat::zeros(6, 6, CV_64F);
    P1.copyTo(M(Rect(0,0,4,3)));//Rect(x,y,width,height),矩阵的row对应y和height，col对应x和width
    P2.copyTo(M(Rect(0,3,4,3)));
    
    M(0,4)=-point1.x;
    M(1,4)=-point1.y;
    M(2,4)=-1;
    M(3,5)=-point2.x;
    M(4,5)=-point2.y;
    M(5,5)=-1;
    Mat_<double> X;
    SVD::solveZ(M, X);
    X=X/X(3);
    point.x=X(0);
    point.y=X(1);
    point.z=X(2);
}
*/
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
Mat_<double> trianglulateAndFindCameraMatrix(const Mat_<double>&points1,const Mat_<double>&points2,const Mat&cameraMatrix,const Mat&P1,const Mat&F,Mat&P2)
{
    
    Mat_<double> E=cameraMatrix.t()*F*cameraMatrix;
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
        triangulate(cameraMatrix.inv()*points1, cameraMatrix.inv()*points2, P1,p2, points3dVector[i]);
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
/*
Scalar reprojectError(const Mat_<double>& points2d,const Mat& cameraMatrix,const Mat&P,const Mat_<double> points3d)
{
    vector<double> errors;
    //cout<<P<<endl;
    //将3d点投影到2d
    Mat_<double> points3d_projected=cameraMatrix*P*points3d;
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
        //cout<<e<<endl;
        errors.push_back(e);
    }
    return mean(errors);
}*/
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
/*
void testTwoView(const vector<Mat>&images,const Mat&cameraMatrix,int nframes,const vector<vector<Point2d> > &keypoints,const vector<Mat> &descriptors,vector<vector<vector<DMatch>> >&matches,const vector<vector<Mat> > &fmatrices)
{
    for (int i=0; i<nframes; i++) {
        for (int j=i+1; j<nframes; j++) {
            Mat_<double>points0,points1;
            getMatchPoints(keypoints[i],keypoints[j],matches[i][j] , points0, points1);
            cout<<"对图像"<<i<<"和图像"<<j<<"的匹配点进行三角化"<<endl;
            
            // Mat K=Mat(Matx33d( 2394,0,932,0,2398,628,0,0,1));
            
            Mat P1=Mat(Matx34d(1,0,0,0,0,1,0,0,0,0,1,0));
            Mat_<double> P2;
            Mat_<double> points=trianglulateAndFindCameraMatrix(points0, points1, cameraMatrix, P1, fmatrices[i][j], P2);
            cout<<"计算图像"<<i<<"和图像"<<j<<"三角化的误差"<<endl;
            cout<<reprojectError(points0,cameraMatrix,P1,points)<<endl;
            cout<<reprojectError(points1,cameraMatrix,P2,points)<<endl;
            
 
             saveMat(points0, "/Users/liuji/Documents/workspace/SFM/x1.txt");
             saveMat(points1, "/Users/liuji/Documents/workspace/SFM/x2.txt");
             saveMat(P1, "/Users/liuji/Documents/workspace/SFM/P1.txt");
             saveMat(P2, "/Users/liuji/Documents/workspace/SFM/P2.txt");
             saveMat(points, "/Users/liuji/Documents/workspace/SFM/X.txt");
            Mat_<double> new_P1;
            Mat_<double> new_P2;
            Mat_<double> new_points;
            opencv_twoview_sba(cameraMatrix,P1,P2,points,points0,points1,new_P1,new_P2,new_points);
            cout<<reprojectError(points0,cameraMatrix,new_P1,new_points)<<endl;
            cout<<reprojectError(points1,cameraMatrix,new_P2,new_points)<<endl;
            Mat outImg1=images[i].clone();
            Mat outImg2=images[j].clone();
            plotCircle(outImg1, points0, Scalar(255,0,0));
            plotCircle(outImg1, cameraMatrix*new_P1*new_points, Scalar(0,0,255));
            plotCircle(outImg2, points1, Scalar(255,0,0));
            plotCircle(outImg2, cameraMatrix*new_P2*new_points, Scalar(0,0,255));
            showImage("outImg1", outImg1);
            showImage("outImg2", outImg2);
            waitKey();
        }
    }
}*/
void reprojectError(const Mat_<double>&cameraMatrix,const vector<CloudPoint>&cloudPoints,const vector<Mat_<double> >&pmatrices,int nframes)
{
    for (int frame=0; frame<nframes; frame++) {
        Mat_<double> p=pmatrices[frame];
        double error=0;
        double n=0;
        for (int i=0; i<cloudPoints.size(); i++) {
            CloudPoint cp=cloudPoints[i];
            Point2d point2d;
            if(cp.getPointInFrame(frame, point2d))
            {
                n++;
                Mat_<double> point(4,1);
                point(0)=cp.x;
                point(1)=cp.y;
                point(2)=cp.z;
                point(3)=1;
                Mat_<double> projected=cameraMatrix*p*point;
                projected=projected/projected(2);
                double ex=projected(0)-point2d.x;
                double ey=projected(1)-point2d.y;
                error+=sqrt(ex*ex+ey*ey);
            }
        }
        cout<<"视图"<<frame<<"投影点数为："<<n<<" 平均投影误差为："<<error/n<<endl;
    }
}
void initialReconstruct(const vector<Mat>&images,const Mat&cameraMatrix,int nframes,const vector<vector<Point2d> > &keypoints,const vector<Mat> &descriptors,vector<vector<vector<DMatch>> >&matches,const vector<vector<Mat> > &fmatrices, vector<CloudPoint>&cloudPoints,vector<Mat_<double> >&pmatrices)
{
    Mat_<double>points0,points1;
    getMatchPoints(keypoints[0],keypoints[1],matches[0][1] , points0, points1);
    cout<<"对图像0和图像1的匹配点进行三角化"<<endl;
    Mat P1=Mat(Matx34d(1,0,0,0,0,1,0,0,0,0,1,0));
    Mat_<double> P2;
    Mat_<double> points=trianglulateAndFindCameraMatrix(points0, points1, cameraMatrix, P1, fmatrices[0][1], P2);
    for (int i=0; i<points.cols; i++) {
        double x=points(0,i);
        double y=points(1,i);
        double z=points(2,i);
        CloudPoint cp(x,y,z,nframes,keypoints);
        cp.setPointIndex(0, matches[0][1][i].queryIdx);
        cp.setPointIndex(1, matches[0][1][i].trainIdx);
        cloudPoints.push_back(cp);
    }
    pmatrices[0]=P1;
    pmatrices[1]=P2;
    //reprojectError(cameraMatrix, cloudPoints, pmatrices, 2);
    cout<<"三角化完成，开始bundle adjustment"<<endl;
    nviewSba(cameraMatrix, pmatrices, cloudPoints, 2);
    reprojectError(cameraMatrix, cloudPoints, pmatrices, 2);
    /*
    Mat_<double> new_P1;
    Mat_<double> new_P2;
    Mat_<double> new_points;
    opencv_twoview_sba(cameraMatrix,P1,P2,points,points0,points1,new_P1,new_P2,new_points);
    cout<<"重建误差为:"<<reprojectError(points0,cameraMatrix,new_P1,new_points)<<" "<<reprojectError(points1,cameraMatrix,new_P2,new_points)<<endl;
    pmatrices[0]=new_P1;
    pmatrices[1]=new_P2;
    for (int i=0; i<new_points.cols; i++) {
        double x=new_points(0,i);
        double y=new_points(1,i);
        double z=new_points(2,i);
        CloudPoint cp(x,y,z,nframes,keypoints);
        cp.setPointIndex(0, matches[0][1][i].queryIdx);
        cp.setPointIndex(1, matches[0][1][i].trainIdx);
        cloudPoints.push_back(cp);
    }*/
}
void printCloudPoints(const vector<CloudPoint>&cloudPoints)
{
    for (vector<CloudPoint>::const_iterator it=cloudPoints.begin(); it!=cloudPoints.end(); it++) {
        cout<<*it<<endl;
    }
}
void plotCloudPoints(const Mat_<double>&cameraMatrix,const vector<Mat>&images,const vector<CloudPoint>&cloudPoints,const vector<Mat_<double> >&pmatrices,int nframes)
{
    for (int frame=0;frame<nframes; frame++) {
        Mat image=images[frame].clone();
        Mat_<double> p=pmatrices[frame];
        for (int i=0; i<cloudPoints.size(); i++) {
            Point2d point;
            CloudPoint cp=cloudPoints[i];
            if(cp.getPointInFrame(frame, point))
            {
                circle(image, point, 3, Scalar(0,0,255),-1);
                Mat_<double> point3d(4,1);
               
                point3d(0)=cp.x;
                point3d(1)=cp.y;
                point3d(2)=cp.z;
                point3d(3)=1;
                Mat_<double> projected=cameraMatrix*p*point3d;
                projected=projected/projected(2);
                circle(image, Point2d(projected(0),projected(1)), 3, Scalar(255,0,0),-1);
            }
        }
        string name;
        stringstream ss;
        ss<<"view:"<<frame;
        ss>>name;
        showImage(name, image);
        waitKey();
        
    }
}
void findKnownPoints(vector<CloudPoint>&cloudPoints,int frame,int prevFrame,const vector<DMatch>&match,bool *known,Mat_<double>&points)
{
    //在前一个视图中查找已知的三维点
    //设当前视图为frame，则前一视图为prevFrame=frame-1
    //前一视图的点是a(索引)，当前视图的匹配点为b
    //若在cloudsPoint中存在点cp，并且cp.pointIndicesInFrame[prevFrame]==a
    //则当前视图中的点a与点cp对应
    for(int i=0;i<match.size();i++){
        int a=match[i].queryIdx;
        int b=match[i].trainIdx;
        for(vector<CloudPoint>::iterator it=cloudPoints.begin();it!=cloudPoints.end();it++)
        {
            if (it->getPointIndex(prevFrame)==a) {
                it->setPointIndex(frame, b);
                known[i]=true;
                points(0,i)=it->x;
                points(1,i)=it->y;
                points(2,i)=it->z;
                points(3,i)=1;
                break;
            }
        }
    }
}
Mat_<double> findPmatrixByKnownPoints(const vector<DMatch>&match,bool* known,const vector<vector<Point2d> > &keypoints,const Mat_<double>&cameraMatrix,const Mat_<double>&points,int frame)
{
    //根据已知的三位点确定投影矩阵
    vector<Point3f> objectPoints;
    vector<Point2f> imagePoints;
    Mat_<double> rvec;
    Mat_<double> tvec;
    for(int i=0;i<match.size();i++)
    {
        if(known[i])
        {
            Point3f p3d(points(0,i),points(1,i),points(2,i));
            objectPoints.push_back(p3d);
            int index=match[i].trainIdx;
            Point2f p2d=keypoints[frame][index];
            imagePoints.push_back(p2d);
        }
    }
    solvePnPRansac(objectPoints, imagePoints, cameraMatrix, noArray(), rvec, tvec);
    Mat_<double> rmat;
    Rodrigues(rvec, rmat);
    Mat_<double> p;
    p.create(3, 4);
    rmat.copyTo(p(Rect(0,0,3,3)));
    p(0,3)=tvec(0);
    p(1,3)=tvec(1);
    p(2,3)=tvec(2);
    return p;
}
void reconstructByKnownPoints(int nframes,int frame,int prevFrame,const vector<vector<Point2d> > &keypoints,const vector<DMatch>&match,const Mat_<double>&cameraMatrix,bool *known,Mat_<double> &points,vector<CloudPoint>&cloudPoints,vector<Mat_<double> >&pmatrices)
{
    //重建未知的三维点
    //Mat_<double>points0,points1;
    //getMatchPoints(keypoints[prevFrame],keypoints[frame],match , points0, points1);
    //Mat_<double> kpoints0=cameraMatrix.inv()*points0;
    //Mat_<double> kpoints1=cameraMatrix.inv()*points1;
    Mat_<double> invCameraMatrix=cameraMatrix.inv();
    for (int i=0; i<match.size(); i++) {
        if(!known[i])
        {
            int a=match[i].queryIdx;
            int b=match[i].trainIdx;
            Mat_<double> point1(3,1);
            Mat_<double> point2(3,1);
            point1(0)=keypoints[prevFrame][a].x;
            point1(1)=keypoints[prevFrame][a].y;
            point1(2)=1;
            point2(0)=keypoints[frame][a].x;
            point2(1)=keypoints[frame][b].y;
            point2(2)=1;
            point1=invCameraMatrix*point1;
            point2=invCameraMatrix*point2;
            Mat_<double> point(4,1);
            triangulatePoint(point1,point2,pmatrices[prevFrame],pmatrices[frame],point);
            CloudPoint cp(point(0),point(1),point(2),nframes,keypoints);
            cp.setPointIndex(prevFrame, a);
            cp.setPointIndex(frame, b);
            cloudPoints.push_back(cp);
        }
    }

    /*
    cout<<"重建误差为:"<<reprojectError(points0,cameraMatrix,pmatrices[prevFrame],points)<<" "<<reprojectError(points1,cameraMatrix,pmatrices[frame],points)<<endl;
    cout<<"开始bundle adjustment"<<endl;
    //把known[i]==true的点放在前面不进行优化
    int nconstpts3D=0;//表示从第几个点开始优化
    Mat_<double> newPoints0(3,(int)match.size());
    Mat_<double> newPoints1(3,(int)match.size());
    Mat_<double> newPoints(4,(int)match.size());
    
    for(int i=0;i<match.size();i++)
    {
        if (known[i]) {
            points0.col(i).copyTo(newPoints0.col(nconstpts3D));
            points1.col(i).copyTo(newPoints1.col(nconstpts3D));
            points.col(i).copyTo(newPoints.col(nconstpts3D));
            nconstpts3D++;
        }
    }
    //把known[i]==false的点放在后面进行优化
    int j=nconstpts3D;
    for(int i=0;i<match.size();i++)
    {
        if (!known[i]) {
            points0.col(i).copyTo(newPoints0.col(j));
            points1.col(i).copyTo(newPoints1.col(j));
            points.col(i).copyTo(newPoints.col(j));
            j++;
        }
    }
    Mat_<double> refinedP1;
    Mat_<double> refinedP2;
    Mat_<double> refinedPoints;
    opencv_twoview_sba(cameraMatrix,pmatrices[prevFrame],pmatrices[frame],newPoints,newPoints0,newPoints1,refinedP1,refinedP2,refinedPoints,nconstpts3D);
    //将优化后的点拷贝回去
    j=nconstpts3D;
    for(int i=0;i<match.size();i++)
    {
        if (!known[i]) {
            refinedPoints.col(j).copyTo(points.col(i));
            j++;
        }
    }
    pmatrices[frame]=refinedP2;
    cout<<"重建误差为:"<<reprojectError(points0,cameraMatrix,pmatrices[prevFrame],points)<<" "<<reprojectError(points1,cameraMatrix,pmatrices[frame],points)<<endl;
    for(int i=0;i<match.size();i++)
    {
        if (!known[i]) {
            CloudPoint cp(points(0,i),points(1,i),points(2,i),nframes,keypoints);
            cp.setPointIndex(prevFrame, match[i].queryIdx);
            cp.setPointIndex(frame, match[i].trainIdx);
            cloudPoints.push_back(cp);
        }
    }
    
    Mat outImg1=images[prevFrame].clone();
    Mat outImg2=images[frame].clone();
    plotCircle(outImg1, points0, Scalar(255,0,0));
    plotCircle(outImg1, cameraMatrix*pmatrices[prevFrame]*points, Scalar(0,0,255));
    plotCircle(outImg2, points1, Scalar(255,0,0));
    plotCircle(outImg2, cameraMatrix*pmatrices[frame]*points, Scalar(0,0,255));
    showImage("outImg1", outImg1);
    showImage("outImg2", outImg2);
    waitKey();*/

}
void addView(int nframes,int frame,const Mat_<double>&cameraMatrix,const vector<vector<Point2d> > &keypoints,vector<vector<vector<DMatch>> >& matches,vector<CloudPoint>&cloudPoints,vector<Mat_<double> >&pmatrices)
{
    int prevFrame=frame-1;
    vector<DMatch> match=matches[prevFrame][frame];
    bool *known=new bool[(int)match.size()]();
    Mat_<double> points(4,(int)match.size()); //根据当前视图和前一视图重建的三维点
    findKnownPoints(cloudPoints,frame,prevFrame,match,known,points);
    pmatrices[frame]=findPmatrixByKnownPoints(match,known,keypoints,cameraMatrix,points,frame);
    reconstructByKnownPoints(nframes,frame,prevFrame,keypoints,match,cameraMatrix,known,points,cloudPoints,pmatrices);
    //reprojectError(cameraMatrix, cloudPoints, pmatrices, frame+1);
    nviewSba(cameraMatrix, pmatrices, cloudPoints, frame+1);
    reprojectError(cameraMatrix, cloudPoints, pmatrices, frame+1);
           /*
     Mat_<double> new_P1;
     Mat_<double> new_P2;
     Mat_<double> new_points;
     opencv_twoview_sba(cameraMatrix,pmatrices[prevFrame],pmatrices[frame],points,points0,points1,new_P1,new_P2,new_points,250,10);
     cout<<"重建误差为:"<<reprojectError(points0,cameraMatrix,new_P1,new_points)<<" "<<reprojectError(points1,cameraMatrix,new_P2,new_points)<<endl;
     */
    delete known;

}
void saveCloudPoints(const vector<CloudPoint>&cloudPoints,const char* fileName)
{
    ofstream out(fileName);
    Mat_<double> m((int)cloudPoints.size(),3);
    for(int i=0;i<cloudPoints.size();i++)
    {
        m(i,0)=cloudPoints[i].x;
        m(i,1)=cloudPoints[i].y;
        m(i,2)=cloudPoints[i].z;
    }
    out<<cv::format(m, "csv");
    out.close();
}
void saveCloudPointsToPly(int nframes,const vector<Mat>&images,const vector<CloudPoint>&cloudPoints,const char* fileName)
{
    ofstream out(fileName);
    out<<"ply"<<endl;
    out<<"format ascii 1.0"<<endl;
    out<<"element vertex "<<cloudPoints.size()<<endl;
    out<<"property float x"<<endl;
    out<<"property float y"<<endl;
    out<<"property float z"<<endl;
    out<<"property uchar diffuse_red"<<endl;
    out<<"property uchar diffuse_green"<<endl;
    out<<"property uchar diffuse_blue"<<endl;
    out<<"end_header"<<endl;
    for(int i=0;i<cloudPoints.size();i++)
    {
        const CloudPoint& cp=cloudPoints[i];
        out<<cp.x<<" "<<cp.y<<" "<<cp.z<<" ";
        for(int frame=0;frame<nframes;frame++)
        {
            Point2d p;
            if(cp.getPointInFrame(frame, p))
            {
                Vec3b color=images[frame].at<Vec3b>(p.y,p.x);
                out<<(int)color[2]<<" "<<(int)color[1]<<" "<<(int)color[0]<<endl;
                break;
            }
        }
    }
    out.close();
}
int main(int argc, const char * argv[]) {
    vector<CloudPoint> cloudPoints;
    //加载相机矩阵
    FileStorage fs(argv[1],FileStorage::READ);
    //Mat K=Mat(Matx33d(8.6454312993195060e+02,0.0,4.9950000000000000e+02,0.0,8.6454312993195060e+02,3.7450000000000000e+02,0.0,0.0,1.0));
    Mat cameraMatrix;
    Mat distortion;
    fs["Camera_Matrix"]>>cameraMatrix;
    fs["Distortion_Coefficients"]>>distortion;
    SIFT sift(0, 3, 0.04, 10, 1.6);
    int nframes; //图片数量
    vector<Mat> images; //图片
    
    cout<<"加载图像"<<endl;
    loadImages(images,argv[2]);
    nframes=(int)images.size();
    vector<vector<Point2d> > keypoints; //keypoints[i][j]表示第i张第j个图片的特征点
    vector<Mat> descriptors; //descriptors[i]表示第i张图片特征点的描述
    vector<vector<Mat> > fmatrices(nframes,vector<Mat>(nframes)); //fmatrices[i][j]表示第i张图片和第j张图片间的基础矩阵
    vector<vector<vector<DMatch>> > matches(nframes,vector<vector<DMatch>>(nframes)); //matches[i][j]表示第i张图片和第j张图片间特征点的匹配
    vector<Mat_<double> > pmatrices(nframes); //pmatrices[i]表示第i张图片的投影矩阵
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
    cout<<"初始重建(图像0和图像1)"<<endl;
    //用视图1和视图2进行重建
    initialReconstruct(images,cameraMatrix,nframes,keypoints,descriptors,matches,fmatrices, cloudPoints, pmatrices);
    
    //增加视图
    for (int frame=2; frame<nframes; frame++) {
       
        cout<<"增加视图"<<frame<<endl;
        addView(nframes,frame,cameraMatrix,keypoints,matches,cloudPoints,pmatrices);
    }
    plotCloudPoints(cameraMatrix, images, cloudPoints, pmatrices, nframes);
    //printCloudPoints(cloudPoints);
    saveCloudPoints(cloudPoints, "/Users/liuji/clouds.txt");
    saveCloudPointsToPly(nframes, images, cloudPoints, "/Users/liuji/clouds.ply");
    //testTwoView(images,K,nframes,keypoints,descriptors,matches,fmatrices);
    return 0;
}
