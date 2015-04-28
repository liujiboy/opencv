//
//  Cloud.cpp
//  sfm
//
//  Created by  刘骥 on 15/4/20.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#include "Cloud.h"
Cloud::Cloud(const Mat&cameraMatrix,const vector<Mat_<double> >&pmatrices,int nframes):cameraMatrix(cameraMatrix),pmatrices(pmatrices),nframes(nframes)
{
    
}
void Cloud::addPoint(const CloudPoint &cp)
{
    this->cloudPoints.push_back(cp);
}
void Cloud::reprojectError(int frameNum)
{
    
    for (int frame=0; frame<frameNum; frame++) {
        Mat_<double> p=pmatrices[frame]; //投影矩阵
        double error=0; //误差
        double n=0; //点数
        for (int i=0; i<cloudPoints.size(); i++) {
            CloudPoint cp=cloudPoints[i];
            Point2d point2d;
            if(cp.getPointInFrame(frame, point2d))
            {
                n++;
                Point2d projected=cp.project(cameraMatrix, p);
                //计算误差
                double ex=projected.x-point2d.x;
                double ey=projected.y-point2d.y;
                error+=sqrt(ex*ex+ey*ey);
            }
        }
        cout<<"视图"<<frame<<"投影点数为："<<n<<" 平均投影误差为："<<error/n<<endl;
    }

}
void Cloud::plotCloudPoints(const vector<Mat>&images)
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
                Point2d projected=cp.project(cameraMatrix, p);
                circle(image, projected, 3, Scalar(255,0,0),-1);
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
void Cloud::saveCloudPoints(const char *fileName)
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
void Cloud::saveCloudPointsToPly(const vector<Mat>&images,const char* fileName)
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
void Cloud::findKnownPoints(int frame,int prevFrame,const vector<DMatch>&match,bool *known,Mat_<double>&points)
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
int Cloud::getPointSize()
{
    return (int)cloudPoints.size();
}
CloudPoint&Cloud::getPoint(int i)
{
    return cloudPoints[i];
}
