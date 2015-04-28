//
//  CloudPoint.cpp
//  sfm
//
//  Created by  刘骥 on 15/4/18.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#include "CloudPoint.h"
CloudPoint::CloudPoint(double x,double y,double z,int nframes,const vector<vector<Point2d> >& keypoints):pointIndices(nframes,-1),keypoints(keypoints)
{
    this->x=x;
    this->y=y;
    this->z=z;
}
void CloudPoint::setPointIndex(int frame, int index)
{
    pointIndices[frame]=index;
}
int CloudPoint::getPointIndex(int frame)const
{
    return pointIndices[frame];
}
bool CloudPoint::getPointInFrame(int frame,Point2d&point)const
{
    
    int index=pointIndices[frame];
    if (index==-1) {
        return false;
    }else{
        point=keypoints[frame][index];
        return true;
    }
}
Point2d CloudPoint::project(const Mat &cameraMatrix,const Mat_<double>& pMatrix)
{
    //转换为齐次坐标
    Mat_<double> point3d(4,1);
    point3d(0)=x;
    point3d(1)=y;
    point3d(2)=z;
    point3d(3)=1;
    //投影到二维
    Mat_<double> projected=cameraMatrix*pMatrix*point3d;
    projected=projected/projected(2);
    return Point2d(projected(0),projected(1));
}
ostream&operator<<(ostream&out,const CloudPoint&cp)
{
    out<<cp.x<<" "<<cp.y<<" "<<cp.z<<" ";
    for (vector<int>::size_type frame=0; frame<cp.pointIndices.size(); frame++) {
        Point2d point;
        if(cp.getPointInFrame((int)frame, point))
        {
            out<<frame<<" "<<point.x<<" "<<point.y<<" ";
        }
    }
    return out;
}