//
//  CloudPoint.cpp
//  sfm
//
//  Created by  刘骥 on 15/4/18.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#include "CloudPoint.h"
CloudPoint::CloudPoint(double x,double y,double z,int nframe,const vector<vector<Point2d> >& pointsInFrame):pointIndicesInFrame(nframe,-1),pointsInFrame(pointsInFrame)
{
    this->x=x;
    this->y=y;
    this->z=z;
}
void CloudPoint::setPointIndex(int frame, int index)
{
    pointIndicesInFrame[frame]=index;
}
int CloudPoint::getPointIndex(int frame)
{
    return pointIndicesInFrame[frame];
}
bool CloudPoint::getPointInFrame(int frame,Point2d&point)const
{
    
    int index=pointIndicesInFrame[frame];
    if (index==-1) {
        return false;
    }else{
        point=pointsInFrame[frame][index];
        return true;
    }
}
ostream&operator<<(ostream&out,const CloudPoint&cp)
{
    out<<cp.x<<" "<<cp.y<<" "<<cp.z<<" ";
    for (vector<int>::size_type frame=0; frame<cp.pointIndicesInFrame.size(); frame++) {
        Point2d point;
        if(cp.getPointInFrame((int)frame, point))
        {
            out<<frame<<" "<<point.x<<" "<<point.y<<" ";
        }
    }
    return out;
}