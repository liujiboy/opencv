//
//  CloudPoint.cpp
//  sfm
//
//  Created by  刘骥 on 15/4/18.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#include "CloudPoint.h"
CloudPoint::CloudPoint(int nframe,const vector<vector<Point2d> >& pointsInFrame):pointIndexInFrame(nframe,-1),pointsInFrame(pointsInFrame)
{
    
}
void CloudPoint::setPointIndex(int frame, int index)
{
    pointIndexInFrame[frame]=index;
}
bool CloudPoint::getPointInFrame(int frame,Point2d&point)const
{
    
    int index=pointIndexInFrame[frame];
    if (index==-1) {
        return false;
    }else{
        point=pointsInFrame[frame][index];
        return true;
    }
}