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
#include "opencv2/opencv.hpp"
using namespace cv;
using namespace std;
class CloudPoint {
private:
    double x;
    double y;
    double z;
    vector<int> pointIndexInFrame;
    const vector<vector<Point2d> >& pointsInFrame;
public:
    CloudPoint(int nframe,const vector<vector<Point2d> >& pointsInFrame);
    void setPointIndex(int frame,int index);
    bool getPointInFrame(int frame,Point2d&point)const;
};

#endif /* defined(__sfm__CloudPoint__) */
