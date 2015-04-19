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
    vector<int> pointIndicesInFrame;
    const vector<vector<Point2d> >& pointsInFrame;
public:
    double x;
    double y;
    double z;
    CloudPoint(double x,double y,double z,int nframe,const vector<vector<Point2d> >& pointsInFrame);
    void setPointIndex(int frame,int index);
    int  getPointIndex(int frame);
    bool getPointInFrame(int frame,Point2d&point)const;
    friend ostream&operator<<(ostream&out,const CloudPoint&cp);
};

#endif /* defined(__sfm__CloudPoint__) */
