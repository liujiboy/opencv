//
//  Cloud.h
//  sfm
//
//  Created by  刘骥 on 15/4/20.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#ifndef __sfm__Cloud__
#define __sfm__Cloud__

#include <iostream>
#include <fstream>
#include <vector>
#include "CloudPoint.h"
#include <opencv2/opencv.hpp>
#include "OpencvTools.h"
using namespace cv;
using namespace std;
using namespace opencvTool;
class Cloud{
private:
    vector<CloudPoint> cloudPoints;
    const Mat&cameraMatrix;
    const vector<Mat_<double> >&pmatrices;
    int nframes;
public:
    Cloud(const Mat&cameraMatrix,const vector<Mat_<double> >&pmatrices,int nframes);
    void addPoint(const CloudPoint&cp);
    void reprojectError(int frameNum);
    void plotCloudPoints(const vector<Mat>&images);
    void saveCloudPoints(const char* fileName);
    void saveCloudPointsToPly(const vector<Mat>&images,const char* fileName);
    void findKnownPoints(int frame,int prevFrame,const vector<DMatch>&match,bool *known,Mat_<double>&points);
    int getPointSize();
    CloudPoint&getPoint(int i);
};
#endif /* defined(__sfm__Cloud__) */
