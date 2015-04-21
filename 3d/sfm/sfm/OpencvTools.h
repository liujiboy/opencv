//
//  OpencvTools.h
//  sfm
//
//  Created by  刘骥 on 15/4/20.
//  Copyright (c) 2015年  刘骥. All rights reserved.
//

#ifndef __sfm__OpencvTools__
#define __sfm__OpencvTools__

#include <stdio.h>
#include <string>
#include <opencv2/opencv.hpp>
using namespace cv;
using namespace std;
namespace opencvTool{
    void showImage(const string&name,const Mat&img,int flags = WINDOW_NORMAL);
}
#endif /* defined(__sfm__OpencvTools__) */
