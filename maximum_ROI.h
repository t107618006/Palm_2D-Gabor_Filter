#pragma once
#include<opencv2\core\core.hpp>
#include<opencv2\opencv.hpp>
#include<iostream>
#include<stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>  
#include<string.h>


#define HIGH 255
#define LOW 0
#define L_BASE 100

#define ROI 0
#define extraction 1


using namespace cv;
using namespace std;
using namespace cv::ml;
Mat capture_plam_image(Mat &a_timage_in, Mat a_timage_plam, Mat &strong_img, string hand);
Mat labeling_replace(Mat a_timage_in);
int labeling(Mat a_timage_in, Mat &imagelabel);