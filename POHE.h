#pragma once
#include"2D-Gabor_ROI.h"
bool POHE2013(const cv::Mat &src, cv::Mat &dst, const cv::Size blockSize, const cv::Mat &sum, const cv::Mat &sqsum);