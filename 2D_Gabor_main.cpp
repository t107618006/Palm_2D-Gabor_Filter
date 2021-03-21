#include "2D-Gabor_ROI.h"
//cv::Mat mkKernel(int ks, double sig, double th, double lm, double ps);
void Gabor_CR_real_imaginary_image(Mat a_timage_in, Mat &a_tConv_real, Mat &a_tConv_imag, int theat, double SD, double CF, int sub_region_length);




double Gabor_standard_deviation(Mat a_timage_in)//評估標準差
{
	Mat ROI_SD, ROI_mean;
	meanStdDev(a_timage_in, ROI_mean, ROI_SD);
	double reg;
	if (ROI_SD.at<double>(0,0) <= 1)
	{
		reg = 1;
	}
	else if (1<ROI_SD.at<double>(0,0) && ROI_SD.at<double>(0, 0) <= 1.4)
	{
		reg = sqrt(2);
	}
	else if (1.4<ROI_SD.at<double>(0,0) && ROI_SD.at<double>(0, 0) <= 2.8)
	{
		reg = 2 * sqrt(2);
	}
	else if (ROI_SD.at<double>(0,0)>2.8)
	{
		reg = 4 * sqrt(2);
	}
	return reg;
}

double Gabor_centeral_frequency(Mat a_timage_in)
{
	double SD = Gabor_standard_deviation(a_timage_in);
	if (SD==1)
	{
		return 0;
	}
	else if(SD==sqrt(2))
	{
		return 0.12;
	}
	else if(SD==2*sqrt(2))
	{
		return 0.8;
	}
	else if(SD==4*sqrt(2))
	{
		return 2;
	}
}


double Gabor_GenvXY(Mat a_timage_in, double SD, int X, int Y, int X0, int Y0)
{
	return  1 / (2 * M_PI*pow(SD, 2))*exp(-1 * ((pow(X - X0, 2) + (pow(Y - Y0, 2))) / (2 * pow(SD, 2))));
}

Mat Lk_image(Mat a_timage_in, double theta, Point center)//產生LK線條的影像，輸入sub-region,要產生的角度只限0,30,60,90,12,150度
{
	int q=6-(theta/30);
	if (theta == 0)
		q = 0;
	double slope_arry[6] = { 0,sqrt(3) / 3,sqrt(3),90,-1 * sqrt(3),-1 * sqrt(3) / 3 };
	Mat a_timage_out(a_timage_in.rows, a_timage_in.cols, CV_8U, Scalar(0));
	//	Point center(a_timage_out.rows / 2, a_timage_out.cols / 2);
	Mat kernel_dil = getStructuringElement(MORPH_RECT, Size(5, 5));
	Mat kernel_ero = getStructuringElement(MORPH_RECT, Size(3, 3));
	double slope_k = slope_arry[q];
	if (slope_k != 90)
	{
		for (int x = 0; x <a_timage_out.cols - center.x; x++)
		{
			int y = (slope_k*(x - center.x));
		
			a_timage_out.at<uchar>(y + center.y, x + center.x) = 255;
	
		}
		for (int x = 0; x >-center.x - 1; x--)
		{
			int y = (slope_k*(x - center.x));
	//		if (y + center.y < 0)
	//			y = -center.y;
	//		if (y + center.y > a_timage_out.rows - 1)
	//			y = a_timage_out.rows - 1 - center.y;
			a_timage_out.at<uchar>(center.y + y, center.x + x) = 255;
/*int y = (slope_k*(x - center.x)) + center.y;
			if (y + center.y < 0)
				y = -center.y;
			if (y + center.y > a_timage_out.rows - 1)
				y = a_timage_out.rows - 1 - center.y;
			a_timage_out.at<uchar>(center.y + y, center.x + x) = 255;*/
		}
		dilate(a_timage_out, a_timage_out, kernel_dil, Point(2, 2), 1);
		erode(a_timage_out, a_timage_out, kernel_ero, Point(2, 2), 1);
	}
	else
	{
		for (int y = 0; y < a_timage_in.rows; y++)
		{
			a_timage_out.at<uchar>(y, center.x) = 255;
		}
		dilate(a_timage_out, a_timage_out, kernel_dil, Point(1, 1), 1);
		erode(a_timage_out, a_timage_out, kernel_ero, Point(1, 1), 1);
	}
	return a_timage_out;
}

double Gabor_Gaussian_Radon_I(Mat a_timage_in, Point center, double SD_GR, double theta)//只針對Lk線段做
{
	Mat a_tmask_Lk = Lk_image(a_timage_in, theta,center);
	double total_Genv = 0;
	double total_gray = 0;
	for (int i = 0; i < a_timage_in.rows; i++)
	{
		for (int j = 0; j < a_timage_in.cols; j++)
		{
			if(a_tmask_Lk.at<uchar>(i,j)==255)
				total_Genv += Gabor_GenvXY(a_tmask_Lk, SD_GR,j,i,center.x,center.y);
		}
	}
	for (int i = 0; i < a_timage_in.rows; i++)
	{
		for (int j = 0; j < a_timage_in.cols; j++)
		{
			if (a_tmask_Lk.at<uchar>(i, j) == 255)
				total_gray += a_timage_in.at<uchar>(i, j)*((Gabor_GenvXY(a_tmask_Lk, SD_GR, j, i, center.x, center.y) / total_Genv));
		}
	}
	return total_gray;
}

double Gabor_Dkxy(Mat a_timage_in, Point center, double SD_GR,double theta)
{
	
//	double output_arry[6];
	double min_output = 1000000;
	int angle_arry[6] = { 0,30,60,90,120,150 };
	for (int theta = 0; theta < 6; theta++)
	{
		min_output = Gabor_Gaussian_Radon_I(a_timage_in, center, SD_GR, angle_arry[theta]);
		//	output_arry[theta] = Gabor_Gaussian_Radon_I(a_timage_in, center, SD_GR, theta_arry[theta]);
	}
//	for (int i = 0; i < 6; i++)
	{
	//	if (output_arry[i] < min_output)
		{
		//	min_output = output_arry[i];
		}
	}
	return min_output;
}

int calculateOrientations(Mat a_timage_in, int theta_num)//回傳輸入圖片的最大角度
{
	Mat gradientX;
	Mat gradientY;
	Sobel(a_timage_in, gradientX, CV_32F, 1, 0, 3);
	Sobel(a_timage_in, gradientY, CV_32F, 0, 1, 3);
	// Create container element
	Mat orientation = Mat(gradientX.rows, gradientX.cols, CV_32F);
	
	vector<int> theta_vector;
	double Max_theta = 99999.9;
	int theta_output = 0;
	double theat_arry[6];
	for (int i = 0; i < theta_num; i++)
	{
		theta_vector.push_back(180 / theta_num * i);
		theat_arry[i] = 0;
	}
	// Calculate orientations of gradients --> in degrees
	// Loop over all matrix values and calculate the accompagnied orientation
	normalize(gradientX, gradientX, 1, 0, NORM_L2);
	normalize(gradientY, gradientY, 1, 0, NORM_L2);
	for (int i = 0; i < gradientX.rows; i++) {
		for (int j = 0; j < gradientX.cols; j++) {
			// Retrieve a single value
			float valueX = gradientX.at<float>(i, j);
			float valueY = gradientY.at<float>(i, j);
			// Calculate the corresponding single direction, done by applying the arctangens function
			float result = fastAtan2(valueY, valueX);
			// Store in orientation matrix element
			orientation.at<float>(i, j) = result;
		}
	}
	for (int i = 0; i < orientation.rows; i++)
	{
		for (int j = 0; j < orientation.cols; j++)
		{
			double diff_arry[6];
			if (orientation.at<float>(i, j) > 165+ DBL_EPSILON)
			{
				orientation.at<float>(i, j) = orientation.at<float>(i, j) - 180;
			}
			for (int x = 0; x < theta_num; x++)
			{
				diff_arry[x] = abs(orientation.at<float>(i, j) - theta_vector[x]);
				if (diff_arry[x] < Max_theta)
				{
					Max_theta= diff_arry[x];
					theta_output = theta_vector[x];
				}
			}
			double temple=orientation.at<float>(i, j) - theta_output;
			int position = theta_output / 30;
			if (i > 0 & i < orientation.rows - 1 & j>0 & j < orientation.cols - 1)
			{
				if (temple == 0.0)
				{
					theat_arry[int(position)]++;
				}
				else
				{
					if (temple >= 0)
					{
						theat_arry[position] += 1-(temple / 30);
					//	cout << 1 - (temple / 30) << endl;
						theat_arry[position + 1] += 1 - ((theta_output + 30 - orientation.at<float>(i, j)) / 30);
					//	cout << 1 - ((theta_output + 30 - orientation.at<float>(i, j)) / 30) << endl;
					}
					else
					{
						theat_arry[position] +=1-( abs(temple) / 30);
						//cout << 1 - (abs(temple) / 30) << endl;
						theat_arry[position - 1] += ((theta_output - orientation.at<float>(i, j)) / 30);
						//cout << ((theta_output-orientation.at<float>(i, j) ) / 30) << endl;
					}
				}
			}
			//orientation.at<float>(i, j) = theta_output;
			//theat_arry[int(position)]++;
			Max_theta = 99999;
		}
	}

/*	for (int i = 1; i < orientation.rows-1; i++)
	{
		for (int j = 1; j < orientation.cols-1; j++)
		{
			float x = orientation.at<float>(i, j);
			x /= 30;
			theat_arry[int(x)]++;
		}
	}*/
	int theta_counter = 0;
	for (int x = 0; x < theta_num; x++)
	{
		if (theat_arry[x] > theta_counter)
		{
			theta_output = x * 180/theta_num;
			theta_counter = theat_arry[x];
		}
	}
	return theta_output;
}

int calculateOrientations_radon(Mat a_timage_in)
{
	double t = getTickCount();
	Mat dst = a_timage_in.clone();
	dst.convertTo(dst, CV_32FC1);
	int angle = 180;
	Mat radon_image = Mat(dst.rows, angle, CV_32FC1);
	int center = dst.rows / 2;

	float shift0[] = { 1, 0, -center,
		0, 1, -center,
		0, 0, 1 };
	float shift1[] = { 1, 0, center,
		0, 1, center,
		0, 0, 1 };
	Mat m0 = Mat(3, 3, CV_32FC1, shift0);
	Mat m1 = Mat(3, 3, CV_32FC1, shift1);
	float *theta = new float[angle];//旋转角度
	for (int t = 0; t<angle; t++)
	{
		theta[t] = t*CV_PI / angle;
		float R[] = { cos(theta[t]), sin(theta[t]), 0,
			-sin(theta[t]), cos(theta[t]), 0,
			0, 0, 1 };
		Mat mR = Mat(3, 3, CV_32FC1, R);
		Mat rotation = m1*mR*m0;
		Mat rotated;
		warpPerspective(dst, rotated, rotation, Size(dst.rows, dst.cols), WARP_INVERSE_MAP);
		for (int j = 0; j<rotated.cols; j++)
		{
			float *p1 = radon_image.ptr<float>(j);
			for (int i = 0; i<rotated.rows; i++)
			{
				float *p2 = rotated.ptr<float>(i);
				p1[t] += p2[j];
			}
		}


	}
	cout << (double)(getTickCount() - t) / getTickFrequency() << endl;
	normalize(radon_image, radon_image, 0, 179, CV_MINMAX);
	imshow("My Radon Transform", radon_image);

	return 0;
}

int Orientation_Matrix(Mat a_timage_in,int &theta)
{
	Mat a_timage_OM(a_timage_in.rows, a_timage_in.cols, CV_32F, Scalar(0));
	double SD_GR = Gabor_standard_deviation(a_timage_in);
//	double OM_max, OM_min;
//	int output_max = 0,max_counter;
	int output_angle;
	int output_min = 100000000000, min_counter;
	Scalar image_SUM = 0;
	int theta_num = 6;//總共要做幾次角度比較	
	double output_arry[6];
	int Orientations = calculateOrientations(a_timage_in, theta_num);
	/*
	for (int theta = 0; theta < 6; theta++)
	{
		for (int i = 0; i < a_timage_in.rows; i++)
		{
			for (int j = 0; j < a_timage_in.cols; j++)
			{
				a_timage_OM.at<float>(i, j) = Gabor_Dkxy(a_timage_in, Point(16, 16), SD_GR,theta_arry[theta]);
			}
		}
		image_SUM = sum(a_timage_OM);
		output_arry[theta] = image_SUM[0];
	}
	for (int q = 0; q < 6; q++)
	{
		if (output_arry[q] < output_min)
		{
			output_min = output_arry[q];
			min_counter = q;
		}
	}
	output_angle = theta_arry[min_counter];
	theta = output_angle;*/
	return output_min;
}


void Gabor_sub_region_parameter(Mat a_timage_in,int &main_orientation,double &SD_GR,double &center_frequency, int sub_region_counter, int theta_num)
{
	
	SD_GR = Gabor_standard_deviation(a_timage_in);
	center_frequency = Gabor_centeral_frequency(a_timage_in);
//	main_orientation = Orientation_Matrix(a_timage_in,theta);
	main_orientation = calculateOrientations(a_timage_in, theta_num);
//	calculateOrientations_radon(a_timage_in);
}

double Gabor_gSD(int X0, int Y0,int SD)
{
	return (1 / (2 * M_PI*pow(SD, 2)))*exp((-1 * (pow(X0, 2) + pow(Y0, 2))) / (2 * pow(SD, 2)));
}

double Gabor_real_part(int i,int j, int theta, double SD_GR, double CF, uchar vaule)
{
	return vaule*cos(2 * M_PI*CF*(j*cos(theta*CV_PI / 180) + i*sin(theta*CV_PI / 180)));
}

double Gabor_imaginary_part(int i, int j, int theta, double SD_GR, double CF, uchar vaule)
{
	return vaule*sin(2 * M_PI*CF*(j*cos(theta*CV_PI / 180) + i*sin(theta*CV_PI / 180)));
}

/*Mat Convolution(Mat a_timage_in, Mat a_timage_mask,Mat a_timage_length)
{
	Point center(a_timage_mask.cols/2, a_timage_mask.rows/2);
	Mat a_timage_out_float(a_timage_length.rows, a_timage_length.cols, CV_32F, Scalar(0));
	Mat a_timage_out_binary(a_timage_length.rows, a_timage_length.cols, CV_8U, Scalar(0));
	double reg_gray = 0;
	for (int i = 0; i < a_timage_out_float.rows; i++)
	{
		for (int j = 0; j < a_timage_out_float.cols; j++)
		{
			Mat a_timage_reg(a_timage_in.rows, a_timage_in.cols, CV_32F, Scalar(0));
			for (int y = 0; y < a_timage_mask.rows; y++)
			{
				for (int x = 0; x < a_timage_mask.cols; x++)
				{
					reg_gray += a_timage_in.at<uchar>(i +  y, j + x)*a_timage_mask.at<float>(y,x);
					a_timage_reg.at<float>(i + y, j + x )=255;
				}
			}
			a_timage_out_float.at<float>(i, j) = reg_gray;
			if (reg_gray >= 0)
				a_timage_out_binary.at<uchar>(i, j) = 255;
			else
			{
				a_timage_out_binary.at<uchar>(i, j) = 0;
			}
			reg_gray = 0;

		}
	}
	return a_timage_out_binary;
}*/

/*double Gabor_Gaussian_Radon(Mat a_timage_in, Point center)//全sub-region做一遍
{
	double SD_GR = Gabor_standard_deviation(a_timage_in);
	double total_Genv = 0;
	double total_gray = 0;
	for (int i = 0; i < a_timage_in.rows; i++)
	{
		for (int j = 0; j < a_timage_in.cols; j++)
		{
			total_Genv += Gabor_GenvXY(a_timage_in, SD_GR, j, i, center.x, center.y);
		}
	}
	for (int i = 0; i < a_timage_in.rows; i++)
	{
		for (int j = 0; j < a_timage_in.cols; j++)
		{
			total_gray += a_timage_in.at<uchar>(i, j)*(Gabor_GenvXY(a_timage_in, SD_GR, j, i, center.x, center.y)) / total_Genv;
		}
	}
	return total_gray;
}

double Orientation_Matrix_2(Mat a_timage_in)//全sub-region做一遍
{
	Mat a_timage_out(a_timage_in.rows, a_timage_in.cols, CV_32F,Scalar(0));
	for (int i = 0; i < a_timage_in.rows; i++)
	{
		for (int j = 0; j < a_timage_in.cols; j++)
		{
			a_timage_out.at<float>(i, j) = Gabor_Gaussian_Radon(a_timage_in, Point(j, i));
		}
	}
	return 1.1;
}
*/

/*Mat convvv(Mat a_timage_in, Mat mask)
{
	Mat a_timage_out(a_timage_in.rows, a_timage_in.cols, CV_32F, Scalar(0));
	double reg = 0;
	for (int i = 0; i <mask.rows; i++)
	{
		for (int j = 0; j < mask.cols; j++)
		{
			for (int x = 0; x < a_timage_in.rows; x++)
			{
				for (int y = 0; y < a_timage_in.cols; y++)
				{
					reg += a_timage_in.at<uchar>(y, x)*mask.at<float>(i, j);
				}
			}
			if (reg >= 0)
			{
				reg = 255;
			}
			else
			{
				reg = 0;
			}
			a_timage_out.at<float>(i, j) = reg;
			reg = 0;
		}
	}
	return a_timage_out;
}*/
/*
Mat Gabor_remove_DC(Mat a_timage_palm, Mat a_timage_Gabor)
{
	Mat tmp_m, tmp_sd;
	double image_mean, SD;
	meanStdDev(a_timage_Gabor, tmp_m, tmp_sd);
	image_mean = tmp_m.at<double>(0, 0);
	SD = tmp_sd.at<double>(0, 0);
	Mat a_timage_out(a_timage_Gabor.rows, a_timage_Gabor.cols, a_timage_Gabor.type(), Scalar(0));
	Mat a_tGabor_rmDC(a_timage_Gabor.size(), CV_32F, Scalar(0));
	for (int i = 0; i < a_timage_palm.rows; i++)
	{
		for (int j = 0; j < a_timage_palm.cols; j++)
		{
			a_tGabor_rmDC.at<float>(i, j) = a_timage_Gabor.at<float>(i, j) - image_mean;
		}
	}
	filter2D(a_timage_palm, a_timage_out, CV_32F, a_tGabor_rmDC);
	return a_timage_out;
}

Mat Gabor_real_image(Mat a_timage_in, int theta, int SD, int CF)
{
	Mat test_gau;
	cv::GaussianBlur(a_timage_in,test_gau,Size(3,3),SD,SD);
	Mat a_tGabor_filter(a_timage_in.size(), CV_32F, Scalar(0));
	for (int i = 0; i < a_timage_in.rows; i++)
	{
		for (int j = 0; j < a_timage_in.cols; j++)
		{
			a_tGabor_filter.at<float>(i, j) = Gabor_real_part(i,j, theta, SD, CF, (test_gau.at<uchar>(i, j)));
		}
	}
	Mat Gabor_filter_real = Gabor_remove_DC(a_timage_in, a_tGabor_filter);
	Mat a_timage_out(a_timage_in.rows, a_timage_in.cols, CV_32F, Scalar(0));
	Mat a_timage_Border;
	double reg = 0;
	copyMakeBorder(a_timage_in, a_timage_Border, a_timage_in.rows / 2, a_timage_in.rows / 2, a_timage_in.cols / 2, a_timage_in.cols / 2, BORDER_CONSTANT, Scalar(0));
	for (int i = 0; i < a_timage_in.rows; i++)
	{
		for (int j = 0; j < a_timage_in.cols; j++)
		{
			for (int y = -a_timage_in.rows / 2; y < a_timage_in.rows / 2; y++)
			{
				for (int x = -a_timage_in.cols / 2; x < a_timage_in.cols / 2; x++)
				{
					reg = reg + a_timage_Border.at<uchar>(i + y + a_timage_in.rows / 2, j + x + a_timage_in.cols / 2)*Gabor_filter_real.at<float>(i,j);
				}
			}
			a_timage_out.at<float>(i, j) = reg;
			reg = 0;
		}
	}
	return a_timage_out;
}

Mat Gabor_imaginary_image(Mat a_timage_in, int theta, int SD, int CF)
{
	Mat test_gau;
	cv::GaussianBlur(a_timage_in, test_gau, Size(3, 3), SD, SD);
	Mat a_tGabor_filter(a_timage_in.size(), CV_32F, Scalar(0));
	for (int i = 0; i < a_timage_in.rows; i++)
	{
		for (int j = 0; j < a_timage_in.cols; j++)
		{
			a_tGabor_filter.at<float>(i, j) = Gabor_imaginary_part(i, j, theta, SD, CF, (test_gau.at<uchar>(i, j)));
		}
	}
	Mat Gabor_filter_imaginary = Gabor_remove_DC(a_timage_in, a_tGabor_filter);
	Mat a_timage_out(a_timage_in.rows, a_timage_in.cols, CV_32F, Scalar(0));
	Mat a_timage_Border;
	double reg = 0;
	copyMakeBorder(a_timage_in, a_timage_Border, a_timage_in.rows / 2, a_timage_in.rows / 2, a_timage_in.cols / 2, a_timage_in.cols / 2, BORDER_CONSTANT, Scalar(0));
	for (int i = 0; i < a_timage_in.rows; i++)
	{
		for (int j = 0; j < a_timage_in.cols; j++)
		{
			for (int y = -a_timage_in.rows / 2; y < a_timage_in.rows / 2; y++)
			{
				for (int x = -a_timage_in.cols / 2; x < a_timage_in.cols / 2; x++)
				{
					reg = reg + a_timage_Border.at<uchar>(i + y + a_timage_in.rows / 2, j + x + a_timage_in.cols / 2)*Gabor_filter_imaginary.at<float>(i, j);
				}
			}
			a_timage_out.at<float>(i, j) = reg;
			reg = 0;
		}
	}
	return a_timage_out;
}

*/

void drawHistImg(const Mat &src, Mat &dst) {
	int histSize = 256;
	float histMaxValue = 0;
	for (int i = 0; i<histSize; i++) {
		float tempValue = src.at<float>(i);
		if (histMaxValue < tempValue) {
			histMaxValue = tempValue;
		}
	}

	float scale = (0.9 * 256) / histMaxValue;
	for (int i = 0; i<histSize; i++) {
		int intensity = static_cast<int>(src.at<float>(i)*scale);
		line(dst, Point(i, 255), Point(i, 255 - intensity), Scalar(0));
	}
}

Mat cal_histogram(Mat &src)
{
	int histSize = 256;
	float range[] = { 0, 255 };
	const float* histRange = { range };
	Mat histImg;
	calcHist(&src, 1, 0, Mat(), histImg, 1, &histSize, &histRange);

	Mat showHistImg(256, 256, CV_8UC1, Scalar(255));  //把直方圖秀在一個256*256大的影像上
	drawHistImg(histImg, showHistImg);
	return showHistImg;
}

int main()
{
	char filename[100];
	char savename[100];
	char Hand_Region[100];
	int quantity_num = 100;
	int finger_count = 1;
	int finger_image_num = 1;
	
	string sdf;
	int spectrum = 850;
	
#if ROI
	string spectrum_select = "850nm";
	int A = 8;
	stringstream temp;
	temp << A;
	string size_num;
	temp >> size_num;
	string original_path = "C:\\Image\\Plam_CASIA\\Image\\850nm\\original\\";
	string ROI_path = "C:\\Image\\Plam_CASIA\\Image\\850nm\\maximum\\ROI\\";
	String log_path = "C:\\Image\\Plam_CASIA\\Image\\850nm\\maximum\\Log\\";

//	save_real = save_real + "\\Gabor\\new_angle\\CLAHE_9_" + size_num + "\\real\\";
	//save_imaginary = save_imaginary + "\\Gabor\\new_angle\\CLAHE_9_" + size_num + "\\imaginary\\";

	vector<String> image;
	glob(original_path, image, false);

	//	string templet = "C:\\Image\\Plam_CASIA\\Good_ROI\\norm_old\\850nm\\POHE\\";
//#pragma omp parallel for
	for (size_t i = 0; i < image.size(); i++)
	{
		//	Mat g_utimage_in = imread("Barbara.png", 0);
		Mat g_utimage_in = imread(image[i], 0);
		/*size_t sdf = image[i].find("\\ROI\\");
		string name = image[i].substr(sdf + 7);*/
		size_t sdf = image[i].find("\\original\\");
		string name = image[i].substr(sdf + 4 + 6, 12);
		string savename = ROI_path + "Maximum_ROI_" + name + ".png";
		string Hand_Region= log_path + "Maximum_ROI_Log_" + name + ".png";
		Mat che = imread(savename, 0);
		if (che.empty() == 0)
		{
			continue;
			waitKey(10);
		}
		cout << image[i] << endl;

		size_t zxc = name.find("_850_");
		string hand_check = name.substr(zxc-1,1);
	Mat g_utimage_OTSU;
	double thresh = 0.85;
	Mat g_utimage_gaussian;
	GaussianBlur(g_utimage_in, g_utimage_gaussian, Size(5, 5), 0, 0);
	threshold(g_utimage_gaussian, g_utimage_OTSU, thresh, 255, CV_THRESH_OTSU);
	for (int i = 0; i < g_utimage_OTSU.rows; i++)
	{
		for (int j = 0; j < g_utimage_OTSU.cols; j++)
		{
			if (g_utimage_OTSU.at<uchar>(i, j) == 255)
			{
				g_utimage_OTSU.at<uchar>(i, j) = 0;
			}
			else if (g_utimage_OTSU.at<uchar>(i, j) == 0)
			{
				g_utimage_OTSU.at<uchar>(i, j) = 255;
			}
		}
	}
	Mat g_utimage_label;
	medianBlur(g_utimage_label, g_utimage_label, 9);
	labeling(g_utimage_OTSU, g_utimage_label);
	g_utimage_label = labeling_replace(g_utimage_label);
	//		medianBlur(g_utimage_label, g_utimage_label, 9);
			//	imwrite("007_l_940_03_binary.png", g_utimage_label);
	Mat g_utstrong_img;
	

	Mat g_utimage_ROI = capture_plam_image(g_utimage_label, g_utimage_in, g_utstrong_img, hand_check);
	hconcat(g_utimage_label, g_utimage_in, g_utimage_label);
	cvtColor(g_utimage_label, g_utimage_label, CV_GRAY2BGR);
	vconcat(g_utstrong_img, g_utimage_label, g_utimage_label);
	imwrite(savename, g_utimage_ROI);
	imwrite(Hand_Region, g_utimage_label);
#endif // ROI
#include "POHE.h"
#if extraction
	string spectrum_select = "850nm";
	int A = 32;
	stringstream temp;
	temp << A;
	string size_num;
	temp  >> size_num;
	string ROI_path = "C:\\Image\\Plam_CASIA\\\Image\\norm_old\\"+ spectrum_select;
	
	//ROI_path = "C:\\Users\\raychen\\Downloads\\vein";
	string save_real = "C:\\Image\\Plam_CASIA\\Good_ROI\\norm_old\\" + spectrum_select;
	string save_imaginary = "C:\\Image\\Plam_CASIA\\Good_ROI\\norm_old\\" + spectrum_select;
	
	save_real = save_real + "\\Gabor\\new_angle\\CLAHE_9_" + size_num + "\\real\\";
	save_imaginary = save_imaginary + "\\Gabor\\new_angle\\CLAHE_9_" + size_num + "\\imaginary\\";

	string user = "187\\";
	ROI_path = "C:\\Image\\templet\\" + user;
	save_real = "C:\\Image\\templet\\"+ user;
	save_imaginary = "C:\\Image\\templet\\"+ user;

	vector<String> image;
	glob(ROI_path, image, false);

//	string templet = "C:\\Image\\Plam_CASIA\\Good_ROI\\norm_old\\850nm\\POHE\\";
	
	for (size_t i = 0; i < image.size(); i++)
	{
		//	Mat g_utimage_in = imread("Barbara.png", 0);
		
		/*size_t sdf = image[i].find("\\ROI\\");
		string name = image[i].substr(sdf + 7);*/
		size_t sdf = image[i].find("Norm_ROI");
		string name = image[i].substr(sdf + 9);
		Mat che = imread(save_imaginary + "norm_gabor_imaginary_" + name, 0);
		if (che.empty() == 0)
		{
			//continue;
			waitKey(10);
		}
		cout << image[i] << endl;
		Mat g_utimage_in = imread(image[i], 0);
		
		resize(g_utimage_in, g_utimage_in, Size(160, 160));
		//	int A = pow(2, Z+1);
		
		int sub_region_length = A;
			//sub_region_length = 32;
		if (sub_region_length > g_utimage_in.rows)
			return 0;
		Mat g_utimage_sub_region;
		Mat g_utimage_in_CLAHE;
		Mat g_utconvolution_real;
		Mat g_utconvolution_imaginary;
		Mat g_utreal_horizontal;//實部水平累加
		Mat g_utreal_vertical;//實部垂直累加
		Mat g_utimaginary_horizontal;//虛部水平累加
		Mat g_utimaginary_vertical;//虛部垂直累加
		int sub_region_counter = g_utimage_in.rows / sub_region_length;
		double SD, center_frequency;
		int main_orientation;
		Mat sum;
		Mat sqsum;
		Mat g_utimage_POHE;
	//	copyMakeBorder(g_utimage_in, g_utimage_in, 1, 1, 1, 1, BORDER_REPLICATE);
	//	Mat one_ori_h = cal_histogram(g_utimage_in);
	//	Mat das;
	//	equalizeHist(g_utimage_in, das);
	//	Mat one_dsa_h=cal_histogram(das);
	
	//	Mat one_POHE_h = cal_histogram(g_utimage_POHE);

		Ptr<CLAHE> clahe = createCLAHE();
		clahe->setClipLimit(9);
		clahe->apply(g_utimage_in, g_utimage_POHE);
	//	POHE2013(g_utimage_in, g_utimage_POHE, Size(g_utimage_in.cols -1, g_utimage_in.rows -1), sum, sqsum);
	//	imwrite(templet + "norm_ROI_POHE_" + name, g_utimage_POHE);

	//	g_utimage_in = g_utimage_in(cv::Rect(1, 1, g_utimage_in.cols - 2, g_utimage_in.rows - 2));
	//	g_utimage_POHE = g_utimage_POHE(cv::Rect(1, 1, g_utimage_POHE.cols - 2, g_utimage_POHE.rows - 2));

		g_utimage_in_CLAHE = g_utimage_POHE.clone();
	//	Mat one_CLAHE_h= cal_histogram(g_utimage_in_CLAHE);
		//					imwrite("Maximum_ROI\\007_r_850_05_ROI_CLAHE.png", g_utimage_in_CLAHE);
		Mat g_utsub_region_getedge;
		Mat reg;
		Mat g_utimage_parameter(sub_region_counter, sub_region_counter, CV_32FC3, Scalar(0));//存參數圖片，順序為角度.SD.CF
		Mat g_utimage_theta(sub_region_counter, sub_region_counter, CV_8U, Scalar(0));//同g_utimage_parameter，只是儲存theta
		int theta_num = 6;//180度內取幾次做角度
		int c = 0, b = 0;
		//			copyMakeBorder(g_utimage_in_CLAHE, g_utsub_region_getedge, sub_region_length / 2, sub_region_length / 2, sub_region_length / 2, sub_region_length / 2, BORDER_CONSTANT, Scalar(0));//將圖片擴大好做convolution
		for (int i = 0; i < sub_region_counter; i++)
		{
			for (int j = 0; j < sub_region_counter; j++)
			{
				reg = g_utimage_in_CLAHE(Rect(j*sub_region_length, i*sub_region_length, sub_region_length, sub_region_length));//這裡的ROI是指子區域(sub-region)
				g_utimage_sub_region = reg.clone();
				//		copyMakeBorder(g_utimage_sub_region, g_utsub_region_getedge, sub_region_length / 2, sub_region_length / 2, sub_region_length / 2, sub_region_length / 2, BORDER_CONSTANT, Scalar(0));//將圖片擴大好做convolution
				Gabor_sub_region_parameter(g_utimage_sub_region, main_orientation, SD, center_frequency, sub_region_counter, theta_num);
				g_utimage_parameter.at<Vec3f>(b, c)[0] = main_orientation;
				g_utimage_parameter.at<Vec3f>(b, c)[1] = SD;
				g_utimage_parameter.at<Vec3f>(b, c)[2] = center_frequency;
				
				Gabor_CR_real_imaginary_image(g_utimage_sub_region, g_utconvolution_real, g_utconvolution_imaginary, main_orientation, SD, center_frequency, sub_region_length);
				
				if (j == 0)
				{
					g_utreal_horizontal = g_utconvolution_real.clone();
					g_utimaginary_horizontal = g_utconvolution_imaginary.clone();
				}
				else
				{
					hconcat(g_utreal_horizontal, g_utconvolution_real, g_utreal_horizontal);
					hconcat(g_utimaginary_horizontal, g_utconvolution_imaginary, g_utimaginary_horizontal);
				}
				
				c++;
			}
			b++;
			c = 0;
			
			if (i == 0)
			{
				g_utreal_vertical = g_utreal_horizontal.clone();
				g_utimaginary_vertical = g_utimaginary_horizontal.clone();
			}
			else
			{
				vconcat(g_utreal_vertical, g_utreal_horizontal, g_utreal_vertical);
				vconcat(g_utimaginary_vertical, g_utimaginary_horizontal, g_utimaginary_vertical);
			}
		}
		
		medianBlur(g_utreal_vertical, g_utreal_vertical, 5);
		medianBlur(g_utimaginary_vertical, g_utimaginary_vertical, 5);
		imwrite(save_real + "norm_gabor_real_" + name, g_utreal_vertical);
		imwrite(save_imaginary + "norm_gabor_imaginary_" + name, g_utimaginary_vertical);
		
#endif // extraction




	}
	return 0;
}
//			}
//		}
//	}
//	
//	
//}

void Gabor_Kernel(Mat &a_tGabor_real,Mat &a_tGabor_imaginary,int Theta,double SD,double CF,int sub_region_length)//葛伯慮波核心產生器
{
	int gabor_size = sub_region_length;
	gabor_size = 5;
	int x, y;
	double xtmp, ytmp, tmp1, tmp2, tmp3;
	double re, im;
	a_tGabor_real.create(gabor_size, gabor_size, CV_32F);
	a_tGabor_imaginary.create(gabor_size, gabor_size, CV_32F);
	double th = Theta*CV_PI / 180;
	for (int i = 0; i < gabor_size; i++) //定义模版大小  
	{
		for (int j = 0; j < gabor_size; j++)
		{
			x = j - gabor_size / 2;
			y = i - gabor_size / 2;

			xtmp = (double)x*cos(th) + (double)y*sin(th);
			ytmp = (double)y*cos(th) - (double)x*sin(th);

			tmp1 = exp(-(pow(xtmp, 2) + pow(ytmp, 2)) / (2 * pow(SD, 2)));
			tmp2 = cos(2 * CV_PI*xtmp*CF*CV_PI / 180);
			tmp3 = sin(2 * CV_PI*xtmp*CF*CV_PI / 180);

			re = tmp1*tmp2;
			im = tmp1*tmp3;

			a_tGabor_real.at<float>(i, j) = re;
			a_tGabor_imaginary.at<float>(i, j) = im;
			//printf("%f,%f\n",re,im);  
		}
	}
}

Mat Gabor_remove_DC( Mat a_timage_Gabor)
{
	Mat tmp_m, tmp_sd;
/*	double image_mean, SD;
	meanStdDev(a_timage_Gabor, tmp_m, tmp_sd);
	image_mean = tmp_m.at<double>(0, 0);
	SD = tmp_sd.at<double>(0, 0);*/
	double reg = 0;
	for (int i = 0; i < a_timage_Gabor.rows; i++)
	{
		for (int j = 0; j < a_timage_Gabor.cols; j++)
		{
			reg += a_timage_Gabor.at<float>(i, j);
		}
	}
	double image_mean = reg / pow(a_timage_Gabor.rows, 2);
	Mat a_timage_out(a_timage_Gabor.rows, a_timage_Gabor.cols, a_timage_Gabor.type(), Scalar(0));
	for (int i = 0; i < a_timage_Gabor.rows; i++)
	{
		for (int j = 0; j < a_timage_Gabor.cols; j++)
		{
			a_timage_out.at<float>(i, j) = a_timage_Gabor.at<float>(i, j) - image_mean;
		}
	}
	return a_timage_out;
}


int km()
{
	const int MAX_CLUSTERS = 5;
	Scalar colorTab[] =
	{
		Scalar(0, 0, 255),
		Scalar(0,255,0),
		Scalar(255,100,100),
		Scalar(255,0,255),
		Scalar(0,255,255)
	};

	Mat img(500, 500, CV_8UC3);
	RNG rng(12345);

	for (;;)
	{
		int k, clusterCount = rng.uniform(2, MAX_CLUSTERS + 1);
		int i, sampleCount = rng.uniform(1, 1001);
		Mat points(sampleCount, 1, CV_32FC2), labels;

		clusterCount = MIN(clusterCount, sampleCount);
		std::vector<Point2f> centers;

		/* generate random sample from multigaussian distribution */
		for (k = 0; k < clusterCount; k++)
		{
			Point center;
			center.x = rng.uniform(0, img.cols);
			center.y = rng.uniform(0, img.rows);
			Mat pointChunk = points.rowRange(k*sampleCount / clusterCount,
				k == clusterCount - 1 ? sampleCount :
				(k + 1)*sampleCount / clusterCount);
			rng.fill(pointChunk, RNG::NORMAL, Scalar(center.x, center.y), Scalar(img.cols*0.05, img.rows*0.05));
		}

		randShuffle(points, 1, &rng);

		double compactness = kmeans(points, clusterCount, labels,
			TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 10, 1.0),
			3, KMEANS_PP_CENTERS, centers);

		img = Scalar::all(0);

		for (i = 0; i < sampleCount; i++)
		{
			int clusterIdx = labels.at<int>(i);
			Point ipt = points.at<Point2f>(i);
			circle(img, ipt, 2, colorTab[clusterIdx], FILLED, LINE_AA);
		}
		for (i = 0; i < (int)centers.size(); ++i)
		{
			Point2f c = centers[i];
			circle(img, c, 40, colorTab[i], 1, LINE_AA);
		}
		cout << "Compactness: " << compactness << endl;

	//	imshow("clusters", img);

		char key = (char)waitKey();
		if (key == 27 || key == 'q' || key == 'Q') // 'ESC'
			break;
	}
	return 0;
}

void Gabor_CR_real_imaginary_image(Mat a_timage_in,Mat &a_tConv_real, Mat &a_tConv_imag,int theat, double SD, double CF, int sub_region_length)
{
	Mat a_tGabor_real;
	Mat a_tGaobr_imag;
	Gabor_Kernel(a_tGabor_real, a_tGaobr_imag, theat, SD, CF, sub_region_length);
	a_tGabor_real = Gabor_remove_DC(a_tGabor_real);
	a_tGaobr_imag = Gabor_remove_DC(a_tGaobr_imag);
	filter2D(a_timage_in, a_tConv_real, CV_32F, a_tGabor_real);
	filter2D(a_timage_in, a_tConv_imag, CV_32F, a_tGaobr_imag);
	Mat a_real = a_tConv_real.clone();
	Mat a_imag = a_tConv_imag.clone();
/*	Mat temp_a;
	Mat centers;

	double compactness = kmeans(a_tConv_real, 2, temp_a,
		TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 10, 1.0),
		3, KMEANS_PP_CENTERS, centers);
//	km();
	double compactness2 = kmeans(a_tConv_imag, 2, a_imag,
		TermCriteria(TermCriteria::EPS + TermCriteria::COUNT, 10, 1.0),
		3, KMEANS_PP_CENTERS, centers);*/
	for (int i = 0; i < a_tConv_real.rows; i++)
	{
		for (int j = 0;j < a_tConv_real.cols; j++)
		{
			if (a_tConv_imag.at<float>(i, j) >=0)//>=0
			{
				a_tConv_imag.at<float>(i, j) = 255;
			}
			else
			{
				a_tConv_imag.at<float>(i, j) = 0;
			}
			if (a_tConv_real.at<float>(i, j)>=0)//>=0
			{
				a_tConv_real.at<float>(i, j) = 255;
			}
			else
			{
				a_tConv_real.at<float>(i, j) = 0;
			}
		}
	}

}
