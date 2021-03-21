#include "2D-Gabor_ROI.h"
int ncon(int p[9])
{
	int	i, i1, i2;
	int q[9];
	int n = 0;

	for (i = 0; i < 9; i++) {
		if ((p[i] == 1) || (p[i] == -1)) q[i] = 0;
		else q[i] = 1;
	}
	for (i = 1; i < 9; i += 2) {
		i1 = i + 1;
		i2 = i + 2;
		if (i2 == 9) i2 = 1;
		n = n + q[i] - q[i] * q[i1] * q[i2];
	}
	return n;
}


Mat thinning(Mat a_timage_in)
{
	Mat a_timage_out = a_timage_in.clone();//(a_timage_in.rows, a_timage_in.cols, CV_8U);
	int flg = 1;
	int i, j, k, n;
	int p[9];	/* �ϧ�:1�I��:0�A�I���Ը�:-1 */

				/*for (i = 0; i < a_timage_in.rows; i++)
				for (j = 0; j < a_timage_in.cols; j++)
				a_timage_out.at<uchar>(i, j) = a_timage_in.at<uchar>(i, j);*/
	while (flg != 0) {
		flg = 0;
		for (i = 1; i <a_timage_in.rows - 1; i++) {
			for (j = 1; j < a_timage_in.cols - 1; j++) {
				p[0] = a_timage_out.at<uchar>(i, j);
				p[1] = a_timage_out.at<uchar>(i, j + 1);
				p[2] = a_timage_out.at<uchar>(i - 1, j + 1);
				p[3] = a_timage_out.at<uchar>(i - 1, j);
				p[4] = a_timage_out.at<uchar>(i - 1, j - 1);
				p[5] = a_timage_out.at<uchar>(i, j - 1);
				p[6] = a_timage_out.at<uchar>(i + 1, j - 1);
				p[7] = a_timage_out.at<uchar>(i + 1, j);
				p[8] = a_timage_out.at<uchar>(i + 1, j + 1);
				for (k = 0; k < 9; k++) {
					if (p[k] == 255) p[k] = 1;
					else if (p[k] == 0) p[k] = 0;
					else                   p[k] = -1;
				}
				/* ����1:�ϧΪ��@���� */
				if (p[0] != 1) continue;
				/*����2:�Ҭɹ���(4�ӾF�񹳯���1�ӥH�W�O�I��)*/
				if (p[1] * p[3] * p[5] * p[7] != 0) continue;
				/*����3:�O�d���I(8�ӾF�񹳯���2�ӥH�W�O�ϧ�)*/
				n = 0;
				for (k = 1; k < 9; k++) if (p[k] != 0) n++;
				if (n < 2) continue;
				/*����4:�O�d�W���I(8�ӾF�񹳯���1�ӥH�W�ϧ�)*/
				n = 0;
				for (k = 1; k < 9; k++) if (p[k] == 1) n++;
				if (n < 1) continue;
				/*����5:�O�d�s����(8�ӳs���Ƭ�1)*/
				if (ncon(p) != 1) continue;
				/*����6:�u���e�׬�2�ɡA�u�h������(8�ӾF�񹳯��O�_�ëD�����O-1�A�p�G�O-1�ɡA��Ȭ�0�A8�ӳs
				���Ƭ�1)*/
				n = 0;
				for (k = 1; k < 9; k++) {
					if (p[k] != -1) n++;
					else if (p[k] == -1) {
						p[k] = 0;
						if (ncon(p) == 1) n++;
						p[k] = -1;
					}
				}
				if (n < 8) continue;
				/*����1-6���������ɪ��d����H */
				a_timage_out.at<uchar>(i, j) = 128;
				flg++;
			}
		}
		for (i = 1; i <a_timage_in.rows - 1; i++)
			for (j = 1; j < a_timage_in.cols - 1; j++)
				if (a_timage_out.at<uchar>(i, j) == 128) a_timage_out.at<uchar>(i, j) = 0;
	}
	return a_timage_out;
}
int Track_local_point(Point &P_in, Point &P_old, Mat a_timage_in)//�M�w�U�@�ӭn�p��ڦ��Z�����I
{
	Mat a_image_ROI = a_timage_in(Rect(P_in.x - 1, P_in.y - 1, 3, 3));
	vector<int> P_x;
	vector<int> P_y;
	vector<int> P_Judge_x;//�i�઺�U�@�I�}�Cx
	vector<int> P_Judge_y;//�i�઺�U�@�I�}�Cy
	Point P_corner[4];
	int a = 0;
	for (int i = -1; i < 2; i += 2)//��X3*3���|�Ө��y��
	{
		for (int j = -1; j < 2; j += 2)
		{
			P_corner[a].x = P_in.x + i;
			P_corner[a].y = P_in.y + j;
			a++;
		}
	}
	Point P_new(0, 0);
	for (int i = 0; i < a_image_ROI.rows; i++)//�έp�`�@���h�֭ӥi�઺�U�@�I
	{
		for (int j = 0; j < a_image_ROI.cols; j++)
		{
			if (a_image_ROI.at<uchar>(i, j) == 255)
			{
				P_x.push_back(j - 1 + P_in.x);
				P_y.push_back(i - 1 + P_in.y);
			}
		}
	}

	for (int x = 0; x < P_x.size(); x++)//�⤤���I�M�e�@�I�簣�i��o�U�@�I
	{
		if (P_x[x] != P_in.x || P_y[x] != P_in.y)
			if (P_x[x] != P_old.x || P_y[x] != P_old.y)
			{
				P_Judge_x.push_back(P_x[x]);
				P_Judge_y.push_back(P_y[x]);
			}
	}

	if (P_Judge_x.size() > 1)
	{
		int watch = 0;
		for (int i = 0; i < P_Judge_x.size(); i++)//�簣�b�e�@�Ӱl���I���䪺���i���I
		{
			for (int j = 0; j < 4; j++)
			{
				if ((abs(P_Judge_x[i] - P_old.x) + abs(P_Judge_y[i] - P_old.y)) == 1)
				{
					//P_new.x = P_Judge_x[i];
					//P_new.y = P_Judge_y[i];
					P_Judge_x.erase(P_Judge_x.begin() + i);
					P_Judge_y.erase(P_Judge_y.begin() + i);
					watch++;
				}
				if (P_Judge_x.size() < 2)
					break;
			}
		}
		if (watch > 0)
		{
			P_new.x = P_Judge_x[0];
			P_new.y = P_Judge_y[0];
		}
		else if (P_Judge_x.size() != 0)
		{
			for (int i = 0; i < P_Judge_x.size(); i++)//�q���l���I�̻������I��ܤU�@�I
			{
				for (int j = 0; j < 4; j++)
				{
					if (P_Judge_x[i] == P_corner[j].x&&P_Judge_y[i] == P_corner[j].y)
					{
						P_new.x = P_Judge_x[i];
						P_new.y = P_Judge_y[i];
						break;
					}
				}
			}
		}

	}
	else if (P_Judge_x.size() != 0)//�u���@�ӥi���I
	{
		P_new.x = P_Judge_x[0];
		P_new.y = P_Judge_y[0];
	}
	P_Judge_x.clear();
	P_Judge_y.clear();
	P_x.clear();
	P_y.clear();
	vector<int>().swap(P_x);
	vector<int>().swap(P_y);
	vector<int>().swap(P_Judge_x);
	vector<int>().swap(P_Judge_y);
	if (P_new.x == 0 && P_new.y == 0)//�����S���U�@�I
		return-1;
	else
	{
		P_old = P_in;
		P_in = P_new;
		return 1;
	}
}
double Euclidean_distance(Point P_ref, Point detection)//�p��ڦ��Z���AP_ref�����I
{
	double d = 0;
	d = sqrt(pow((P_ref.x - detection.x), 2) + pow((P_ref.y - detection.y), 2));
	return d;
}
int wave_peak_valley(Mat a_timage_in, Point &check_point, int check_length)//�ˬd�o���I�O�_���i�p�Ϊi���A�i�p�^��1�A�i���^��2�A�����O�^��0
{
	double center = a_timage_in.at<float>(0, check_point.x);
	int counter_peak = 0;
	int counter_valley = 0;
	int reg = 0;

	for (int i = check_point.x + check_length / 2 * -1; i < check_point.x + check_length / 2; i++)
	{
		reg = i;
		if (i < 0)	continue;
		if (a_timage_in.at<float>(0, i) <= center)
			//		if (a_timage_in.at<float>(0, i + check_point.x) <= center)
		{
			counter_peak++;
		}
		if (a_timage_in.at<float>(0, i) >= center)
		{
			counter_valley++;
		}
		i = reg;
	}
	if (counter_peak == check_length)
	{
		return 1;
	}
	else if (counter_valley == check_length)
	{
		return 2;
	}
	else
	{
		return 0;
	}
}

void find_start_track_point_main(Mat &a_timage_canny, Point &P_ref, Point &detection, int head, int left)//��M�̤W���̥k�䪺�I��@��l��M�I
{
	int watchdog = 0;
	int count = 2;
	while (watchdog == 0)
	{
		for (int i = 0; i < a_timage_canny.cols / 3 + 50; i++)
		{
			int j = P_ref.x - count;
			if (a_timage_canny.at<uchar>(i, j) == 255)
			{
				detection.x = j;
				detection.y = i;
				watchdog++;
				break;
			}
		}
		count++;
	}


	for (int i = detection.y + 1; i < detection.y + 5; i++)//�N����I�W���U��M���קK����
	{
		a_timage_canny.at<uchar>(i, detection.x) = 0;
	}
}

void track_plamedge2image_main(int &track_check, int &watchdog, Point &detection, Point &detection_old, Mat &a_timage_canny, Mat &euc_distance, Mat &reg, Point &P_ref, int &euc_counter)
//�p���Ӥ�x���ڦ��Z���A�æs�J�Ϥ��}�C��
{
//	fstream file;
//	file.open("Fpsfinal.csv", ios::out);
	while (track_check == 1 && watchdog == 0)
	{
		track_check = Track_local_point(detection, detection_old, a_timage_canny);
		if (track_check == 1)
		{
			euc_distance.at<float>(0, euc_counter) = Euclidean_distance(P_ref, detection);
			euc_distance.at<float>(1, euc_counter) = detection.x;
			euc_distance.at<float>(2, euc_counter) = detection.y;
	//		file << euc_distance.at<float>(0, euc_counter) << endl;
			euc_counter++;
			reg.at<uchar>(detection) = 255;
		}
		if (euc_counter == 3038)
		{
			int fdadadde = 0;
			fdadadde++;
		}

	}
//	file.close();
}

void near_point_find_main(Mat a_timage_canny, int detection_near_num, Point &detection, Point &detection_old, int &watchdog, int &track_check, Mat &euc_distance, int &euc_counter, Mat &reg, Point P_ref)
//��M���_�I�̪񪺨����I�A�ϥΤ��_�I�H�U���I�M��A�íp��X�U�I���ڦ��Z���A�̪��I�Y�P�_�����U�@�I�C�åB�]�t�F����ڦ��Z������M����
{
	Mat near_euc_length(3, 2500, CV_32F);
	int near_x = 0;
	int nead_y = 0;
	int near_counter = 0;
	float near_euc_num = 500000.0;
	int near_euc_position = 0;
	int a;//�N��i
	Point ori_detection = detection;
	for (int i = detection.y + 1; i < detection.y + detection_near_num; i++)
	{
		for (int j = detection.x - detection_near_num / 2; j < detection.x + detection_near_num / 2; j++)
		{
			a = i;
			if (i >= a_timage_canny.rows)
				a = a_timage_canny.rows - 1;
			detection_old.x = j;
			detection_old.y = a;
			if (a_timage_canny.at<uchar>(detection_old) == 255 && detection != detection_old)
			{
				near_euc_length.at<float>(0, near_counter) = Euclidean_distance(detection, detection_old);
				near_euc_length.at<float>(1, near_counter) = detection_old.x;
				near_euc_length.at<float>(2, near_counter) = detection_old.y;
				if (near_euc_length.at<float>(0, near_counter) <= near_euc_num)
				{
					near_euc_num = near_euc_length.at<float>(0, near_counter);
					near_euc_position = near_counter;
				}
				near_counter++;
			}
		}
	}
	detection.x = near_euc_length.at<float>(1, near_euc_position);
	detection.y = near_euc_length.at<float>(2, near_euc_position);
	if (detection.x > 0 && detection.y > 0)
	{
		for (int i = detection.y - 2; i < detection.y + 2; i++)
		{
			for (int j = detection.x + 1; j < detection.x + 4; j++)
			{
				a_timage_canny.at<uchar>(i, j) = 0;
			}
		}
		watchdog = 0;
		track_check = 1;
		detection_old = detection;
		track_plamedge2image_main(track_check, watchdog, detection, detection_old, a_timage_canny, euc_distance, reg, P_ref, euc_counter);
	}
}

void find_peak_vally_point_main(int &check_length, int &euc_counter, Point &P_reg, Mat &euc_distance, float *peak_x, float *peak_y, int &peak_counter, float *valley_x, float *valley_y, int &valley_counter)
//��X�ڦ��Z�����A��x�ڦ��Z�����̰��M�̧C�I(RDF�Ϲ�)
{
	int not_find_near_distance = 30;
	for (int z = check_length / 2; z < euc_counter - check_length / 2; z++)
	{
		P_reg.x = z;
		int receive = wave_peak_valley(euc_distance, P_reg, check_length);
		if (receive == 1)
		{
			peak_x[peak_counter] = euc_distance.at<float>(1, z);
			peak_y[peak_counter] = euc_distance.at<float>(2, z);
			peak_counter++;
			if (abs(peak_x[peak_counter - 1] - peak_x[peak_counter - 2]) < not_find_near_distance && abs(peak_y[peak_counter - 1] - peak_y[peak_counter - 2]) < not_find_near_distance)//�קK��쪺�̧C�I���۾F�I
			{
				peak_x[peak_counter] = 0;
				peak_y[peak_counter] = 0;
				peak_counter--;
			}
		}
		if (receive == 2)
		{
			valley_x[valley_counter] = euc_distance.at<float>(1, z);
			valley_y[valley_counter] = euc_distance.at<float>(2, z);
			valley_counter++;
			if (abs(valley_x[valley_counter - 1] - valley_x[valley_counter - 2]) < not_find_near_distance && abs(valley_y[valley_counter - 1] - valley_y[valley_counter - 2]) < not_find_near_distance)//�קK��쪺�̧C�I���۾F�I
			{
				valley_x[valley_counter] = 0;
				valley_y[valley_counter] = 0;
				valley_counter--;
			}
		}
	}
}

static int aoiGravityCenter(Mat in, Point &center)//��X�Ϥ�����
{

	IplImage* img = new IplImage(in);
	double m00, m10, m01;
	CvMoments moment;
	cvMoments(img, &moment, 1);
	m00 = cvGetSpatialMoment(&moment, 0, 0);
	if (m00 == 0)
		return 1;
	m10 = cvGetSpatialMoment(&moment, 1, 0);
	m01 = cvGetSpatialMoment(&moment, 0, 1);
	center.x = (int)(m10 / m00);
	center.y = (int)(m01 / m00);
	return 0;
}



	
Mat capture_plam_image(Mat &a_timage_in, Mat a_timage_plam, Mat &strong_img, string hand)/*��J���b���G��ƹϤ��A�æ^��ROI�Ϥ��Ain=�G��ƹϤ�,plam=��x���,hand_check=��J�{�b�O�����٬O�k��Ϥ�,strong_img�i�N�Ϥ��^�ǥD�{�����Ȧs��*/
{
	Mat kernal(5, 5, CV_8U);
	erode(a_timage_in, a_timage_in, kernal);
	erode(a_timage_in, a_timage_in, kernal);
	erode(a_timage_in, a_timage_in, kernal);
	erode(a_timage_in, a_timage_in, kernal);
	dilate(a_timage_in, a_timage_in, kernal);
	dilate(a_timage_in, a_timage_in, kernal);
	dilate(a_timage_in, a_timage_in, kernal);
	dilate(a_timage_in, a_timage_in, kernal);

	float head = 1000000, leg = 0, right = 0, left = 1000000;//�ϥγo�|�ӰѼƧ�X�G��ƹϤ����@�ӥ���Τj�p
	for (int i = 0; i < a_timage_in.rows; i++)
	{
		for (int j = 0; j < a_timage_in.cols; j++)
		{
			if (a_timage_in.at<uchar>(i, j) == 255)
			{
				if (i < head)
				{
					head = i;
				}
				if (i > leg)
				{
					leg = i;
				}
				if (j<left)
				{
					left = j;
				}
				if (j > right)
				{
					right = j;
				}
			}
		}
	}
	
	Point center;
	aoiGravityCenter(a_timage_in, center);
	Point middle(abs(left - right) * 7.5 / 10 + left, abs(leg - head) / 2 + head);
	Mat a_timage_noWrist = a_timage_in.clone();
	for (int i = 0; i < a_timage_in.rows; i++)//�N��ó����h����for�j��
	{
		for (int j = 0; j < a_timage_in.cols; j++)
		{
			if (j>middle.x)
			{
				a_timage_noWrist.at<uchar>(i, j) = 0;
			}
		}
	}
	Point P = center;
	Point P_ref = middle;
	Mat a_timage_canny = a_timage_in.clone();
	//	Canny(a_timage_noWrist, a_timage_canny, 50, 50, 3);
	//	a_timage_canny = thinning(a_timage_canny);
	Mat a_canny_test;
	Mat a_canny_test_2;
	medianBlur(a_timage_noWrist, a_canny_test, 5);
	Canny(a_canny_test, a_canny_test_2, 50, 50, 3);
	a_canny_test_2 = thinning(a_canny_test_2);
	a_timage_canny = a_canny_test_2.clone();
	int patch_start = 0, patch_end = 0;
	//	imwrite("binary hand.bmp", a_canny_test);

	for (int i = 0; i < a_timage_canny.rows; i += a_timage_canny.rows - 1)
	{
		for (int j = 0; j < a_timage_canny.cols; j++)//���F��x�W�X�Ϥ���t�ɡA���������ϰ���I�A�����N�I�ɤW�h�A�õL�S�O��k�A�æ첾��U��@��pixel�B
		{
			if (a_timage_canny.at<uchar>(i, j) == 255 && a_timage_canny.at<uchar>(i, j - 1) == 0 && a_timage_canny.at<uchar>(i, j + 1) == 0)
			{
				patch_start = j;
				j++;
				while (a_timage_canny.at<uchar>(i, j) == 0)
				{
					j++;
					if (a_timage_canny.at<uchar>(i, j) == 255)
					{
						patch_end = j + 1;
						for (int a = patch_start; a < patch_end; a++)
						{
							if (i == 0)
							{
								a_timage_canny.at<uchar>(i + 1, a) = 255;
								a_timage_canny.at<uchar>(i, a) = 0;
							}
							else if (i == a_timage_canny.rows - 1)
							{
								a_timage_canny.at<uchar>(i - 1, a) = 255;
								a_timage_canny.at<uchar>(i, a) = 0;
							}
						}
						break;
					}
					if (j > a_timage_canny.cols - 2)
					{
						break;
					}
				}
			}
		}
	}
	for (int j = 0; j < a_timage_canny.cols; j += a_timage_canny.cols - 1)
	{
		for (int i = 0; i < a_timage_canny.rows; i++)//���F��x�W�X�Ϥ���t�ɡA���������ϰ���I�A�����N�I�ɤW�h�A�õL�S�O��k�A�æ첾��U��@��pixel�B
		{
			if (a_timage_canny.at<uchar>(i, j) == 255 && a_timage_canny.at<uchar>(i - 1, j) == 0 && a_timage_canny.at<uchar>(i + 1, j) == 0)
			{
				patch_start = i;
				i++;
				while (a_timage_canny.at<uchar>(i, j) == 0)
				{
					i++;
					if (a_timage_canny.at<uchar>(i, j) == 255)
					{
						patch_end = i + 1;
						for (int a = patch_start; a < patch_end; a++)
						{
							if (j == 0)
							{
								a_timage_canny.at<uchar>(a, j + 1) = 255;
								a_timage_canny.at<uchar>(a, j) = 0;
							}
							else if (i == a_timage_canny.cols - 1)
							{
								a_timage_canny.at<uchar>(a, j - 1) = 255;
								a_timage_canny.at<uchar>(a, j) = 0;
							}
						}
						break;
					}
					if (i > a_timage_canny.rows - 2)
					{
						break;
					}
				}
			}
		}
	}
	for (int j = 0; j < a_timage_canny.cols; j++)
	{
		if (a_timage_canny.at<uchar>(0, j) == 255)
		{
			a_timage_canny.at<uchar>(0, j) = 0;
			a_timage_canny.at<uchar>(1, j) = 255;
		}
		if (a_timage_canny.at<uchar>(a_timage_canny.rows - 1, j) == 255)
		{
			a_timage_canny.at<uchar>(a_timage_canny.rows - 1, j) = 0;
			a_timage_canny.at<uchar>(a_timage_canny.rows - 1 - 1, j) = 255;;
		}
	}
	for (int i = 0; i < a_timage_canny.rows; i++)
	{
		if (a_timage_canny.at<uchar>(i, 0) == 255)
		{
			a_timage_canny.at<uchar>(i, 0) = 0;
			a_timage_canny.at<uchar>(i, 1) = 255;
		}
		if (a_timage_canny.at<uchar>(i, a_canny_test.cols - 1) == 255)
		{
			a_timage_canny.at<uchar>(i, a_canny_test.cols - 1) = 0;
			a_timage_canny.at<uchar>(i, a_canny_test.cols - 2) = 255;
		}
	}




	for (int i = 0; i < a_timage_canny.rows; i++)//�NP_ref�V��2�I�k�䪺canny�Ϲ��R��
	{
		for (int j = 0; j < a_timage_canny.cols; j++)
		{
			if (j > P_ref.x - 2)
			{
				a_timage_canny.at<uchar>(i, j) = 0;
			}
		}
	}
	Point detection(P_ref.x, P_ref.y);

	find_start_track_point_main(a_timage_canny, P_ref, detection, head, left);
	for (int i = 0; i < a_timage_canny.rows; i++)//�N�_�l�H�I�k�䪺canny�Ϲ��R��
	{
		for (int j = 0; j < a_timage_canny.cols; j++)
		{
			if (j > detection.x)
			{
				a_timage_canny.at<uchar>(i, j) = 0;
			}
		}
	}
	Mat reg(a_timage_canny.rows, a_timage_canny.cols, a_timage_canny.type(), Scalar(0));
	Mat euc_distance(3, a_timage_canny.rows*a_timage_canny.cols / 10, CV_32F, Scalar(0));//�Ψ��x�s�ڦ��Z������ơA�Ĥ@�q�D�s�ڦ��Z���A�ĤG�q�D�sx�A�ĤT�q�D�sy
	int euc_counter = 0;
	int watchdog = 0;
	int track_check = 1;
	Point detection_old = detection;
	a_timage_canny.at<uchar>(detection) = 128;
	track_plamedge2image_main(track_check, watchdog, detection, detection_old, a_timage_canny, euc_distance, reg, P_ref, euc_counter);

/*	try
	{
		track_plamedge2image_main(track_check, watchdog, detection, detection_old, a_timage_canny, euc_distance, reg, P_ref, euc_counter);
	}
	catch (const char* str)
	{
		cout << "���~:" << str << endl;
		fstream  fp;  //�� �i�@��fstream�� ��
		fp.open("error_log.txt", ios::out);
		if (!fp) //�p�G�}���ɮץ��ѡAfp��0�F���\�Afp���D0     
		{
			cout << "Fail to open error_log " << endl;
		}
	//	fp << name << endl;
		fp.close();
	}*/
	int detection_near_num = 180;//�˴�����X��pxiel�|�X�{�U�@�I
								 /*while (euc_counter<1200)
								 {
								 near_point_find_main(a_timage_canny, detection_near_num, detection, detection_old, watchdog, track_check, euc_distance, euc_counter, reg, P_ref);
								 }*/
	Point P_reg(1, 1);
	float peak_x[100];
	float peak_y[100];
	float valley_x[100];
	float valley_y[100];
	int peak_counter = 0;
	int valley_counter = 0;
	int check_length = 100;//�ˬd�e���`�@�X��pixel;


	find_peak_vally_point_main(check_length, euc_counter, P_reg, euc_distance, peak_x, peak_y, peak_counter, valley_x, valley_y, valley_counter);
	Point P1, P2;
	if (valley_counter <= 3)
	{
		near_point_find_main(a_timage_canny, detection_near_num, detection, detection_old, watchdog, track_check, euc_distance, euc_counter, reg, P_ref);
		find_peak_vally_point_main(check_length, euc_counter, P_reg, euc_distance, peak_x, peak_y, peak_counter, valley_x, valley_y, valley_counter);
	}
	
	for (int i = 0; i < 10; i++)
	{
		for (int j = i + 1; j < 10; j++)
		{
			if (valley_x[i] == valley_x[j] && valley_y[i] == valley_y[j])
			{
				valley_x[j] = NULL;
				valley_y[j] = NULL;
			}
			if (peak_x[i] == peak_x[j] && peak_y[i] == peak_y[j])
			{
				peak_x[j] = NULL;
				peak_y[j] = NULL;
			}
			/*if (detection.x - 10 <= valley_x[i])
			{
			if (valley_x[i] <= detection.x + 10)
			{
			valley_x[j] = NULL;
			valley_y[j] = NULL;
			}
			}
			if (detection.y - 10 <= valley_y[i] <= detection.y + 10)
			{
			if (valley_y[i] <= detection.y + 10)
			{
			valley_x[j] = NULL;
			valley_y[j] = NULL;
			}
			}*/
		}
	}
	for (int i = 0; i < 10; i++)
	{
		for (int j = i + 1; j < 10; j++)
		{
			if (valley_x[j] == NULL&&valley_y[j] == NULL)
			{
				if (j + 1 < 10)
				{
					valley_x[j] = valley_x[j + 1];
					valley_y[j] = valley_y[j + 1];
					valley_x[j + 1] = NULL;
					valley_y[j + 1] = NULL;
				}
			}
		}
	}
	

	string hand_check = hand;
	vector <int> axis_x;
	vector <int> axis_y;
	for (int i = 0; i < 100; i++)
	{
		if (valley_x[i] > 0)
			axis_x.push_back(valley_x[i]);
		else
		{
			break;
		}
	}
//	imwrite("wirst.png", a_timage_in);
//	imwrite("nowirst.png", a_timage_noWrist);
	for (int i = 0; i < 100; i++)
	{
		if (valley_y[i] > 0)
			axis_y.push_back(valley_y[i]);
		else
		{
			break;
		}
	}
	if (axis_x.size() < 4)
	{
		//waitKey(0);
		axis_x.clear();
		vector <int>().swap(axis_x);
		axis_y.clear();
		vector <int>().swap(axis_y);
		throw "�����4�Ө��I";
		return Mat();
	}
	else if (axis_y.size() < 4)
	{
		//waitKey(0);
		axis_x.clear();
		vector <int>().swap(axis_x);
		axis_y.clear();
		vector <int>().swap(axis_y);
		throw "�����4�Ө��I";
		return Mat();
	}
	else
	{
		Mat reg_train(1, 8, CV_32F, Scalar(0));
		for (int i = 0; i < 8; i += 2)
		{
			reg_train.at<float>(0, i) = axis_x.at(i / 2);
			reg_train.at<float>(0, i + 1) = axis_y.at(i / 2);
		}
	//	Ptr<SVM> svm = SVM::load("hand_check.sav");
	//	int check = svm->predict(reg_train);
		
		axis_x.clear();
		vector <int>().swap(axis_x);
		axis_y.clear();
		vector <int>().swap(axis_y);
		/*if (check == 0)
		{
			hand_check = "r";
		}
		else if (check == 1)
			hand_check = "l";
		if (hand_check != hand)
			waitKey(10);*/


		/******************************************�U�������w���|���u���T���I�����p�X�{*************************************/


		if (hand_check == "l")
		{
			P2.x = valley_x[2];
			P2.y = valley_y[2];
			P1.x = valley_x[3];
			P1.y = valley_y[3];
		}
		if (hand_check == "r")
		{
			if (valley_x[3] == 0 && valley_y[3] == 0)
			{
				P1.x = valley_x[0];
				P1.y = valley_y[0];
				P2.x = valley_x[2];
				P2.y = valley_y[2];
			}
			else
			{
				P1.x = valley_x[0];
				P1.y = valley_y[0];
				P2.x = valley_x[1];
				P2.y = valley_y[1];
			}
		}
	}
		Mat a_draw_point(a_canny_test.rows, a_canny_test.cols, CV_8UC3);//= a_canny_test.clone();
		cvtColor(a_canny_test, a_draw_point, CV_GRAY2BGR);
		//	a_draw.convertTo(a_draw, CV_32FC3, 1, 0);
		Point center_1;
		for (int i = 0; i < 5; i++)
		{
			center_1.x = valley_x[i];
			center_1.y = valley_y[i];
			circle(a_draw_point, center_1, 5, Scalar(0, 255, 0), 3, 8, 0);
		}
		for (int i = 0; i < 5; i++)
		{
			center_1.x = peak_x[i];
			center_1.y = peak_y[i];
			circle(a_draw_point, center_1, 5, Scalar(0, 0, 255), 3, 8, 0);
		}
	//	imwrite("peak_vally.bmp", a_draw_point);
		//	hconcat(g_utimage_label, g_utimage_in, g_utimage_label);
		Mat free = a_timage_plam.clone();
		free.at<uchar>(P2) = 255;
		free.at<uchar>(P1) = 255;
		double xp2 = P2.x, xp1 = P1.x, yp2 = P2.y, yp1 = P1.y;
		Mat a_draw_line(a_canny_test.rows, a_canny_test.cols, CV_8UC3);
		cvtColor(a_timage_plam, a_draw_line, CV_GRAY2BGR);
		circle(a_draw_line, P1, 5, Scalar(0, 0, 255), 3, 8, 0);
		circle(a_draw_line, P2, 5, Scalar(0, 0, 255), 3, 8, 0);
		line(a_draw_line, P1, P2, Scalar(0, 255, 0), 2, 8, 0);
		int length = Euclidean_distance(P1, P2);
		if (hand_check == "r")
		{
			center_1.x = P1.x - length;
			center_1.y = P1.y;
			line(a_draw_line, P1, center_1, Scalar(255, 255, 255), 2, 8, 0);
		}
		if (hand_check == "l")
		{
			center_1.x = P1.x - length;
			center_1.y = P1.y;
			line(a_draw_line, P1, center_1, Scalar(255, 255, 255), 2, 8, 0);
		}
		double theta;
		if (hand_check == "l")
		{
			theta = (xp2 - xp1) / (yp2 - yp1);
		}
		else
		{
			theta = (xp1 - xp2) / (yp1 - yp2);
		}
		reg.at<uchar>(P1) = 128;
		reg.at<uchar>(P2) = 128;
		theta = atan(theta) * 180 / M_PI;
		int scale = 1;
		Mat a_rotate_plam;
		Mat rot_mat;
		if (hand_check == "l")
		{
			theta -= 90;
			theta *= -1;
			rot_mat = getRotationMatrix2D(P1, theta, scale);
		}
		else if (hand_check == "r")
		{
			theta += 90;
			theta *= -1;
			rot_mat = getRotationMatrix2D(P1, theta, scale);
		}
		
	Mat a_timage_rot;//�G�Ȥƪ���x����Ϥ�
	//imwrite("draw_point.png", a_draw_point);
	warpAffine(a_timage_plam, a_rotate_plam, rot_mat, a_timage_plam.size());
	warpAffine(a_timage_in, a_timage_rot, rot_mat, a_timage_plam.size());
	//warpAffine(a_draw_point, a_draw_point, rot_mat, a_timage_plam.size());
	int ROI_size = 0;//��l�ƪ�size
	int ROI_size_watchddog = 0;
	int ROI_right;
	int ROI_left;
	if (hand_check == "r")
	{
		P1.y += 10;//�]�i��|���P���٦���L�P����0�I�s�b�A�]���������C
	}
	if (hand_check == "l")
	{
		P1.y -= 10;//�]�i��|���P���٦���L�P����0�I�s�b�A�]���������C
	}
	for (int j = P1.x; j < a_timage_rot.cols; j++)
	{
		if (a_timage_rot.at<uchar>(P1.y, j) == 0)
		{
			ROI_right = j - 1;
			break;
		}
	}
	for (int j = P1.x; j > 0; j--)
	{
		if (a_timage_rot.at<uchar>(P1.y, j) == 0)
		{
			ROI_left = j + 1;
			break;
		}
	}
	
	Mat a_timage_sub_region;
	int sub_region_length = ROI_right-ROI_left;
	watchdog = 0;
	int disp_watchdog = 0;
	int disp_testing_frequency = 5;
	int disp_arry[5] = { 1,5,10,15,20 };//�x�}���ץ����Mdisp_testing_frequency�ۦP
	int disp_counter_arry[5];//�x�}���ץ����Mdisp_testing_frequency�ۦP
	int initialization_display;//�@�}�l�������q
	int counter = 1;
	int disp_counter = 1;
	int left_watchdog = 0, right_watchdog = 0;
	int left_counter = 1, right_counter = 1;
	Scalar image_mean;
	Mat a_timage_out;
//	Mat ROI_SD, ROI_mean;
	strong_img = a_rotate_plam.clone();
//	imwrite("draw_point_totation.png", a_draw_point);
//	imwrite("draw_line.png", a_draw_line);
	int check_width = 80;
	int check_hight = 20;
	Mat check_low_point =  a_timage_rot(Rect(P1.x- check_width/2, P1.y- check_hight/2, check_width, check_hight));//(���Wx,y,w,h)
	Mat check_sum(check_hight, 1, CV_32F, Scalar(0));
	for(int i=0;i<check_low_point.rows;i++)
		for (int j = 0; j < check_low_point.cols; j++)
		{
			check_sum.at<float>(i) += check_low_point.at<uchar>(i, j);
		}
	for (int i = 0; i < check_sum.rows; i++)
		if (check_sum.at<float>(i)==255* check_width)
		{
			if (check_sum.at<float>(i)>check_sum.rows/2)
			{
				P1.y = P1.y + (i - (check_sum.rows / 2));
				break;
			}
		}
	if (hand_check == "r")
	{
		while (left_watchdog == 0)
		{
			a_timage_sub_region = a_timage_rot(Rect(P1.x - left_counter, P1.y, left_counter , left_counter ));//�o�̪�ROI�O���l�ϰ�(sub-region)
			image_mean = mean(a_timage_sub_region);
			//		meanStdDev(a_timage_sub_region, ROI_mean, ROI_SD);
			if (image_mean[0] != 255)
			{
				left_watchdog = 1;
				left_counter--;
				Mat whitetest(left_counter, left_counter, CV_8U, Scalar(255));
				//	a_timage_out = a_rotate_plam(Rect(P1.x - left_counter, P1.y - left_counter, left_counter, left_counter));
				a_timage_sub_region = strong_img(Rect(P1.x - left_counter, P1.y, left_counter, left_counter));
				whitetest.copyTo(a_timage_sub_region);
		//		imwrite("left_Iteration.png", strong_img);
				break;
			}
			left_counter++;

		}
		strong_img = a_rotate_plam.clone();
		while (right_watchdog == 0)
		{
			a_timage_sub_region = a_timage_rot(Rect(P1.x, P1.y , right_counter, right_counter));
			image_mean = mean(a_timage_sub_region);
			//	meanStdDev(a_timage_sub_region, ROI_mean, ROI_SD);
			if (image_mean[0] != 255)
			{
				right_watchdog = 1;
				right_counter--;
				Mat whitetest(right_counter, right_counter, CV_8U, Scalar(255));
				//	a_timage_out = a_rotate_plam(Rect(P1.x, P1.y - right_counter, right_counter, right_counter));
				a_timage_sub_region = strong_img(Rect(P1.x, P1.y , right_counter, right_counter));
				whitetest.copyTo(a_timage_sub_region);
		//		imwrite("right_Iteration.png", strong_img);
				break;
			}
			right_counter++;
		}
		counter = left_counter + right_counter;
		strong_img = a_rotate_plam.clone();
		while (watchdog == 0)
		{
			a_timage_sub_region = a_timage_rot(Rect(P1.x - left_counter, P1.y , counter, counter));
			image_mean = mean(a_timage_sub_region);
			if (image_mean[0] == 255)
			{
				watchdog = 1;
				counter--;
				Mat whitetest(counter, counter, CV_8U, Scalar(255));
				a_timage_out = a_rotate_plam(Rect(P1.x - left_counter, P1.y, counter, counter));
				a_timage_sub_region = strong_img(Rect(P1.x - left_counter, P1.y, counter, counter));
				whitetest.copyTo(a_timage_sub_region);
				break;
			}
			counter--;
		}
		strong_img = a_rotate_plam.clone();
		initialization_display = disp_arry[disp_watchdog];
		while (disp_watchdog != disp_testing_frequency)//X��V�����˴�
		{
			a_timage_sub_region = a_timage_rot(Rect(P1.x - left_counter + disp_counter + initialization_display , P1.y, counter+disp_counter, counter+disp_counter));
			image_mean = mean(a_timage_sub_region);
			if (image_mean[0] != 255)
			{
				disp_watchdog++;
				disp_counter--;
				disp_counter_arry[disp_watchdog - 1] = disp_counter;
				a_timage_sub_region = strong_img(Rect(P1.x - left_counter + disp_counter + initialization_display, P1.y, counter + disp_counter, counter + disp_counter));
				if (disp_watchdog == disp_testing_frequency)
				{
					int max_counter = 0;
					for (int x = 0; x < disp_testing_frequency; x++)
					{
						if (disp_counter_arry[x] > max_counter)
						{
							max_counter = disp_counter_arry[x];
							initialization_display = x;
						}
					}
					if (initialization_display > (disp_testing_frequency - 1))
						initialization_display = 0;
					else
						initialization_display = disp_arry[initialization_display];
					disp_counter = max_counter;
					a_timage_sub_region = strong_img(Rect(P1.x - left_counter + disp_counter + initialization_display, P1.y, counter + disp_counter, counter + disp_counter));
					if (a_timage_out.rows < a_timage_sub_region.rows)
					{
						Mat whitetest(counter + disp_counter, counter + disp_counter, CV_8U, Scalar(255));
						a_timage_out = a_rotate_plam(Rect(P1.x - left_counter + disp_counter + initialization_display, P1.y, counter + disp_counter, counter + disp_counter));
						a_timage_sub_region = strong_img(Rect(P1.x - left_counter + disp_counter + initialization_display, P1.y, counter + disp_counter, counter + disp_counter));
						whitetest.copyTo(a_timage_sub_region);
						
						break;
					}
					else
					{
						Mat whitetest(counter, counter, CV_8U, Scalar(255));
						a_timage_out = a_rotate_plam(Rect(P1.x - left_counter, P1.y, counter, counter));
						a_timage_sub_region = strong_img(Rect(P1.x - left_counter, P1.y, counter, counter));
						whitetest.copyTo(a_timage_sub_region);
					}
					break;
				}
				initialization_display = disp_arry[disp_watchdog];
				disp_counter = 0;//�]�U���|+1�A�]���o�̪�l�Ƭ�0
			}
			disp_counter++;
		}
	}
	if (hand_check == "l")
	{
		while (left_watchdog==0)
		{
			a_timage_sub_region = a_timage_rot(Rect(P1.x - left_counter , P1.y - left_counter , left_counter , left_counter ));//�o�̪�ROI�O���l�ϰ�(sub-region)
			image_mean = mean(a_timage_sub_region);
	//		meanStdDev(a_timage_sub_region, ROI_mean, ROI_SD);
			if (image_mean[0] != 255)
			{
				left_watchdog = 1;
				left_counter--;
				Mat whitetest(left_counter, left_counter, CV_8U, Scalar(255));
			//	a_timage_out = a_rotate_plam(Rect(P1.x - left_counter, P1.y - left_counter, left_counter, left_counter));
				a_timage_sub_region = strong_img(Rect(P1.x - left_counter, P1.y - left_counter, left_counter, left_counter));
				whitetest.copyTo(a_timage_sub_region);
		//		imwrite("left_Iteration.png", strong_img);
				break;
			}
			left_counter++; 
			
		}
		strong_img = a_rotate_plam.clone();
		while (right_watchdog==0)
		{
			a_timage_sub_region = a_timage_rot(Rect(P1.x, P1.y - right_counter, right_counter, right_counter));
			image_mean = mean(a_timage_sub_region);
		//	meanStdDev(a_timage_sub_region, ROI_mean, ROI_SD);
			if (image_mean[0] != 255)
			{
				right_watchdog = 1;
				right_counter--;
				Mat whitetest(right_counter, right_counter, CV_8U, Scalar(255));
			//	a_timage_out = a_rotate_plam(Rect(P1.x, P1.y - right_counter, right_counter, right_counter));
				a_timage_sub_region = strong_img(Rect(P1.x, P1.y - right_counter, right_counter, right_counter));
				whitetest.copyTo(a_timage_sub_region);
			//	imwrite("right_Iteration.png", strong_img);
				break;
			}
			right_counter++;
		}
		counter = left_counter + right_counter;
		strong_img = a_rotate_plam.clone();
		while (watchdog==0)
		{
			if (P1.y<counter)
			{
				counter = P1.y;
			}
			a_timage_sub_region = a_timage_rot(Rect(P1.x - left_counter, P1.y - counter, counter, counter));
			image_mean = mean(a_timage_sub_region);
			if (image_mean[0] == 255)
			{
				watchdog = 1;
				counter--;
				Mat whitetest(counter, counter, CV_8U, Scalar(255));
				a_timage_out = a_rotate_plam(Rect(P1.x - left_counter, P1.y - counter, counter, counter));
				a_timage_sub_region = strong_img(Rect(P1.x - left_counter, P1.y - counter, counter, counter));
				whitetest.copyTo(a_timage_sub_region);
				break;
			}
			counter--;
		}
		strong_img = a_rotate_plam.clone();
		initialization_display = disp_arry[disp_watchdog];
		while (disp_watchdog!= disp_testing_frequency)//X��V�����˴�
		{
			a_timage_sub_region = a_timage_rot(Rect(P1.x - left_counter + disp_counter+ initialization_display, P1.y - counter-disp_counter, counter+ disp_counter, counter+ disp_counter));
			image_mean = mean(a_timage_sub_region);
			if (image_mean[0] != 255)
			{
				disp_watchdog++;
				disp_counter--;
				disp_counter_arry[disp_watchdog - 1] = disp_counter;
				a_timage_sub_region = strong_img(Rect(P1.x - left_counter + disp_counter, P1.y - counter - disp_counter, counter + disp_counter, counter + disp_counter));
				if (disp_watchdog == disp_testing_frequency)
				{
					int max_counter = 0;
					for (int x = 0; x < disp_testing_frequency; x++)
					{
						if (disp_counter_arry[x] > max_counter)
						{
							max_counter = disp_counter_arry[x];
							initialization_display = x;
						}
					}
					if (initialization_display > (disp_testing_frequency - 1))
						initialization_display = 0;
					else
						initialization_display = disp_arry[initialization_display];
					disp_counter = max_counter;
					a_timage_sub_region = strong_img(Rect(P1.x - left_counter + disp_counter, P1.y - counter - disp_counter, counter + disp_counter, counter + disp_counter));
					if (a_timage_out.rows < a_timage_sub_region.rows)
					{
						Mat whitetest(counter + disp_counter, counter + disp_counter, CV_8U, Scalar(255));
						a_timage_out = a_rotate_plam(Rect(P1.x - left_counter + disp_counter + initialization_display, P1.y - counter - disp_counter, counter + disp_counter, counter + disp_counter));
						a_timage_sub_region = strong_img(Rect(P1.x - left_counter + disp_counter + initialization_display, P1.y - counter - disp_counter, counter + disp_counter, counter + disp_counter));
						whitetest.copyTo(a_timage_sub_region);
						break;
					}
					else
					{
						Mat whitetest(counter, counter, CV_8U, Scalar(255));
						a_timage_out = a_rotate_plam(Rect(P1.x - left_counter, P1.y - counter, counter, counter));
						a_timage_sub_region = strong_img(Rect(P1.x - left_counter, P1.y - counter, counter, counter));
						whitetest.copyTo(a_timage_sub_region);
					}
					break;
				}
				initialization_display = disp_arry[disp_watchdog];
				disp_counter = 0;//�]�U���|+1�A�]���o�̪�l�Ƭ�0
			}
			disp_counter++;
		}
	}
//	imwrite("total_Iteration.png", strong_img);
	a_timage_in = strong_img.clone();
	hconcat(a_draw_line, a_draw_point, strong_img);
	return a_timage_out;
}



Mat labeling_replace(Mat a_timage_in)
{
	Mat a_timage_out = a_timage_in.clone();
	int label = L_BASE;
	int number_counter = 1;
	int total_pixel = a_timage_in.rows*a_timage_in.cols;
	while (number_counter != 0)
	{
		for (int i = 0; i < a_timage_in.rows; i++)
		{
			for (int j = 0; j < a_timage_in.cols; j++)
			{
				if (a_timage_in.at<uchar>(i, j) == label)	number_counter++;
			}
		}
		if (number_counter < total_pixel / 100 * 3)
		{
			for (int i = 0; i < a_timage_in.rows; i++)
			{
				for (int j = 0; j < a_timage_in.cols; j++)
				{
					if (a_timage_in.at<uchar>(i, j) == label)
					{
						a_timage_out.at<uchar>(i, j) = 0;
					}
				}
			}
		}
		else if (number_counter>total_pixel / 100 * 3)
		{
			for (int i = 0; i < a_timage_in.rows; i++)
			{
				for (int j = 0; j < a_timage_in.cols; j++)
				{
					if (a_timage_in.at<uchar>(i, j) == label)
					{
						a_timage_out.at<uchar>(i, j) = 255;
					}
				}
			}
		}
		label++;
		if (number_counter != 1)
		{
			number_counter = 1;
		}
		else if (number_counter == 1)
		{
			number_counter = 0;
		}
	}
	return a_timage_out;
}
/*--- labelset --- ��۳s�����������[���Ҹ��X ------------------------------------------------------------------
image: �v���}�C
xs, ys: �_�l��m
label: ���ҧǸ�
----------------------------------------------------------------------------------------------------------------*/
void labelset(Mat a_timage_in, int xs, int ys, int label)
{
	int	i, j, cnt, im, ip, jm, jp;

	a_timage_in.at<uchar>(ys, xs) = label;
	for (;;) {
		cnt = 0;
		for (i = 0; i < a_timage_in.rows; i++)
			for (j = 0; j < a_timage_in.cols; j++)
				if (a_timage_in.at<uchar>(i, j) == label) {
					im = i - 1; ip = i + 1; jm = j - 1; jp = j + 1;
					if (im < 0) im = 0; if (ip >= a_timage_in.rows) ip = a_timage_in.rows - 1;
					if (jm < 0) jm = 0; if (jp >= a_timage_in.cols) jp = a_timage_in.cols - 1;
					if (a_timage_in.at<uchar>(i, jp) == HIGH) {
						a_timage_in.at<uchar>(i, jp) = label; cnt++;
					}
					if (a_timage_in.at<uchar>(im, jp) == HIGH) {
						a_timage_in.at<uchar>(im, jp) = label; cnt++;
					}
					if (a_timage_in.at<uchar>(im, j) == HIGH) {
						a_timage_in.at<uchar>(im, j) = label; cnt++;
					}
					if (a_timage_in.at<uchar>(im, jm) == HIGH) {
						a_timage_in.at<uchar>(im, jm) = label; cnt++;
					}
					if (a_timage_in.at<uchar>(i, jm) == HIGH) {
						a_timage_in.at<uchar>(i, jm) = label; cnt++;
					}
					if (a_timage_in.at<uchar>(ip, jm) == HIGH) {
						a_timage_in.at<uchar>(ip, jm) = label; cnt++;
					}
					if (a_timage_in.at<uchar>(ip, j) == HIGH) {
						a_timage_in.at<uchar>(ip, j) = label; cnt++;
					}
					if (a_timage_in.at<uchar>(ip, jp) == HIGH) {
						a_timage_in.at<uchar>(ip, jp) = label; cnt++;
					}
				}
		if (cnt == 0) break;
	}
}
int labeling(Mat a_timage_in, Mat &imagelabel)
{
	int	i, j, label;

	for (i = 0; i < a_timage_in.rows; i++)
		for (j = 0; j < a_timage_in.cols; j++)
			a_timage_in.at<uchar>(i, j) = 255 - a_timage_in.at<uchar>(i, j);
	imagelabel = a_timage_in.clone();
	label = L_BASE;
	for (i = 0; i < a_timage_in.rows; i++)
		for (j = 0; j < a_timage_in.cols; j++) {
			if (imagelabel.at<uchar>(i, j) == HIGH) {
				if (label >= HIGH) return -1;
				labelset(imagelabel, j, i, label); label++;
			}
		}
	return label - L_BASE;
}