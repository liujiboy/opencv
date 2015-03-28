#include<opencv2/opencv.hpp>
#include<opencv2/features2d/features2d.hpp>
#include<opencv2/nonfree/nonfree.hpp>
#include<vector>
#include<iostream>
using namespace std;
using namespace cv;
//直接比较特征描述
void bruteForce(const Mat&im1,const Mat&im2,const vector<KeyPoint>&kp1,const vector<KeyPoint>&kp2,const Mat&des1,const Mat&des2,Mat&outImg)
{
	//true表示进行双向比较
	BFMatcher matcher(NORM_L2,true);
	vector<DMatch> matches;
	matcher.match(des1,des2,matches);
	drawMatches(im1,kp1,im2,kp2,matches,outImg);
}
//用Ratio Test过滤错误的特征描述
void ratioTest(const Mat&im1,const Mat&im2,const vector<KeyPoint>&kp1,const vector<KeyPoint>&kp2,const Mat&des1,const Mat&des2,Mat&outImg)
{
	BFMatcher matcher;
	vector<vector<DMatch> > matchesVector;
	//找到前两个最近的匹配	
	matcher.knnMatch(des1,des2,matchesVector,2);
	vector<DMatch>	matches;
	for(vector<vector<DMatch> >::iterator it=matchesVector.begin(); it!=matchesVector.end();it++)
	{
		DMatch m1=(*it)[0];
		DMatch m2=(*it)[1];
		//两个最近匹配的distance比值小于0.8
		if(m1.distance/m2.distance<0.8)
			matches.push_back(m1);
	}
	drawMatches(im1,kp1,im2,kp2,matches,outImg);
}
//用基础矩阵加Ransac
void ransac(const Mat&im1,const Mat&im2,const vector<KeyPoint>&kp1,const vector<KeyPoint>&kp2,const Mat&des1,const Mat&des2,Mat&outImg)
{
	BFMatcher matcher(NORM_L2,true);
	vector<DMatch> matches;
	matcher.match(des1,des2,matches);
	vector<Point2f> points1,points2;
	for(vector<DMatch>::iterator it=matches.begin();it!=matches.end();it++)
	{
		points1.push_back(kp1[it->queryIdx].pt);
		points2.push_back(kp2[it->trainIdx].pt);

	}
	vector<uchar> mask;
	Mat fundamentalMat=cv::findFundamentalMat(points1,points2,FM_RANSAC,3,0.99,mask);
	vector<DMatch> trueMatches;
	for(vector<DMatch>::size_type i=0;i<matches.size();i++)
	{
		if(mask[i])
			trueMatches.push_back(matches[i]);
	}
	drawMatches(im1,kp1,im2,kp2,trueMatches,outImg);
	//cout<<trueMatches.size()<<endl;
//	cout<<fundamentalMat<<endl;
}
//Ratio Test+Ransac
void ratioTestAndRansac(const Mat&im1,const Mat&im2,const vector<KeyPoint>&kp1,const vector<KeyPoint>&kp2,const Mat&des1,const Mat&des2,Mat&outImg)
{
	//首先进行Ratio Test
	BFMatcher matcher;
	vector<vector<DMatch> > matchesVector;
	//找到前两个最近的匹配	
	matcher.knnMatch(des1,des2,matchesVector,2);
	vector<DMatch>	matches;
	for(vector<vector<DMatch> >::iterator it=matchesVector.begin(); it!=matchesVector.end();it++)
	{
		DMatch m1=(*it)[0];
		DMatch m2=(*it)[1];
		//这里取0.8可以包含更多的匹配点
		if(m1.distance/m2.distance<0.8)
			matches.push_back(m1);
	}
	//再进行ransac
	vector<Point2f> points1,points2;
	for(vector<DMatch>::iterator it=matches.begin();it!=matches.end();it++)
	{
		points1.push_back(kp1[it->queryIdx].pt);
		points2.push_back(kp2[it->trainIdx].pt);

	}
	vector<uchar> mask;
	Mat fundamentalMat=cv::findFundamentalMat(points1,points2,FM_RANSAC,3,0.99,mask);
	vector<DMatch> trueMatches;
	for(vector<DMatch>::size_type i=0;i<matches.size();i++)
	{
		if(mask[i])
			trueMatches.push_back(matches[i]);
	}
	drawMatches(im1,kp1,im2,kp2,trueMatches,outImg);


}
int main(int argc,char* argv[])
{

	if(argc!=3)
	{
		cout<<"使用方法：程序名 [图像1路径] [图像2路径]"<<endl;
		return 0;
	}
	Mat im1=imread(argv[1]);
	Mat im2=imread(argv[2]);
	SIFT sift;
	vector<KeyPoint> kp1,kp2;
	Mat des1,des2;
	sift(im1,Mat(),kp1,des1);
	sift(im2,Mat(),kp2,des2);
	Mat outImg1,outImg2,outImg3,outImg4;
	bruteForce(im1,im2,kp1,kp2,des1,des2,outImg1);
	ratioTest(im1,im2,kp1,kp2,des1,des2,outImg2);
	ransac(im1,im2,kp1,kp2,des1,des2,outImg3);
	ratioTestAndRansac(im1,im2,kp1,kp2,des1,des2,outImg4);
	namedWindow("BruteForce",WINDOW_NORMAL);
	imshow("BruteForce",outImg1);
	namedWindow("Ratio Test",WINDOW_NORMAL);
	imshow("Ratio Test",outImg2);
	namedWindow("Ransac",WINDOW_NORMAL);
	imshow("Ransac",outImg3);
	namedWindow("Ratio Test And Ransac",WINDOW_NORMAL);
	imshow("Ratio Test And Ransac",outImg4);
	while(waitKey()!=27)
		;
	destroyAllWindows();
	return 0;
}
