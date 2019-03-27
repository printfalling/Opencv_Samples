#pragma once
#include <opencv2/core.hpp>
#include <opencv2/core/utility.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <iostream>
#include <sstream>

#define w 400

using namespace std;
using namespace cv;

namespace Tutorial {



	class Tutorial
	{
	public:
		Tutorial();
		~Tutorial();
		int how_to_scan_images(const int argc, const char* argv[]);
		int mask_operation(const int argc, const char* argv[]);
		int add_image(const int argc, const char* argv[]);
		int change_contract(const int argc, const char* argv[]);
		int discrete_fourier_transformation(const int argc, const char* argv[]);
		int xml_yaml_in_out(const int argc, const char* argv[]);
		int basic_draw(const int argc, const char* argv[]);
		int draw_random_lines(const int argc, const char* argv[]);

	private:
		Mat& ScanImageAndReduceC(Mat& I, const uchar* const table);
		Mat& ScanImageAndReduceIterator(Mat& I, const uchar* const table);
		Mat& ScanImageAndReduceRandomAccess(Mat& I, const uchar* const table);

		void Sharpen(const Mat& myImg, Mat& Result);

		// basic draw
		void MyEllipse(Mat img, double angle);
		void MyFilledCircle(Mat img, Point center);
		void MyPolygon(Mat img);
		void MyLine(Mat img, Point start, Point end);

	};

	class MyData
	{
	public:
		MyData() : A(0), X(0), id()
		{}
		explicit MyData(int _value) : A(_value), X(CV_PI), id("mydata123") // explicit to avoid implicit conversion
		{}
		void write(FileStorage& fs) const                        //Write serialization for this class
		{
			fs << "{" << "A" << A << "X" << X << "id" << id << "}";
		}
		void read(const FileNode& node)                          //Read serialization for this class
		{
			A = (int)node["A"];
			X = (double)node["X"];
			id = (string)node["id"];
		}
	public:   // Data Members
		int A;
		double X;
		string id;
	};

}

