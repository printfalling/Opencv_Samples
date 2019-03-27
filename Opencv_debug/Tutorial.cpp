#include "Tutorial.h"

namespace Tutorial {



	Tutorial::Tutorial()
	{
	}


	Tutorial::~Tutorial()
	{
	}



	static void how_to_scan_images_help()
	{
		cout
			<< "\n--------------------------------------------------------------------------" << endl
			<< "This program shows how to scan image objects in OpenCV (cv::Mat). As use case"
			<< " we take an input image and divide the native color palette (255) with the " << endl
			<< "input. Shows C operator[] method, iterators and at function for on-the-fly item address calculation." << endl
			<< "Usage:" << endl
			<< "./how_to_scan_images <imageNameToUse> <divideWith> [G]" << endl
			<< "if you add a G parameter the image is processed in gray scale" << endl
			<< "--------------------------------------------------------------------------" << endl
			<< endl;
	}

	int Tutorial::how_to_scan_images(const int argc, const char* argv[]) {

		how_to_scan_images_help();
		if (argc < 3)
		{
			cout << "Not enough parameters" << endl;
			return -1;
		}

		Mat I, J;
		if (argc == 4 && !strcmp(argv[3], "G"))
			I = imread(argv[1], IMREAD_GRAYSCALE);
		else
			I = imread(argv[1], IMREAD_COLOR);

		if (I.empty())
		{
			cout << "The image" << argv[1] << " could not be loaded." << endl;
			return -1;
		}

		//! [dividewith]
		int divideWith = 0; // convert our input string to number - C++ style
		stringstream s;
		s << argv[2];
		s >> divideWith;
		if (!s || !divideWith)
		{
			cout << "Invalid number entered for dividing. " << endl;
			return -1;
		}

		uchar table[256];
		for (int i = 0; i < 256; ++i)
			table[i] = (uchar)(divideWith * (i / divideWith));
		//! [dividewith]

		const int times = 10;
		double t;

		t = (double)getTickCount();

		for (int i = 0; i < times; ++i)
		{
			cv::Mat clone_i = I.clone();
			J = ScanImageAndReduceC(clone_i, table);
		}

		t = 1000 * ((double)getTickCount() - t) / getTickFrequency();
		t /= times;

		cout << "Time of reducing with the C operator [] (averaged for "
			<< times << " runs): " << t << " milliseconds." << endl;

		t = (double)getTickCount();

		for (int i = 0; i < times; ++i)
		{
			//cout << i << endl;
			cv::Mat clone_i = I.clone();
			J = ScanImageAndReduceIterator(clone_i, table);
		}

		t = 1000 * ((double)getTickCount() - t) / getTickFrequency();
		t /= times;

		cout << "Time of reducing with the iterator (averaged for "
			<< times << " runs): " << t << " milliseconds." << endl;

		t = (double)getTickCount();

		for (int i = 0; i < times; ++i)
		{
			cv::Mat clone_i = I.clone();
			ScanImageAndReduceRandomAccess(clone_i, table);
		}

		t = 1000 * ((double)getTickCount() - t) / getTickFrequency();
		t /= times;

		cout << "Time of reducing with the on-the-fly address generation - at function (averaged for "
			<< times << " runs): " << t << " milliseconds." << endl;

		//! [table-init]
		Mat lookUpTable(1, 256, CV_8U);
		uchar* p = lookUpTable.ptr();
		for (int i = 0; i < 256; ++i)
			p[i] = table[i];
		//! [table-init]

		t = (double)getTickCount();

		for (int i = 0; i < times; ++i)
			//! [table-use]
			LUT(I, lookUpTable, J);
		//! [table-use]

		t = 1000 * ((double)getTickCount() - t) / getTickFrequency();
		t /= times;

		cout << "Time of reducing with the LUT function (averaged for "
			<< times << " runs): " << t << " milliseconds." << endl;
		getchar();
		return 0;

	}

	static void mask_operate_help(const char* progName)
	{
		cout << endl
			<< "This program shows how to filter images with mask: the write it yourself and the"
			<< "filter2d way. " << endl
			<< "Usage:" << endl
			<< progName << " [image_path -- default ../data/lena.jpg] [G -- grayscale] " << endl << endl;
	}

	int Tutorial::mask_operation(const int argc, const char* argv[]) {
		mask_operate_help(argv[0]);
		if (argc < 2) return -1;
		const char* filename = argv[1];

		Mat src, dst0, dst1;
		if (argc >= 3 && strcmp("G", argv[2]))
			src = imread(filename, IMREAD_GRAYSCALE);
		else
			src = imread(filename, IMREAD_COLOR);

		if (src.empty())
		{
			cerr << "Can't open image [" << filename << "]" << endl;
			return -1;
		}

		//namedWindow("Input", WINDOW_AUTOSIZE);
		//namedWindow("Output", WINDOW_AUTOSIZE);

		//imshow("Input", src);
		double t = (double)getTickCount();

		Sharpen(src, dst0);

		t = ((double)getTickCount() - t) / getTickFrequency();
		cout << "Hand written function time passed in seconds: " << t << endl;

		//imshow("Output", dst0);
		waitKey();

		Mat kernel = (Mat_<uchar>(3, 3) << 0, -1, 0, -1, 5, -1, 0, -1, 0);
		t = (double)getTickCount();

		filter2D(src, dst1, src.depth(), kernel);
		t = ((double)getTickCount() - t) / getTickFrequency();
		cout << "Built-in filter2D time passed in seconds:     " << t << endl;

		//imshow("Output", dst1);

		waitKey();
		getchar();
		return 0;

	}

	int Tutorial::add_image(const int argc, const char* argv[]) {
		double a = 0.5;
		double b;
		double input;

		Mat src1, src2, dst;

		cout << " Simple Linear Blender " << endl;
		cout << "-----------------------" << endl;
		cout << "* Enter alpha [0.0-1.0]: ";
		cin >> input;

		if (input > 0 && input < 1) a = input;
		if (argc < 3) return -1;

		src1 = imread(argv[1]);
		src2 = imread(argv[2]);

		if (src1.empty()) {
			cout << "Error loading src1" << endl;
			getchar();
			return -1;
		}
		if (src2.empty()) {
			cout << "Error loading src2" << endl;
			getchar();
			return -1;
		}

		b = (1 - a);
		addWeighted(src1, a, src2, b, 0.0, dst);
		imshow("Linear Blend", dst);
		waitKey(0);

		return 0;
	}

	int Tutorial::change_contract(const int argc, const char * argv[])
	{
		CommandLineParser parser(argc, argv, "{@input | linux_logo.jpg | input image}");
		Mat image = imread(parser.get<String>("@input"));
		if (image.empty()) {
			cout << "Could not open or find the image!\n" << endl;
			cout << "Usage: " << argv[0] << " <Input image>" << endl;
			return -1;
		}
		Mat new_image = Mat::zeros(image.size(), image.type());
		double alph = 1;
		int bet = 0;
		cout << " Basic Linear Transforms " << endl;
		cout << "-------------------------" << endl;
		cout << "* Enter the alpha value [1.0-3.0]: "; cin >> alph;
		cout << "* Enter the beta value [0-100]: ";    cin >> bet;

		for (int x = 0; x < image.rows; x++)
			for (int y = 0; y < image.cols; y++)
				for (int c = 0; c < image.channels(); c++)
					new_image.at<Vec3b>(y, x)[c] = saturate_cast<uchar>(alph*image.at<Vec3b>(y, x)[c] + bet);
		imshow("Original Image", image);
		imshow("New Image", new_image);
		waitKey();
		return 0;
	}
	//These write and read functions must be defined for the serialization in FileStorage to work
	static void write(FileStorage& fs, const std::string&, const MyData& x)
	{
		x.write(fs);
	}

	static void read(const FileNode& node, MyData& x, const MyData& default_value = MyData()) {
		if (node.empty())
			x = default_value;
		else
			x.read(node);
	}

	// This function will print our custom class to the console
	static ostream& operator<<(ostream& out, const MyData& m)
	{
		out << "{ id = " << m.id << ", ";
		out << "X = " << m.X << ", ";
		out << "A = " << m.A << "}";
		return out;
	}

	static void xml_yaml_in_out_help(const char** av)
	{
		cout << endl
			<< av[0] << " shows the usage of the OpenCV serialization functionality." << endl
			<< "usage: " << endl
			<< av[0] << " outputfile.yml.gz" << endl
			<< "The output file may be either XML (xml) or YAML (yml/yaml). You can even compress it by "
			<< "specifying this in its extension like xml.gz yaml.gz etc... " << endl
			<< "With FileStorage you can serialize objects in OpenCV by using the << and >> operators" << endl
			<< "For example: - create a class and have it serialized" << endl
			<< "             - use it to read and write matrices." << endl;
	}

	int Tutorial::xml_yaml_in_out(const int argc, const char* argv[]) {
		if (argc < 2) {
			xml_yaml_in_out_help(argv);
			return -1;
		}
		string filename = argv[1];
		{	// write
			Mat R = Mat_<uchar>::eye(3, 3),
				T = Mat_<double>::zeros(3, 1);
			MyData m(1);
			FileStorage fs(filename, FileStorage::WRITE);
			fs << "iterationNr" << 100;
			fs << "strings" << "[";                              // text - string sequence
			fs << "image1.jpg" << "Awesomeness" << "windows_logo.jpg";
			fs << "]";                                           // close sequence

			fs << "Mapping";                              // text - mapping
			fs << "{" << "One" << 1;
			fs << "Two" << 2 << "}";

			fs << "R" << R;                                      // cv::Mat
			fs << "T" << T;

			fs << "MyData" << m;                                // your own data structures

			fs.release();                                       // explicit close
			cout << "Write Done." << endl;
		}
		{//read
			cout << endl << "Reading: " << endl;
			FileStorage fs;
			fs.open(filename, FileStorage::READ);

			int itNr;
			//fs["iterationNr"] >> itNr;
			itNr = (int)fs["iterationNr"];
			cout << itNr;
			if (!fs.isOpened())
			{
				cerr << "Failed to open " << filename << endl;
				xml_yaml_in_out_help(argv);
				return 1;
			}

			FileNode n = fs["strings"];                         // Read string sequence - Get node
			if (n.type() != FileNode::SEQ)
			{
				cerr << "strings is not a sequence! FAIL" << endl;
				return 1;
			}

			FileNodeIterator it = n.begin(), it_end = n.end(); // Go through the node
			for (; it != it_end; ++it)
				cout << (string)*it << endl;


			n = fs["Mapping"];                                // Read mappings from a sequence
			cout << "Two  " << (int)(n["Two"]) << "; ";
			cout << "One  " << (int)(n["One"]) << endl << endl;


			MyData m;
			Mat R, T;

			fs["R"] >> R;                                      // Read cv::Mat
			fs["T"] >> T;
			fs["MyData"] >> m;                                 // Read your own structure_

			cout << endl
				<< "R = " << R << endl;
			cout << "T = " << T << endl << endl;
			cout << "MyData = " << endl << m << endl << endl;

			//Show default behavior for non existing nodes
			cout << "Attempt to read NonExisting (should initialize the data structure with its default).";
			fs["NonExisting"] >> m;
			cout << endl << "NonExisting = " << endl << m << endl;
		}
		cout << endl
			<< "Tip: Open up " << filename << " with a text editor to see the serialized data." << endl;
		getchar();
		return 0;
	}


	static void discrete_fourier_transformation_help(void)
	{
		cout << endl
			<< "This program demonstrated the use of the discrete Fourier transform (DFT). " << endl
			<< "The dft of an image is taken and it's power spectrum is displayed." << endl
			<< "Usage:" << endl
			<< "./discrete_fourier_transform [image_name -- default ../data/lena.jpg]" << endl;
	}

	int Tutorial::discrete_fourier_transformation(const int argc, const char * argv[])
	{
		discrete_fourier_transformation_help();
		const char* filename = argc >= 2 ? argv[1] : "windows_logo.jpg";

		Mat I = imread(filename, IMREAD_GRAYSCALE);
		if (I.empty()) {
			cout << "Error opening image" << endl;
			return -1;
		}

		Mat padded;
		int m = getOptimalDFTSize(I.rows);
		int n = getOptimalDFTSize(I.cols);
		copyMakeBorder(I, padded, 0, m - I.rows, 0, n - I.cols, BORDER_CONSTANT, Scalar_<int>::all(0));

		Mat planes[] = { Mat_<float>(padded), Mat::zeros(padded.size(), CV_32F) };
		Mat complexI;
		merge(planes, 2, complexI);
		dft(complexI, complexI);
		// compute the magnitude and switch to logarithmic scale
		// => log(1 + sqrt(Re(DFT(I))^2 + Im(DFT(I))^2))
		split(complexI, planes);
		magnitude(planes[0], planes[1], planes[0]);
		Mat magI = planes[0];

		magI += Scalar::all(1);  // switch to logarithmic scale
		log(magI, magI);
		// crop the spectrum, if it has an odd number of rows or columns
		magI = magI(Rect(0, 0, magI.cols & -2, magI.rows & -2));

		// rearrange the quadrants of Fourier image  so that the origin is at the image center
		int cx = magI.cols / 2;
		int cy = magI.rows / 2;

		Mat q0(magI, Rect(0, 0, cx, cy));   // Top-Left - Create a ROI per quadrant
		Mat q1(magI, Rect(cx, 0, cx, cy));  // Top-Right
		Mat q2(magI, Rect(0, cy, cx, cy));  // Bottom-Left
		Mat q3(magI, Rect(cx, cy, cx, cy)); // Bottom-Right

		Mat tmp;                           // swap quadrants (Top-Left with Bottom-Right)
		q0.copyTo(tmp);
		q3.copyTo(q0);
		tmp.copyTo(q3);

		q1.copyTo(tmp);                    // swap quadrant (Top-Right with Bottom-Left)
		q2.copyTo(q1);
		tmp.copyTo(q2);

		normalize(magI, magI, 0, 1, NORM_MINMAX); // Transform the matrix with float values into a
												// viewable image form (float between values 0 and 1).

		imshow("Input Image", I);    // Show the result
		imshow("spectrum magnitude", magI);
		waitKey();

		return 0;
	}


	int Tutorial::basic_draw(const int argc, const char * argv[])
	{
		char atom_window[] = "Drawing 1: Atom";
		char rook_window[] = "Drawing 2: Rook";

		Mat atom_image = Mat::zeros(w, w, CV_8UC3);
		Mat rook_image = Mat::zeros(w, w, CV_8UC3);


		MyEllipse(atom_image, 90);
		MyEllipse(atom_image, 0);
		MyEllipse(atom_image, 45);
		MyEllipse(atom_image, -45);

		MyFilledCircle(atom_image, Point(w / 2, w / 2));


		MyPolygon(rook_image);

		rectangle(rook_image,
			Point(0, 7 * w / 8),
			Point(w, w),
			Scalar(0, 255, 255),
			FILLED,
			LINE_8);

		MyLine(rook_image, Point(0, 15 * w / 16), Point(w, 15 * w / 16));
		MyLine(rook_image, Point(w / 4, 7 * w / 8), Point(w / 4, w));
		MyLine(rook_image, Point(w / 2, 7 * w / 8), Point(w / 2, w));
		MyLine(rook_image, Point(3 * w / 4, 7 * w / 8), Point(3 * w / 4, w));

		imshow(atom_window, atom_image);
		moveWindow(atom_window, 0, 200);
		imshow(rook_window, rook_image);
		moveWindow(rook_window, w, 200);

		waitKey(0);
		return(0);
	}

	Mat& Tutorial::ScanImageAndReduceC(Mat& I, const uchar* table) {
		CV_Assert(I.depth() == CV_8U);

		const int channels = I.channels();

		int nRows = I.rows;
		int nCols = I.cols * channels;

		if (I.isContinuous()) {
			nCols *= nRows;
			nRows = 1;
		}

		uchar *p;
		for (int i = 0; i < nRows; i++)
		{
			p = I.ptr<uchar>(i);
			for (int j = 0; j < nCols; j++)
			{
				p[j] = table[p[j]];
			}
		}
		return I;
	}

	Mat& Tutorial::ScanImageAndReduceIterator(Mat& I, const uchar* table) {
		CV_Assert(I.depth() == CV_8U);

		const int channels = I.channels();
		switch (channels)
		{
		case 1:
		{
			MatIterator_<uchar> it, end;
			for (it = I.begin<uchar>(), end = I.end<uchar>(); it != end; ++it)
				*it = table[*it];
			break;
		}
		case 3:
		{
			MatIterator_<Vec3b> it, end;
			for (it = I.begin<Vec3b>(), end = I.end<Vec3b>(); it != end; ++it) {
				(*it)[0] = table[(*it)[0]];
				(*it)[1] = table[(*it)[1]];
				(*it)[2] = table[(*it)[2]];
			}
		}
		}
		return I;
	}

	Mat& Tutorial::ScanImageAndReduceRandomAccess(Mat& I, const uchar* const table) {
		CV_Assert(I.depth() == CV_8U);

		const int channels = I.channels();
		switch (channels)
		{
		case 1:
		{
			for (int i = 0; i < I.rows; ++i)
				for (int j = 0; j < I.cols; j++)
					I.at<uchar>(i, j) = table[I.at<uchar>(i, j)];
			break;
		}
		case 3:
		{
			Mat_<Vec3b> _I = I;
			for (int i = 0; i < _I.rows; i++)
				for (int j = 0; j < _I.cols; j++)
				{
					_I(i, j)[0] = table[_I(i, j)[0]];
					_I(i, j)[1] = table[_I(i, j)[1]];
					_I(i, j)[2] = table[_I(i, j)[2]];
				}
			I = _I;
			break;
		}
		}
		return I;
	}


	void Tutorial::Sharpen(const Mat& myImg, Mat &Result) {
		CV_Assert(myImg.depth() == CV_8U);
		const int nChannels = myImg.channels();
		Result.create(myImg.size(), myImg.type());
		for (int i = 1; i < myImg.rows - 1; ++i) {
			const uchar* pre = myImg.ptr<uchar>(i - 1);
			const uchar* cur = myImg.ptr<uchar>(i);
			const uchar* nex = myImg.ptr<uchar>(i + 1);
			uchar * out = Result.ptr<uchar>(i);
			for (int j = 1; j < (myImg.cols - 1)*nChannels; ++j) {
				*out++ = saturate_cast<uchar>(5 * cur[j] - cur[j - 1] - cur[j + 1] - pre[j] - nex[j]);
			}

			Result.row(0).setTo(Scalar(0));
			Result.col(0).setTo(Scalar(0));
			Result.row(Result.rows - 1).setTo(Scalar(0));
			Result.col(Result.cols - 1).setTo(Scalar(0));
		}
	}

	void Tutorial::MyEllipse(Mat img, double angle)
	{
		int thickness = 2;
		int lineType = 8;

		ellipse(img,
			Point(w / 2, w / 2),
			Size(w / 4, w / 16),
			angle,
			0,
			360,
			Scalar(255, 0, 0),
			thickness,
			lineType);
	}

	void Tutorial::MyFilledCircle(Mat img, Point center)
	{
		circle(img,
			center,
			w / 32,
			Scalar(0, 0, 255),
			FILLED,
			LINE_8);
	}

	void Tutorial::MyPolygon(Mat img)
	{
		int lineType = LINE_8;

		Point rook_points[1][20];
		rook_points[0][0] = Point(w / 4, 7 * w / 8);
		rook_points[0][1] = Point(3 * w / 4, 7 * w / 8);
		rook_points[0][2] = Point(3 * w / 4, 13 * w / 16);
		rook_points[0][3] = Point(11 * w / 16, 13 * w / 16);
		rook_points[0][4] = Point(19 * w / 32, 3 * w / 8);
		rook_points[0][5] = Point(3 * w / 4, 3 * w / 8);
		rook_points[0][6] = Point(3 * w / 4, w / 8);
		rook_points[0][7] = Point(26 * w / 40, w / 8);
		rook_points[0][8] = Point(26 * w / 40, w / 4);
		rook_points[0][9] = Point(22 * w / 40, w / 4);
		rook_points[0][10] = Point(22 * w / 40, w / 8);
		rook_points[0][11] = Point(18 * w / 40, w / 8);
		rook_points[0][12] = Point(18 * w / 40, w / 4);
		rook_points[0][13] = Point(14 * w / 40, w / 4);
		rook_points[0][14] = Point(14 * w / 40, w / 8);
		rook_points[0][15] = Point(w / 4, w / 8);
		rook_points[0][16] = Point(w / 4, 3 * w / 8);
		rook_points[0][17] = Point(13 * w / 32, 3 * w / 8);
		rook_points[0][18] = Point(5 * w / 16, 13 * w / 16);
		rook_points[0][19] = Point(w / 4, 13 * w / 16);

		const Point* ppt[1] = { rook_points[0] };
		int npt[] = { 20 };

		fillPoly(img,
			ppt,
			npt,
			1,
			Scalar(255, 255, 255),
			lineType);
	}

	void Tutorial::MyLine(Mat img, Point start, Point end)
	{
		int thickness = 2;
		int lineType = LINE_4;

		line(img,
			start,
			end,
			Scalar(0, 0, 0),
			thickness,
			lineType);
	}
}