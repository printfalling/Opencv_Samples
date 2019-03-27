#include "Tutorial.h"



int main()
{
	Tutorial::Tutorial *tu = new Tutorial::Tutorial();
	const int test_no = 6;

	int argc;
	const char** argv;

	switch (test_no)
	{
	case 0:
		//C:\Users\lenobo\Desktop\外观OK2\A52_2_1.bmp 32 G
		argc = 4;
		argv = new const char*[4];
		argv[0] = "how_to_scan_images";
		argv[1] = "C:\\Users\\lenobo\\Desktop\\外观OK2\\A52_2_1.bmp";
		argv[2] = "32";
		argv[3] = "G";
		
		return tu->how_to_scan_images(argc, argv);

	case 1:
		argc = 3;
		argv = new const char*[3];
		argv[0] = "mask_operate";
		argv[1] = "C:\\Users\\lenobo\\Desktop\\外观OK2\\A52_2_1.bmp";
		argv[2] = "G";

		return tu->mask_operation(argc, argv);

	case 2:
		argc = 3;
		argv = new const char*[3];
		argv[0] = "add_liner";
		argv[1] = "linux_logo.jpg";
		argv[2] = "windows_logo.jpg";

		return tu->add_image(argc, argv);

	case 3:
		argc = 2;
		argv = new const char*[2];
		argv[0] = "change constract";
		argv[1] = "windows_logo.jpg";

		return tu->change_contract(argc, argv);

	case 4:
		argc = 2;
		argv = new const char*[2];
		argv[0] = "Discrete Fourier Transform";
		argv[1] = "imageTextN.png";

		tu->discrete_fourier_transformation(argc, argv);
		argv[1] = "imageTextR.png";
		tu->discrete_fourier_transformation(argc, argv);

	case 5:
		argc = 2;
		argv = new const char*[2];
		argv[0] = "xml_yaml_in_out";
		argv[1] = "I.xml";

		tu->xml_yaml_in_out(argc, argv);

	case 6:
		argc = 1;
		argv = new const char*[1];
		argv[0] = "basic_draw";

		tu->basic_draw(argc, argv);
	default:
		return -1;
	}
}
