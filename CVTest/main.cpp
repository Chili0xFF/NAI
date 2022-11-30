#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <iostream>

int main(int argc, char** argv)
{
    using namespace cv;
    VideoCapture vc(0);
    if (!vc.isOpened() ) return 1;
    cv::Mat frame;
    int i = 0;
    while (waitKey(10) != 27) {
        vc >> frame;
        cv::Mat src = frame;
        cv::Mat dst;
        cv::flip(src,dst,1);
        cv::imshow("Frame",dst);

    }
    return 0;
}