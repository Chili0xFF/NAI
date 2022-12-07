#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <iostream>

int main(int argc, char** argv)
{
    using namespace cv;
    using namespace std;
    VideoCapture vc(0);
    if (!vc.isOpened() ) return 1;
    Mat frame, src, dst, detected, dilated, lined;
    std::vector<int> lower = {110,95,30};
    std::vector<int> upper = {174,180,80};

    namedWindow("pierwsze", WINDOW_AUTOSIZE);
    namedWindow("detected", WINDOW_AUTOSIZE);
    createTrackbar("lh", "detected", &lower[0], 255);
    createTrackbar("ls", "detected", &lower[1], 255);
    createTrackbar("lv", "detected", &lower[2], 255);
    createTrackbar("hh", "detected", &upper[0], 255);
    createTrackbar("hs", "detected", &upper[1], 255);
    createTrackbar("hv", "detected", &upper[2], 255);

    while (waitKey(10) != 27) {
        vc >> frame;
        src = frame;
        flip(src,dst,1);
        imshow("Original",dst);
        inRange(dst, lower,upper, detected);
        imshow("detected",detected);
        auto kernel = getStructuringElement(MORPH_ELLIPSE,Size{5,5});
        dilate(detected, dilated, kernel);
        erode(dilated, dilated, kernel);
        imshow("dilated", dilated);
        vector<vector<Size>> contours,restOfContours;
        vector<Vec4i> hierarchy;
        findContours (dilated, contours, hierarchy, RETR_LIST, CHAIN_APPROX_SIMPLE);
        //cout << contours.size() << endl;
        std::sort(contours.begin(), contours.end(),[](auto a, auto b){
            return contourArea(a) > contourArea(b);
        });

        cout << "Object A: " << contours[0][0] << endl;
        cout << "Object B: " << contours[1][0] << endl;
        lined = dst;
        line(lined,contours[0][0],contours[1][0],Scalar(0,0,255),3,LINE_AA);
        imshow("lined",lined);
    }
    return 0;
}
