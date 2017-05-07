/*
#
#  File        : tutorial.cpp
#                ( C++ source file )
#
#  Description : View the color profile of an image, along the X-axis.
#                This file is a part of the CImg Library project.
#                ( http://cimg.eu )
#
#  Copyright   : David Tschumperle
#                ( http://tschumperle.users.greyc.fr/ )
#
#  License     : CeCILL v2.0
#                ( http://www.cecill.info/licences/Licence_CeCILL_V2-en.html )
#
#  This software is governed by the CeCILL  license under French law and
#  abiding by the rules of distribution of free software.  You can  use,
#  modify and/ or redistribute the software under the terms of the CeCILL
#  license as circulated by CEA, CNRS and INRIA at the following URL
#  "http://www.cecill.info".
#
#  As a counterpart to the access to the source code and  rights to copy,
#  modify and redistribute granted by the license, users are provided only
#  with a limited warranty  and the software's author,  the holder of the
#  economic rights,  and the successive licensors  have only  limited
#  liability.
#
#  In this respect, the user's attention is drawn to the risks associated
#  with loading,  using,  modifying and/or developing or reproducing the
#  software by the user in light of its specific status of free software,
#  that may mean  that it is complicated to manipulate,  and  that  also
#  therefore means  that it is reserved for developers  and  experienced
#  professionals having in-depth computer knowledge. Users are therefore
#  encouraged to load and test the software's suitability as regards their
#  requirements in conditions enabling the security of their systems and/or
#  data to be ensured and,  more generally, to use and operate it in the
#  same conditions as regards security.
#
#  The fact that you are presently reading this means that you have had
#  knowledge of the CeCILL license and that you accept its terms.
#
*/

// Include CImg library file and use its main namespace
#include "CImg.h"
#include <iostream>
#include <iomanip>
#include <fstream>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <map>

using namespace std;
using namespace cimg_library;

#ifndef cimg_imagepath
#define cimg_imagepath "img/HDRsequence/"
#endif

const double EulerConstant = std::exp(1.0);

typedef struct ImageInfo {
	std::string imageName;
	float t;
};

typedef struct rgb {
	float r;
	float g;
	float b;

	rgb() :r(0.f),
		g(0.f),
		b(0.f) {}
	rgb(float _r, float _g, float _b) :r(_r),
		g(_g),
		b(_b) {}
};

typedef map<float, map<int, array<rgb, 256>>> CardinalityMap;



// Main procedure
//----------------
//void main()
//{

//}
vector<ImageInfo> readFile(std::string filename) {
	string string;
	std::ifstream in(filename);
	std::string line;
	if (!in.is_open()) {
		cout << "loadHDRGEN: can not open " << endl;
	}

	vector<ImageInfo> imgVector;
	std::string img;
	float t;
	const int MAX = 256;
	for (std::string img; in >> setw(MAX) >> img; ) {
		in >> setw(MAX) >> t;
		ImageInfo info;
		info.imageName = img;
		info.t = t;
		imgVector.push_back(info);
	}

	in.close();

	/*std::sort(imgVector.begin(), imgVector.end(),
	[](const ImageInfo& lhs, const ImageInfo& rhs) {
	return lhs.t > rhs.t;
	});*/

	return imgVector;
}

void writeToCSV(vector<rgb> I) {
	ofstream myfile;
	string r = "", g = "", b = "";
	myfile.open("responseCurve.csv");
	for (int i = 0; i < I.size(); i++) {
		r = std::to_string(I[i].r);
		g = std::to_string(I[i].g);
		b = std::to_string(I[i].b);
		myfile << std::to_string(i) << "," << r << "," << g << "," << b << "\n";
	}
	myfile.close();
	return;
}

void plotGnuPlot() {
	system("gnuplot -p -e \"load 'responseCurve.p'\"");
}

//weight function
float weight(float y) {
	if (y > 255.f) cout << "big y: " << y << endl;
	if (y <= 0.f || y >= 255.f) return 0;
	else return std::exp(-4.f * (powf(y - 127.5f, 2) / powf(127.5f, 2)));
}

//xj
CImg<float> calculate_irradiance(vector<CImg<float>> images, vector<ImageInfo> imgInfo, vector<rgb> I) {
	CImg<float> result(images[0]); //bild mit gleichen dimensionen anlegen
	result.fill(0);
	vector<float> numerator(images[0].size(), 0.0f);
	vector<float> sum(images[0].size(), 0.0f);
	float Iy, t, w;
	CImg<float> img;

	for (int i = 0; i < images.size(); i++) {
		img = images[i];
		t = 1.0f / imgInfo[i].t; //exposure time

		cimg_foroff(img, j) { //iteriert über bild buffer
			w = weight((float)img[j]);

			//response function value for pixel value
			if (j < img.size() / 3) { //Red values
				Iy = I[img[j]].r;
			}
			else if (j < 2 * img.size() / 3) { //Green values 
				Iy = I[img[j]].g;
			}
			else { //Blue values
				Iy = I[img[j]].b;
			}

			//von paper
			numerator[j] += (w * t * Iy);
			sum[j] += (w * powf(t, 2));
			//von folien
			//result[j] += (w * powf(t, 2) * (Iy / t)) / (w * powf(t, 2));
		}
	}

	for (int y = 0; y < result.size(); y++) {
		result[y] = numerator[y] / sum[y];
	}


	return result;
}

vector<rgb> calculate_response_curve(vector<CImg<float>> images, vector<ImageInfo> imgInfo, CImg<float> x) {
	vector<rgb> I(256, { 0.f,0.f,0.f });
	CImg<float> img;
	float t;

	for (int m = 0; m < I.size(); m++) {
		rgb Card = { 0.f,0.f,0.f };
		for (int i = 0; i < images.size(); i++) {
			img = images[i];
			t = 1.0f / imgInfo[i].t;
			//if (std::equalimgInfo[i].imageName)
			cimg_foroff(img, j) { //iteriert über bild buffer
				if ((int)img[j] == m) {
					//if ((int)x[j] > 100) cout << (int)x[j] << endl;
					if (j < img.size() / 3) { //Red values
						I[m].r += t * x[j];
						Card.r += 1;
					}
					else if (j < 2 * img.size() / 3) { //Green values 
						I[m].g += t * x[j];
						Card.g += 1;
					}
					else { //Blue values
						I[m].b += t * x[j];
						Card.b += 1;
					}
				}
			}
		}
		//cout << "Card.r: " << Card.r << endl;
		//cout << "I[m].r: " << I[m].r << endl;
		I[m].r /= Card.r;
		I[m].g /= Card.g;
		I[m].b /= Card.b;
	}

	//normalisieren sodass I_128 = 1.0f ist
	rgb i128 = I[128];
	for (int i = 0; i < 256; i++) {
		I[i].r /= i128.r;
		I[i].g /= i128.g;
		I[i].b /= i128.b;
	}


	return I;
}

rgb calculate_objective_f(vector<CImg<float>> images, vector<ImageInfo> imgInfo, vector<rgb> I, CImg<float> x) {
	rgb result = { 0.f,0.f,0.f };
	float t, w;
	CImg<float> img;

	for (int i = 0; i < images.size(); i++) {
		t = 1.0f / imgInfo[i].t;
		img = images[i];
		cimg_foroff(img, j) { //iteriert über bild buffer
			if (j < img.size() / 3) { //Red values
				result.r += weight(img[j]) * (powf(I[img[j]].r - (t*x[j]), 2));
			}
			else if (j < 2 * img.size() / 3) { //Green values 
				result.g += weight(img[j]) * ( powf( I[ img[j] ].g - (t*x[j]), 2) );
			}
			else { //Blue values
				result.b += weight(img[j]) * (powf(I[img[j]].g - (t*x[j]), 2));
			}
		}
	}

	return result;
}


bool imageInfoComparison(ImageInfo i, ImageInfo j) {
	cout << (i.t == j.t) << endl;
	return (i.t == j.t);
}

int main(int argc, char **argv) {

	//hdrgen datei einlesen
	// dim: 308 x 200
	// 61600 Buffer Einträge pro Farbkanal = 61600 x 3 = 184800
	vector<ImageInfo> imageNames = readFile("img/HDRsequence/max.hdrgen");
	//vector<ImageInfo> imageNames = readFile("img/fenster/d/d.hdrgen");

	for (std::vector<ImageInfo>::const_iterator i = imageNames.begin(); i != imageNames.end(); ++i)
		std::cout << (*i).imageName << ' ' << (*i).t << "\n";

	std::cout << imageNames.size() << std::endl;

	//Alle bilder einlesen
	std::vector<CImg<float>> images;
	for (int i = 0; i < imageNames.size(); i++) {
		std::string imageName = "img/HDRsequence/" + imageNames[i].imageName;
		//std::string imageName = "img/fenster/d/" + imageNames[i].imageName;
		images.push_back(CImg<>(imageName.c_str(	)));
	}

	//einträge mit gleicher exposure time löschen, um cardinalität anzulegen
	/*std::vector<ImageInfo>::iterator it;
	it = std::unique(imageNames.begin(), imageNames.end(), imageInfoComparison);
	imageNames.resize(std::distance(imageNames.begin(), it));*/

	//initiale response curve mit I_1 = 1/128, I_128 = 1, I_256 = 2f, 
	vector<rgb> I(256);
	for (int i = 0; i < 256; i++) {
		float v = (i) / 128.f;
		I[i] = { v,v,v };
	}
	CImg<float> x = calculate_irradiance(images, imageNames, I);

	rgb o = calculate_objective_f(images, imageNames, I, x);
	cout << "o.rgb : " << o.r << ", " << o.g << ", " << o.b  << endl;

	//TODO convergenz durch delta oder so in while schleife abfragen
	for (int n = 0; n < 5; n++) {
		cout << "Iteration: " << n << endl;
		I = calculate_response_curve(images, imageNames, x);
		x = calculate_irradiance(images, imageNames, I);

		o = calculate_objective_f(images, imageNames, I, x);
		cout << "o.rgb : " << o.r << ", " << o.g << ", " << o.b << endl;
	}
	
	CImg<float> image(x);

	for (int i = 0; i < image.size(); i++) {
		image[i] = x[i];
	}
	image.normalize(0, 255);
	//test ob die if abfragen die richtigen farbkanäle ansprechen
	//
	//CImg<unsigned char> image(images[0]);
	//for (int j = 0; j < image.size(); j++) {
	//	if (j < image.size() / 3) { //Red values
	//		image[j] = images[0][j] * 2;
	//	}
	//	else if (j < 2 * image.size() / 3) { //Green values 
	//		image[j] = images[0][j] * 0.5f;
	//	}
	//	else { //Blue values
	//		image[j] = images[0][j] * 0.25f;
	//	}
	//}

	//test der weight funktion
	//
	/*for (int j = 0; j < 256; j++) {
	cout << j  << ": " << weight(j) << endl;
	}*/

	writeToCSV(I);
	//plotGnuPlot(); wirst du nicht ausführen können, musst gnuplot installieren :S

	//std::string imageName = "img/fenster/d/d0.ppm";
	//image = (CImg<>(imageName.c_str()));

	// Create two display window, one for the image, the other for the color profile.
	CImgDisplay
		main_disp(image, "Color image (Try to move mouse pointer over)", 0),
		draw_disp(500, 400, "Color profile of the X-axis", 0);

	// Define colors used to plot the profile, and a hatch to draw the vertical line
	unsigned int hatch = 0xF0F0F0F0;
	const unsigned char
		red[] = { 255,0,0 },
		green[] = { 0,255,0 },
		blue[] = { 0,0,255 },
		black[] = { 0,0,0 };



	// Enter event loop. This loop ends when one of the two display window is closed or
	// when the keys 'ESC' or 'Q' are pressed.
	while (!main_disp.is_closed() && !draw_disp.is_closed() &&
		!main_disp.is_keyESC() && !draw_disp.is_keyESC() && !main_disp.is_keyQ() && !draw_disp.is_keyQ()) {

		// Handle display window resizing (if any)
		if (main_disp.is_resized()) main_disp.resize().display(image);
		draw_disp.resize();

		if (main_disp.mouse_x() >= 0 && main_disp.mouse_y() >= 0) { // Mouse pointer is over the image

			const int
				xm = main_disp.mouse_x(),                     // X-coordinate of the mouse pointer over the image
				ym = main_disp.mouse_y(),                     // Y-coordinate of the mouse pointer over the image
				xl = xm*draw_disp.width() / main_disp.width(),  // Corresponding X-coordinate of the hatched line
				x = xm*image.width() / main_disp.width(),     // Corresponding X-coordinate of the pointed pixel in the image
				y = ym*image.height() / main_disp.height();   // Corresponding Y-coordinate of the pointex pixel in the image

															  // Retrieve color component values at pixel (x,y)
			const unsigned int
				val_red = image(x, y, 0),
				val_green = image(x, y, 1),
				val_blue = image(x, y, 2);

			// Create and display the image of the intensity profile
			CImg<unsigned char>(draw_disp.width(), draw_disp.height(), 1, 3, 255).
				draw_grid(-50 * 100.0f / image.width(), -50 * 100.0f / 256, 0, 0, false, true, black, 0.2f, 0xCCCCCCCC, 0xCCCCCCCC).
				draw_axes(0, image.width() - 1.0f, 255.0f, 0.0f, black).
				draw_graph(image.get_shared_row(y, 0, 0), red, 1, 1, 0, 255, 1).
				draw_graph(image.get_shared_row(y, 0, 1), green, 1, 1, 0, 255, 1).
				draw_graph(image.get_shared_row(y, 0, 2), blue, 1, 1, 0, 255, 1).
				draw_text(30, 5, "Pixel (%d,%d)={%d %d %d}", black, 0, 1, 16,
					main_disp.mouse_x(), main_disp.mouse_y(), val_red, val_green, val_blue).
				draw_line(xl, 0, xl, draw_disp.height() - 1, black, 0.5f, hatch = cimg::rol(hatch)).
				display(draw_disp);
		}
		else
			// else display a text in the profile display window.
			CImg<unsigned char>(draw_disp.width(), draw_disp.height()).fill(255).
			draw_text(draw_disp.width() / 2 - 130, draw_disp.height() / 2 - 5, "Mouse pointer is outside the image",
				black, 0, 1, 16).display(draw_disp);

		// Temporize event loop
		cimg::wait(20);
	}

	return 0;
}