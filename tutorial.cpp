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
using namespace std;
using namespace cimg_library;

#ifndef cimg_imagepath
#define cimg_imagepath "img/HDRsequence/"
#endif

struct ImageInfo {
	std::string imageName;
	float t;
};

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
	return imgVector;
}

int main(int argc,char **argv) {

	//hdrgen datei einlesen
	vector<ImageInfo> imageNames = readFile("img/HDRsequence/max.hdrgen");

	for (std::vector<ImageInfo>::const_iterator i = imageNames.begin(); i != imageNames.end(); ++i)
		std::cout << (*i).imageName << ' ' << (*i).t << "\n";

	std::cout << imageNames.size() << std::endl;

	//Alle bilder einlesen
	std::vector<CImg<unsigned char>> images;
	for (int i = 0; i < imageNames.size(); i++) {
		std::string imageName = "img/HDRsequence/" + imageNames[i].imageName;
		images.push_back(CImg<>(imageName.c_str()).normalize(0, 255));
	}

	//Bild was angezeigt wird
	const CImg<unsigned char> image = images[0];


  // Create two display window, one for the image, the other for the color profile.
  CImgDisplay
    main_disp(image,"Color image (Try to move mouse pointer over)",0),
    draw_disp(500,400,"Color profile of the X-axis",0);

  // Define colors used to plot the profile, and a hatch to draw the vertical line
  unsigned int hatch = 0xF0F0F0F0;
  const unsigned char
    red[]   = { 255,0,0 },
    green[] = { 0,255,0 },
    blue [] = { 0,0,255 },
    black[] = { 0,0,0 };

  	

    // Enter event loop. This loop ends when one of the two display window is closed or
    // when the keys 'ESC' or 'Q' are pressed.
    while (!main_disp.is_closed() && !draw_disp.is_closed() &&
           !main_disp.is_keyESC() && !draw_disp.is_keyESC() && !main_disp.is_keyQ() && !draw_disp.is_keyQ()) {

      // Handle display window resizing (if any)
      if (main_disp.is_resized()) main_disp.resize().display(image);
      draw_disp.resize();

      if (main_disp.mouse_x()>=0 && main_disp.mouse_y()>=0) { // Mouse pointer is over the image

        const int
          xm = main_disp.mouse_x(),                     // X-coordinate of the mouse pointer over the image
          ym = main_disp.mouse_y(),                     // Y-coordinate of the mouse pointer over the image
          xl = xm*draw_disp.width()/main_disp.width(),  // Corresponding X-coordinate of the hatched line
          x = xm*image.width()/main_disp.width(),     // Corresponding X-coordinate of the pointed pixel in the image
          y = ym*image.height()/main_disp.height();   // Corresponding Y-coordinate of the pointex pixel in the image

        // Retrieve color component values at pixel (x,y)
        const unsigned int
          val_red   = image(x,y,0),
          val_green = image(x,y,1),
          val_blue  = image(x,y,2);

        // Create and display the image of the intensity profile
        CImg<unsigned char>(draw_disp.width(),draw_disp.height(),1,3,255).
          draw_grid(-50*100.0f/image.width(),-50*100.0f/256,0,0,false,true,black,0.2f,0xCCCCCCCC,0xCCCCCCCC).
          draw_axes(0,image.width() - 1.0f,255.0f,0.0f,black).
          draw_graph(image.get_shared_row(y,0,0),red,1,1,0,255,1).
          draw_graph(image.get_shared_row(y,0,1),green,1,1,0,255,1).
          draw_graph(image.get_shared_row(y,0,2),blue,1,1,0,255,1).
          draw_text(30,5,"Pixel (%d,%d)={%d %d %d}",black,0,1,16,
                    main_disp.mouse_x(),main_disp.mouse_y(),val_red,val_green,val_blue).
          draw_line(xl,0,xl,draw_disp.height() - 1,black,0.5f,hatch=cimg::rol(hatch)).
          display(draw_disp);
      } else
        // else display a text in the profile display window.
        CImg<unsigned char>(draw_disp.width(),draw_disp.height()).fill(255).
          draw_text(draw_disp.width()/2 - 130,draw_disp.height()/2 - 5,"Mouse pointer is outside the image",
                    black,0,1,16).display(draw_disp);

      // Temporize event loop
      cimg::wait(20);
    }

    return 0;
}
