/*  =================== File Information =================
	File Name: ppm.cpp
	Description:
	Author: Michael Shah
	Last Modified: 4/2/14

	Purpose: 
	Usage:	
	===================================================== */

#include <FL/gl.h>
#include <iostream>
#include <string>
#include <fstream>
#include "ppm.h"

/*	===============================================
Desc:	Default constructor for a ppm
Precondition: _fileName is the image file name. It is also expected that the file is of type "ppm"
              It is expected that width and height are also set in the constructor. 
Postcondition: The array 'color' is allocated memory according to the image dimensions.
				width and height private members are set based on ppm header information.
=============================================== */ 
ppm::ppm(std::string _fileName){
  /* Algorithm
      Step 1: Parse header of PPM
      Step 2: Read in colors into array
      Step 3: Allocate memory for width and height dimensions
  */
  // Open an input file stream for reading a file
  std::ifstream ppmFile(_fileName.c_str());
  // If our file successfully opens, begin to process it.
  if (ppmFile.is_open()) {
	  // line will store one line of input
	  std::string line;
	  // Our loop invariant is to continue reading input until
	  // we reach the end of the file and it reads in a NULL character
	  std::cout << "Reading in ppm file: " << _fileName << std::endl;
	  // Our delimeter pointer which is used to store a single token in a given
	  // string split up by some delimeter(s) specified in the strtok function
	  char* delimeter_pointer;
	  int iteration = 0;
	  int pos = 0;
	  while (getline(ppmFile, line)) {
		  char* copy = new char[line.length() + 1];
		  strcpy(copy, line.c_str());
		  delimeter_pointer = strtok(copy, " ");

		  // Read in the magic number
		  if (iteration == 0) {
			  magicNumber = delimeter_pointer;
			  std::cout << "Magic Number: " << magicNumber << " ";
			  std::cout << std::endl;
		  }
		  else if (iteration == 1) {
			  std::cout << line << std::endl;
		  }
		  // Read in dimensions
		  else if (iteration == 2) {
			  width = atoi(delimeter_pointer);
			  std::cout << "width: " << width << " ";
			  delimeter_pointer = strtok(NULL, " ");
			  height = atoi(delimeter_pointer);
			  std::cout << "height: " << height << std::endl;
			  // Allocate memory for the color array
			  if (width > 0 && height > 0) {
				  color = new unsigned char[width*height * 3];
				  if (color == NULL) {
					  std::cout << "Unable to allocate memory for ppm" << std::endl;
					  exit(1);
				  }
			  }
			  else {
				  std::cout << "PPM not parsed correctly, width and height dimensions are 0" << std::endl;
				  exit(1);
			  }
		  }
		  else if (iteration == 3) {
			  std::cout << "color range: 0-" << delimeter_pointer << std::endl;
		  }
		  else {
			  // Iterate through the entire line and begin storing values
			  while (delimeter_pointer != NULL) {
				  //std::cout << pos << " ";
				  int value = atoi(delimeter_pointer);
				  color[pos] = (unsigned char)value;
				  pos++;
				  delimeter_pointer = strtok(NULL, " ");
			  }
		  }
		  delete[] copy;
		  iteration++;
	  }
	std::cout << " Finished ppm mapping ";
	  ppmFile.close();
  }
  else{
      std::cout << "Unable to open ppm file" << std::endl;
  } 

}



/*	===============================================
Desc:	Default destructor for a ppm
Precondition: 
Postcondition: 'color' array memory is deleted,
=============================================== */ 
ppm::~ppm(){
  if(color!=NULL){
    delete[] color;
    color=0;
  }	
}

/*	===============================================
Desc:	Draws a 2D quad on the front of the screen at window coordinates x and y
Precondition: 
Postcondition:
=============================================== */ 
void ppm::render(int x, int y){
	 
}

/*  ===============================================
Desc: Sets a pixel in our array a specific color
Precondition: 
Postcondition:
=============================================== */ 
void ppm::setPixel(int x, int y, int r, int g, int b){
  if(x > width || y > height){
    return;
  }
  else{
    /*std::cout << "modifying pixel at " 
              << x << "," << y << "from (" <<
              (int)color[x*y] << "," << (int)color[x*y+1] << "," << (int)color[x*y+2] << ")";
*/
    color[(x*3)+height*(y*3)] = r;
    color[(x*3)+height*(y*3)+1] = g;
    color[(x*3)+height*(y*3)+2] = b;
/*
    std::cout << " to (" << (int)color[x*y] << "," << (int)color[x*y+1] << "," << (int)color[x*y+2] << ")" << std::endl;
*/
  }
}

SceneColor ppm::getPixel(float s, float t){
	// std::cout << "s " << s << " t " << t << std::endl;
	// std::cout << "width " << width << "height " << height << std::endl;
	SceneColor ret;
	if(s > width || t > height){
		ret.r = 0;
		ret.g = 0;
		ret.b = 0;
    return ret;
  }else{
	int s_i = (int)s;
	int t_i = (int)t;
// std::cout << "getting pixel at " 
//               << s << "," << t << "from (" <<
//               (int)color[s_i*t_i] << "," << (int)color[s_i*t_i+1] << "," << (int)color[s_i*t_i+2] << ")";
    ret.r = color[s_i*3+width*t_i*3];
    ret.g = color[s_i*3+width*t_i*3+1];
    ret.b = color[s_i*3+width*t_i*3+2];
	return ret;
  }
}

