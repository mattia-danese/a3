/*  =================== File Information =================
	File Name: ppm.h
	Description:
	Author: Michael Shah
	Last Modified: 4/2/14

	Purpose: 
	Usage:	

	Further Reading resources: 
	PPM Format:
	Specs from: http://netpbm.sourceforge.net/doc/ppm.html
	===================================================== */

#ifndef PPM_H
#define PPM_H

#include <string>
#include "./scene/SceneData.h"

/*
	A ppm is a simple image format.
	It consists of a header with some information about the file to be parsed.
	The rest of the file is a tuple of red(r),green(g), and blue(b) values (you will commonly see 
	abbreviated as RGB.  Other color systems such as Cyan, Magenta, and Yellow may be used by other systems).

	The way OpenGL reads textures, is by taking in an array of character values(that is, integer values 0-255) in a 
	3-tuple(red,green,blue).  Generally OpenGL also likes the dimensions of the image to be in powers of two(32x32, 64x64,etc)
	.

*/
class ppm{
	public:
		/*	===============================================
		Desc:	Default constructor for a ppm
		Precondition: _fileName is the image file name. It is also expected that the file is of type "ppm" 
		Postcondition: The array 'color' is allocated memory according to the image dimensions.
						width and height private members are set based on ppm header information.
		=============================================== */ 
		ppm(std::string _fileName);
		/*	===============================================
		Desc:	Default destructor for a ppm
		Precondition: 
		Postcondition: 'color' array memory is deleted,
		=============================================== */ 
		~ppm();
		/*	===============================================
		Desc:	Draws a 2D quad on the front of the screen at window coordinates x and y
		Precondition: 
		Postcondition:
		=============================================== */ 
		void render(int x, int y);
		/*	===============================================
		Desc:	Sets a pixel in our array a specific color
		Precondition: 
		Postcondition:
		=============================================== */ 
		void setPixel(int x, int y, int r, int g, int b);
		SceneColor getPixel(float s, float t);

		// Getter functions
		int getWidth() { return width;}
		int getHeight() { return height;}
		unsigned char* getPixels() { return color;}
	private:
		std::string magicNumber;	// Used in the header to determine
									// how to parse this file. Example, P3, P6, etc.
		int width;		
		int height;
		unsigned char* color;				// dynamic Array that stores color values.
									// This is a single dimensional array storing a tuple
									// of color values.  We make it a char, because we need to
									// store values of 0-255 (Thus a char, which is generally 1 byte).
									//
									// color[0] = first r value
									// color[1] = first g value
									// color[2] = first b value
									// color[3] = second r value
									// color[4] = second g value
									// color[5] = second b value
									// etc.

		
};

#endif