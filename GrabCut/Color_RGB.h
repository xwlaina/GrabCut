/*
 * GrabCut implementation source code Copyright(c) 2005-2006 Justin Talbot
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 */

#ifndef COLOR_RGB_H
#define COLOR_RGB_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class Color_RGB
{

public:

        Color_RGB(): r(0), g(0), b(0){}
        Color_RGB(double _r, double _g, double _b): r(_r), g(_g), b(_b){}
		Color_RGB(const Color_RGB& RightSides): r(RightSides.r), g(RightSides.g), b(RightSides.b){}


		Color_RGB& operator = (const Color_RGB& RightSides)
		{
			r = RightSides.r;
			g = RightSides.g;
			b = RightSides.b;
			return *this;
		}

		const Color_RGB operator*(double alpha)const
        {
            return Color_RGB(r*alpha,g*alpha,b*alpha);
        }

        double r, g, b;
};

// Compute squared distance between two colors
inline double distance2(const Color_RGB &c1, const Color_RGB &c2)
{
        return ((c1.r - c2.r) * (c1.r - c2.r) + 
		(c1.g - c2.g) * (c1.g - c2.g) + 
		(c1.b - c2.b) * (c1.b - c2.b));
}
#endif //COLOR_RGB_H
