/*
 * GrabCut implementation source code Copyright(c) 2005-2006 Justin Talbot
 *
 * All Rights Reserved.
 * For educational use only; commercial use expressly forbidden.
 * NO WARRANTY, express or implied, for this software.
 */

#ifndef COLOR_Lab_H
#define COLOR_Lab_H

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

class Color_Lab
{

public:

        Color_Lab(): L(0), a(0), b(0){}
        Color_Lab(double _L, double _a, double _b): L(_L), a(_a), b(_b){}
		Color_Lab(const Color_Lab& RightSides): L(RightSides.L), a(RightSides.a), b(RightSides.b){}


		Color_Lab& operator = (const Color_Lab& RightSides)
		{
			L = RightSides.L;
			a = RightSides.a;
			b = RightSides.b;
			return *this;
		}

        double L, a, b;
};

// Compute squared distance between two colors
inline double distance2(const Color_Lab &c1, const Color_Lab &c2)
{
        return ((c1.L - c2.L) * (c1.L - c2.L) + 
		(c1.a - c2.a) * (c1.a - c2.a) + 
		(c1.b - c2.b) * (c1.b - c2.b));
}
#endif //COLOR_Lab_H
