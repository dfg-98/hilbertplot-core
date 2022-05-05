#ifndef HPOINT_H
#define HPOINT_H

#include <iostream>
#include <fstream>
#include <vector>
#include "hilbertdefines.h"


//..............................................................
// Defines the coordinates of a point in the 2D plane of a HilbertCurve
//..............................................................
class HPoint
{
    public:
        //Constructors
        HPoint();
        HPoint(hint xx, hint yy);
        HPoint(hint n);
        HPoint(const HPoint & p);

        //Getters and Setter
        hint X() const;
        hint Y() const;
        hint X(hint xx);
        hint Y(hint yy);
        hfloat DifferenceValue() const;

        //Operator overloading
        HPoint operator+(const HPoint & p);
        HPoint operator-(const HPoint & p);
        HPoint operator/(hint n);
        HPoint operator%(hint n);

        HPoint & operator+=(const HPoint & p);
        HPoint & operator-=(const HPoint & p);
        HPoint & operator=(const HPoint & p);

        HPoint & operator+=(hint n);
        HPoint & operator-=(hint n);
        HPoint & operator=(hint n);
        HPoint & operator*=(hint n);
        HPoint & operator/=(hint n);
        HPoint & operator%=(hint n);

        friend HPoint operator*(const HPoint & p, hint n);
        friend HPoint operator*(hint n, const HPoint & p);

        friend bool operator==(const HPoint & p1, const HPoint & p2);
        friend bool operator!=(const HPoint & p1, const HPoint & p2);
        friend bool operator>(const HPoint & p1, const HPoint & p2);
        friend bool operator>=(const HPoint & p1, const HPoint & p2);
        friend bool operator<(const HPoint & p1, const HPoint & p2);
        friend bool operator<=(const HPoint & p1, const HPoint & p2);
        friend bool indexCmp(const HPoint &p1, const HPoint & p2);

        friend std::ostream & operator<<(std::ostream &, const HPoint & p);
        friend std::istream & operator>>(std::istream &, HPoint & p);
        friend std::ofstream & write(std::ofstream &out, const HPoint &p);
        friend std::ifstream & read(std::ifstream &out, HPoint &p);

        friend class HilbertCurve;
        friend class HilbertPlot;
        friend class QuasiSquare;

    protected:
        //Coordinates
        hint x;
        hint y;
        hfloat difference;
        hint index;
};

bool operator==(const HPoint & p1, const HPoint & p2);
bool operator!=(const HPoint & p1, const HPoint & p2);
bool operator>(const HPoint & p1, const HPoint & p2);
bool operator>=(const HPoint & p1, const HPoint & p2);
bool operator<(const HPoint & p1, const HPoint & p2);
bool operator<=(const HPoint & p1, const HPoint & p2);
bool indexCmp(const HPoint &p1, const HPoint &p2);
std::ostream & operator<<(std::ostream & out, const HPoint & p);
std::istream & operator>>(std::istream & is, HPoint & p);
std::ifstream &read(std::ifstream &in, HPoint &p);
std::ofstream & write(std::ofstream &out, const HPoint &p);


#endif // HPOINT_H
