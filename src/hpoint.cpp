/*!
  \headerfile "hpoit.h"

  \title Hilbert Point Definition

  \brief The hpoint.h header declares the HPoint class.
*/
#include "hpoint.h"




//.....................................................................
//                      Operator overloading
//.....................................................................

// .............................................................................
// Name: operator+=
//
// Synopsis: a=a+b
//
// Parameters:
//			const Pohint2D & p;  ----> second summand
//
// Returns:
//         *this
//
// Exceptions:
//            None
//..............................................................................
/*!
  Sums point \a p and assign it to \c this.
*/
HPoint & HPoint::operator +=(const HPoint & p)
{
    x+=p.x;
    y+=p.y;

    return *this;
}
// .............................................................................
// Name: operator-=
//
// Synopsis: a=a-b
//
// Parameters:
//			const Pohint2D & p;  ----> second summand
//
// Returns:
//         *this
//
// Exceptions:
//            None
//....................................................
/*!
  Substracts point \a p and assign difference to \c this.
*/
HPoint & HPoint::operator -=(const HPoint & p)
{
    x-=p.x;
    y-=p.y;

    return *this;
}
// .............................................................................
// Name: operator=
//
// Synopsis: a=b. Assignment operator
//
// Parameters:
//			const Pohint2D & p;  ----> second summand
//
// Returns:
//         *this
//
// Exceptions:
//            None
//....................................................
/*!
  \brief Assignment operator.
  Assign \a p.
*/
HPoint & HPoint::operator =(const HPoint & p)
{
    x = p.x;
    y = p.y;
    difference = p.difference;
    index = p.index;
    return *this;
}
// .............................................................................
// Name: operator+
//
// Synopsis: a+b. Sum operator
//
// Parameters:
//			const Pohint2D & p;  ----> second summand
//
// Returns:
//         a+b
//
// Exceptions:
//            None
//....................................................
/*!
  Returns the sum of this with \a p.
*/
HPoint HPoint::operator+(const HPoint & p)
{
    HPoint newpohint(*this);
    newpohint+=p;
    return newpohint;
}
// .............................................................................
// Name: operator-
//
// Synopsis: a-b. Diference operator
//
// Parameters:
//			const Pohint2D & p;  ----> second summand
//
// Returns:
//         a-b
//
// Exceptions:
//            None
//....................................................
/*!
  Resturns the differnce between \c this and \a p.
*/
HPoint HPoint::operator-(const HPoint & p)
{
    HPoint newpohint(*this);
    newpohint-=p;
    return newpohint;
}
// .............................................................................
// Name: operator+=
//
// Synopsis: a=a+n. Adds an hinteger to both coordinates
//
// Parameters:
//			hint n;  ----> second summand
//
// Returns:
//         *this
//
// Exceptions:
//            None
//..............................................................................
/*!
  \overload operator+=()
  Add \a n to both components an assign it.
*/
HPoint & HPoint::operator +=(hint n)
{
    x+=n;
    y+=n;
    return *this;
}
// .............................................................................
// Name: operator-=
//
// Synopsis: a=a-n. Substract an hinteger to both coordinates
//
// Parameters:
//          hint n;  ----> second summand
//
// Returns:
//         *this
//
// Exceptions:
//            None
//..............................................................................
/*!
  \overload operator-=()
  Substract \a n from both components and assign it.
*/
HPoint & HPoint::operator -=(hint n)
{
    x-=n;
    y-=n;
    return *this;
}
// .............................................................................
// Name: operator=
//
// Synopsis: a=n. Assign an hinteger to both coordantes.
//
// Parameters:
//			hint n;  ----> second summand
//
// Returns:
//         *this
//
// Exceptions:
//            None
//..............................................................................
/*!
  \overload operator=()
  Assign \a n to both components.
*/
HPoint & HPoint::operator =(hint n)
{
    x=n;
    y=n;
    return *this;
}
// .............................................................................
// Name: operator/=
//
// Synopsis: a=a/n. Divide both coordinate by an hinteger.
//
// Parameters:
//			hint n;  ----> second summand
//
// Returns:
//         *this
//
// Exceptions:
//            None
//..............................................................................
/*!
  Divide both component by \a n and assign it. Did not chek for zero division.
  Standar exception will be rised in that case.
*/
HPoint & HPoint::operator /=(hint n)
{
    x/=n;
    y/=n;
    return *this;
}
// .............................................................................
// Name: operator*=
//
// Synopsis: a=a*n. Multiplies both coordinate by an hinteger.
//
// Parameters:
//			hint n;  ----> second summand
//
// Returns:
//         *this
//
// Exceptions:
//            None
//..............................................................................
/*!
  Multiplies both coordinate by \a n.
*/
HPoint & HPoint::operator*=(hint n)
{
    x*=n;
    y*=n;
    return *this;
}
// .............................................................................
// Name: operator%=
//
// Synopsis: a=a%n. Module operation over each coordinate.
//
// Parameters:
//			hint n;  ----> operand value
//
// Returns:
//         *this
//
// Exceptions:
//            None
//..............................................................................
/*!
 Module operation over each coordinate (a=a% \a n).
*/
HPoint & HPoint::operator%=(hint n)
{
    x%=n;
    y%=n;
    return *this;
}

/*!
  \class HPoint
  \inmodule hilbertlib
  \ingroup hcurve
   \brief Base element of the HilbertCurve.

   The HPoint class represent a point in a 2D space. It has \c x and \c y
   coordinates, an \c index value and a \c difference value. The \c x and \c y
   represent the position on the plane of the curve and th e \c index the linear position.
   So and HPoint could be treated as a numeric ray point or a plane point. The difference
   data member is used for plot the Difference Map (See HilbertPlot).

*/

//.................................................................................
// Constructors
//.................................................................................
/*!
  Default constructs a HPoint at (0, 0) with index value zero and difference map also zero.
*/
HPoint::HPoint():x(0), y(0), difference(0), index(0){}
/*!
  Constructs a HPoint at (\a{xx}, \a{yy}) with index value zero and difference map also zero.
*/
HPoint::HPoint(hint xx, hint yy):x(xx), y(yy), difference(0), index(0){}
/*!
  Constructs a HPoint at (\a n, \a n).
*/
HPoint::HPoint(hint n):x(n), y(n), difference(0), index(0){}
/*!
  Constructs a copy of \a p.
*/
HPoint::HPoint(const HPoint &p):x(p.x), y(p.y), difference(p.difference), index(p.index){}

//---------------------------
// Data memebers accesors
//---------------------------
/*! Returns \c X coordinate.*/
hint HPoint::X() const
{
    return x;
}
/*! Returns \c Y coordinate.*/
hint HPoint::Y() const
{
    return y;
}
/*! Assign and returns \c X coordinate to \a xx.*/
hint HPoint::X(hint xx)
{
    x = xx; return x;
}
/*! Assign and returns \c Y coordinate to \a yy.*/
hint HPoint::Y(hint yy)
{
    y = yy; return y;
}
/*! Returns difference value.*/
hfloat HPoint::DifferenceValue() const
{
    return difference;
}
/*!
  Compare \a p1 with \a p2 based on index.
*/
bool indexCmp(const HPoint &p1, const HPoint &p2)
{
    return p1.index < p2.index;
}
// .............................................................................
// Name: operator/
//
// Synopsis: a/n. Division by an hinteger
//
// Parameters:
//          hint n;  ----> second summand
//
// Returns:
//         a/n
//
// Exceptions:
//            None
//....................................................
/*!
  Division by the hint \a n.
*/
HPoint HPoint::operator/(hint n)
{
    HPoint newpohint(*this);
    newpohint/=n;
    return newpohint;
}
// .............................................................................
// Name: operator%
//
// Synopsis: a/n. Module operation by an hinteger to each coordinate
//
// Parameters:
//          hint n;  ----> second summand
//
// Returns:
//         a%n
//
// Exceptions:
//            None
//....................................................
/*!Module operation by \a n  to each coordinate*/
HPoint HPoint::operator%(hint n)
{
    HPoint newpohint(*this);
    newpohint%=n;
    return newpohint;
}
// .............................................................................
// Name: operator*
//
// Synopsis: a*n. Multiplies by an hinteger
//
// Parameters:
//          hint n;  ----> second summand
//
// Returns:
//         a*n
//
// Exceptions:
//            None
//....................................................
/*!Multiplies \a p by \a n.*/
HPoint operator*(const HPoint & p, hint n)
{
    HPoint pohint(p);
    pohint *= n;
    return pohint;
}
/*!Multiplies \a n by \a p.*/
HPoint operator*(hint n, const HPoint & p)
{
    return p*n;
}
// .............................................................................
// Name: logical operators
//
// Synopsis: 1) a == b  (a.y == b.y && a.x == b.x)
//           2) a != b  (a.y != b.y && a.x != b.x)
//           3) a < b   (a.y < b.y || (a.y == b.y && a.x < b.x))
//           4) a <= b  (! (a > b) )
//           5) a > b   (b < a)
//           6) a >= b  (! (a < b))
//
// Parameters:
//			const Pohint2D & a  ---> operand
//			const Pohint2D & b  ---> operand
//
// Returns:
//         bool value
//
// Exceptions:
//            None
//..............................................................................
/*!Compares \a p1 and \a p2.*/
bool operator==(const HPoint & p1, const HPoint & p2)
{
    return p1.x == p2.x && p1.y == p2.y;
}
//.........
/*!Compares \a p1 and \a p2.*/
bool operator!=(const HPoint & p1, const HPoint & p2)
{
    return !(p1==p2);
}
//........
/*!Compares \a p1 and \a p2.*/
bool operator>(const HPoint & p1, const HPoint & p2)
{
    bool flag=true;
    if(p1.y < p2.y)
        flag=false;
    else if(p1.y==p2.y && p1.x < p2.x)
        flag=false;
    return flag;
}
//.......
/*!Compares \a p1 and \a p2.*/
bool operator>=(const HPoint & p1, const HPoint & p2)
{
    return !(p2 > p1);
}
//.......
/*!Compares \a p1 and \a p2.*/
bool operator<(const HPoint & p1, const HPoint & p2)
{
    return p2>p1;
}
//......
/*!Compares \a p1 and \a p2.*/
bool operator<=(const HPoint & p1, const HPoint & p2)
{
    return !(p2 < p1);
}
// .............................................................................
// Name: operator<<
//
// Synopsis: dump operator
//
// Parameters:
//			std::ostream & os       ---> output stream
//			const Pohint & b  ---> object to be stored
//
// Returns:
//         the output stream
//
// Exceptions:
//            None
//..............................................................................
/*!Dump \a p to the output stream \a out.*/
std::ostream & operator<<(std::ostream & out, const HPoint & p)
{
    out << "(" << p.x << "," << p.y << ")";
    return out;
}
// .............................................................................
// Name: operator>>
//
// Synopsis: restore operator
//
// Parameters:
//			std::ostream & is       ---> input stream
//			const Pohint2D & b  ---> upon return restored object
//
// Returns:
//         the input stream
//
// Exceptions:
//            None
//..............................................................................
/*! Restore \a p from the input stream \a is*/
std::istream & operator>>(std::istream & is, HPoint & p)
{
    is >> p.x;
    is >> p.y;
    return is;
}
// .............................................................................
// Name: read
//
// Synopsis: read from ifstream
//
// Parameters:
//			std::ifstream & in       ---> input stream
//			const Point2D & b        ---> upon return restored object
//
// Returns:
//         the input stream
//
// Exceptions:
//            None
//..............................................................................
/*!Read \a p from an input stream \a in.*/
std::ifstream &read(std::ifstream &in, HPoint &p)
{
    hint x, y;
    in.read (reinterpret_cast<char *> (&x), sizeof(x));
    in.read (reinterpret_cast<char *> (&y), sizeof(y));
    p.x = x;
    p.y = y;
    return in;
}
// .............................................................................
// Name: write
//
// Synopsis: write to ofstream
//
// Parameters:
//			std::ofstream & out       ---> output stream
//			const Point2D & b        ---> Point to write
//
// Returns:
//         the input stream
//
// Exceptions:
//            None
//..............................................................................
/*! Write \a p to the output stream \a out*/
std::ofstream & write(std::ofstream &out, const HPoint &p)
{
    hint x = p.x; hint y =p.y;
    out.write(reinterpret_cast<char *> (&x), sizeof(x));
    out.write(reinterpret_cast<char *> (&y), sizeof(y));
    return out;
}
