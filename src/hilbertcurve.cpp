/*!
  \headerfile hilbertcurve.h
  \title Hilbert Curve Definitions
  \brief The \e{"hilbertcurve.h"} header file defines two classes that defines a HilbertCurve.
*/
#include "hilbertcurve.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <cmath>
#include <future>

#include "threads_utility.h"
#include "parallel_algorithm.h"


/*!
  \class QuasiSquare
  \since 0.1
  \inmodule hilbertlib
  \ingroup hcurve
  \brief Aproximated even partition of a square

     The QuasiSquare class is defined as in \e {"Approximately even partition algorithm for coding the
    Hilbert curve of arbitrary-sized image"} by \b {C.-C. Wu Y.-I. Chang}.
    It represents a \c {p x q} region where \c{|p-q| <= 1}

    The QuasiSquare class is used by HilbertCurve class for generating the corresponding
    curve.

    \sa HilbertCurve
*/
/*!
  \brief General constructor

  Constructs a QuasiSquare with dimensions \a nn x \a mm and orientation \a o starting
  at location given by \a point.
*/
QuasiSquare::QuasiSquare(hsize nn, hsize mm, HPoint point, Orientation o)
{
    n = nn;
    m = mm;
    coord = point;
    oABCD = o;
}
/*!

  \brief Default Constructor
*/
QuasiSquare::QuasiSquare(): n(0), m(0), coord(HPoint()), oABCD(A){}
/*!
  \brief Copy constructor of the class

  Constructs a copy of the given QuasiSquare \a q
*/
QuasiSquare::QuasiSquare(const QuasiSquare &q): n(q.n), m(q.m), coord(q.coord), oABCD(q.oABCD){}

//.....................................................................
//                      Operator overloading
//.....................................................................
/*!
 \fn  QuasiSquare & QuasiSquare::operator =(const QuasiSquare & qs)
 \brief Assign operator

 Copy \a qs and assign it.
 */
QuasiSquare & QuasiSquare::operator =(const QuasiSquare & qs)
{
    n = qs.n;
    m = qs.m;
    coord = qs.coord;
    oABCD = qs.oABCD;
    return *this;
}

/*!
    \brief Perform a QuasiSquare partition.
    Perform a even square partition. Partitioned QuasiSquare are returned as
    reference in \a partition_vec.
*/
std::vector<QuasiSquare> & QuasiSquare::Partition(std::vector<QuasiSquare> & partition_vec)
{
    QuasiSquare newQS;
    hsize n1, n2, m1, m2;
    HPoint point;

    //Even partition of the quasiquare
    n1 = n/2;
    n2 = n-n1;
    m1 = m/2;
    m2 = m - m1;

    if(oABCD == A || oABCD == B)
    {
        if(n1%2 == 1) std::swap(n1, n2);
        if(m1%2 == 1) std::swap(m1 , m2);
    }
    if(oABCD == C || oABCD == D)
    {
        if(n2%2 == 1) std::swap(n1, n2);
        if(m2%2 == 1) std::swap(m1 , m2);
    }

    try
    {
        switch (oABCD)
        {
            case A:
                point.X (coord.X()+m1);
                point.Y (coord.Y());

                newQS = QuasiSquare(n1, m2, point, D);
                partition_vec.push_back (newQS);

                point.X(coord.X()+m1);
                point.Y(coord.Y()+n1);

                newQS = QuasiSquare(n2, m2, point, A);
                partition_vec.push_back (newQS);

                point.X (coord.X());
                point.Y (coord.Y() + n1);

                newQS = QuasiSquare(n2, m1, point, A);
                partition_vec.push_back (newQS);

                point.X (coord.X());
                point.Y (coord.Y());

                newQS = QuasiSquare(n1, m1, point, B);
                partition_vec.push_back (newQS);

                break;

            case B:

                point.X (coord.X());
                point.Y (coord.Y() + n1);

                newQS = QuasiSquare(n2, m1, point, C);
                partition_vec.push_back (newQS);

                point.X(coord.X()+m1);
                point.Y(coord.Y()+n1);

                newQS = QuasiSquare(n2, m2, point, B);
                partition_vec.push_back (newQS);

                point.X (coord.X()+m1);
                point.Y (coord.Y());

                newQS = QuasiSquare(n1, m2, point, B);
                partition_vec.push_back (newQS);

                point.X (coord.X());
                point.Y (coord.Y());

                newQS = QuasiSquare(n1, m1, point, A);
                partition_vec.push_back (newQS);

                break;

            case C:

                point.X (coord.X());
                point.Y (coord.Y() + n1);

                newQS = QuasiSquare(n2, m1, point, B);
                partition_vec.push_back (newQS);

                point.X (coord.X());
                point.Y (coord.Y());

                newQS = QuasiSquare(n1, m1, point, C);
                partition_vec.push_back (newQS);

                point.X (coord.X()+m1);
                point.Y (coord.Y());

                newQS = QuasiSquare(n1, m2, point, C);
                partition_vec.push_back (newQS);

                point.X(coord.X()+m1);
                point.Y(coord.Y()+n1);

                newQS = QuasiSquare(n2, m2, point, D);
                partition_vec.push_back (newQS);

                break;

            case D:

                point.X (coord.X()+m1);
                point.Y (coord.Y());

                newQS = QuasiSquare(n1, m2, point, A);
                partition_vec.push_back (newQS);

                point.X (coord.X());
                point.Y (coord.Y());

                newQS = QuasiSquare(n1, m1, point, D);
                partition_vec.push_back (newQS);

                point.X (coord.X());
                point.Y (coord.Y() + n1);

                newQS = QuasiSquare(n2, m1, point, D);
                partition_vec.push_back (newQS);

                point.X(coord.X()+m1);
                point.Y(coord.Y()+n1);

                newQS = QuasiSquare(n2, m2, point, C);
                partition_vec.push_back (newQS);

                break;
        }
    }
    catch (std::exception &e)
    {
        throw e;
    }
    return partition_vec;
}
/*!
  \brief Build the Hilbert Curve recursively.
  Build the Hilbert curve returning the coordinates as reference in \a coordinates_list.
*/
std::vector<HPoint> &QuasiSquare::BuildCurve(std::vector<HPoint> & coordinates_list, hsize index)
{

    std::vector<QuasiSquare> qsv;
    QuasiSquare qs;
    HPoint p;

    if(n > 2 || m > 2)//QuasiSquare isn't a primitive so need to keep Partitioning
    {
        Partition (qsv);
        //std::future<void> futures[2];
        for(int i=0; i < 4; ++i)
        {
            qs = qsv.back ();
            qsv.pop_back ();
            // Proccess only the half in other threads. This will be happening recursively
            if(i < 2)
            {
                // Wrapping the BuildCurve method in a function object
                std::function<void()> func = std::bind(&QuasiSquare::BuildCurve, qs, std::ref(coordinates_list), index);
                // Push the function to the thread pool
                thread_pool::instance ().push_task(func);
            }
            else
            {
                qs.BuildCurve (coordinates_list, index);
            }
            //qs.BuildCurve (coordinates_list, index);
            index += qs.n * qs.m;
        }
        //futures[0].wait (); futures[1].wait ();
    }
//    std::vector<QuasiSquare> qsv;
//    QuasiSquare qs;
//    HPoint p;

//    if(n > 2 || m > 2)
//    {
//        Partition (qsv);
//        std::vector<HPoint> coordinates[4];
//        std::future<std::vector<HPoint>&> future_coordinates[4];
//        for(int i=0; i < 4; ++i)
//        {
//            qs = qsv.back ();
//            qsv.pop_back ();

//            future_coordinates[i] = std::async(&QuasiSquare::BuildCurve, std::ref(qs), std::ref(coordinates[i]));
//            //qs.BuildCurve (coordinates_list);
//        }
//        for(int i=0; i < 4; ++i)
//        {
//            future_coordinates[i].wait ();
//            coordinates_list.insert (coordinates_list.end (), coordinates[i].begin (), coordinates[i].end ());
//        }
//    }
//    std::cout << "BuildCurve. Current Thread: " << std::this_thread::get_id () << std::endl;
//    std::cout << "Start index: " << index << std::endl;
//    std::cout << "m: " << m << " n: " << n << " coord: " << coord << std::endl;
    try
    {
        if(n == 1 && m == 1)
        {
            coordinates_list[index] = coord;
//            coordinates_list[index].index = index;
        }
        if(n == 1 && m == 2)
        {
            if(oABCD == A|| oABCD == B)
            {
                coordinates_list[index] = coord;
//                coordinates_list[index].index = index;
                p.X(coord.X() + 1);
                p.Y(coord.Y());
//                p.index = index+1;
                coordinates_list[index+1] = p;
            }
            else if(oABCD == C || oABCD == D)
            {
                p.X(coord.X() + 1);
                p.Y(coord.Y());
                coordinates_list[index] = p;
//                coordinates_list[index].index = index;
                coordinates_list[index+1] = coord;
//                coordinates_list[index+1].index = index+1;
            }
        }
        if(n == 2 && m == 1)
        {
            if(oABCD == B || oABCD == A)
            {
                coordinates_list[index] = coord;
//                coordinates_list[index].index = index;
                p.X(coord.X());
                p.Y(coord.Y()+1);
                coordinates_list[index+1] = p;
//                coordinates_list[index+1].index = index+1;

            }
            else if(oABCD == D || oABCD == C)
            {
                p.X(coord.X());
                p.Y(coord.Y()+1);
                coordinates_list[index] = p;
//                coordinates_list[index].index = index;
                coordinates_list[index+1] = coord;
//                coordinates_list[index+1].index = index+1;
            }
        }
        if(n == 2 && m == 2)
        {
            if(oABCD == A)
            {
                coordinates_list[index] = coord;
//                coordinates_list[index].index = index;
                p.X(coord.X());
                p.Y(coord.Y()+1);
                coordinates_list[index+1] = p;
//                coordinates_list[index+1].index = index+1;
                p.X(coord.X()+1);
                coordinates_list[index+2] = p;
//                coordinates_list[index+2].index = index+2;
                p.Y(p.Y()-1);
                coordinates_list[index+3] = p;
//                coordinates_list[index+3].index = index+3;
            }
            else if(oABCD == B)
            {
                coordinates_list[index] = coord;
//                coordinates_list[index].index = index;
                p.X(coord.X()+1);
                p.Y(coord.Y());
                coordinates_list[index+1] = p;
//                coordinates_list[index+1].index = index+1;
                p.Y(coord.Y()+1);
                coordinates_list[index+2] = p;
//                coordinates_list[index+2].index = index+2;
                p.X(p.X()-1);
                coordinates_list[index+3] = p;
//                coordinates_list[index+3].index = index+3;
            }
            else if(oABCD == C)
            {
                p.X(coord.X()+1);
                p.Y(coord.Y()+1);
                coordinates_list[index] = p;
//                coordinates_list[index].index = index;
                p.Y(coord.Y());
                coordinates_list[index+1] = p;
//                coordinates_list[index+1].index = index+1;
                coordinates_list[index+2] = coord;
//                coordinates_list[index+2].index = index+2;
                p.X(coord.X());
                p.Y(coord.Y()+1);
                coordinates_list[index+3] = p;
//                coordinates_list[index+3].index = index+3;

            }
            else if(oABCD == D)
            {
                p.X(coord.X()+1);
                p.Y(coord.Y()+1);
                coordinates_list[index] = p;
//                coordinates_list[index].index = index;
                p.X(coord.X());
                coordinates_list[index+1] = p;
//                coordinates_list[index+1].index = index+1;
                coordinates_list[index+2] = coord;
//                coordinates_list[index+2].index = index+2;
                p.X(coord.X()+1);
                p.Y(coord.Y());
                coordinates_list[index+3] = p;
//                coordinates_list[index+3].index = index+3;
            }
        }

    }
    catch (std::bad_alloc& ba)
    {
        throw HilbertBadAlloc();
    }
//    std::cout << "Printing curve for debuggin" << std::endl;
//    for(auto val : coordinates_list)
//    {
//        std::cout << val << " ";
//    }
//    std::cout << std::endl;
    return coordinates_list;
}
/*!
  \brief Returns the curve lenghts
*/
hsize HilbertCurve::lenght() const
{
    return hsize(m_curve.size ());
}
/*!
    \brief Returns the plot width
*/
hsize HilbertCurve::width() const
{
    return this->m;
}
/*!
 \brief Returns the plot height.
*/
hsize HilbertCurve::height() const
{
    return this->n;
}
/*!
    \brief Returns the curve type.
*/
HilbertCurve::CurveType HilbertCurve::type() const
{
    return m_type;
}
/*!
    \brief Reference operator
    Return a reference to the HPoint store at \a index.

*/
std::vector<HPoint>::reference HilbertCurve::operator[](std::vector<HPoint>::size_type index)
{
    return m_curve.at (index);
}
/*!
  \overload operator[]()

  Returns a constant reference to the HPoint at \a index.
*/
std::vector<HPoint>::const_reference HilbertCurve::operator[](std::vector<HPoint>::size_type index) const
{
    return m_curve.at (index);
}
/*!
    \brief Begin Iterator
    Returns an iterator at the begin of the vector.
*/
HilbertCurve::const_iterator HilbertCurve::begin() const
{
    return m_curve.begin ();
}
/*!
    \brief Const Begin Iterator
    Returns a constant iterator at the begin of the vector.
*/
HilbertCurve::const_iterator HilbertCurve::cbegin() const
{
    return m_curve.cbegin ();
}
/*!
    \brief End Iterator
    Returns an iterator at the end of the vector.
*/
HilbertCurve::const_iterator HilbertCurve::end() const
{
    return m_curve.end ();
}
/*!
    \brief Const End Iterator
    Returns a constant iterator at the end of the vector.
*/
HilbertCurve::const_iterator HilbertCurve::cend() const
{
    return m_curve.cend ();
}
/*!
  \fn void HilbertCurve::SaveSVG(const char *filename, const char *colorName, float stroke_width)
  \brief Export curve as an SVG file

  Export the HilbertCurve using SVG format to \a filename, using \a colorName and \a stroke_width
  to stylish the output curve.
  \note \a colorName must be a valid HTML color.
  \note \a stroke_width must be in range [0,1]
*/
void HilbertCurve::SaveSVG(const char *filename, const char *colorName, float stroke_width)
{
    reflectY ();
    std::ofstream os(filename);
    std::vector<HPoint> coordinates = m_curve;
    reflectY ();
        unsigned int xmax, ymax;

        xmax=std::max_element(coordinates.begin(), coordinates.end(), [](HPoint & p1,HPoint & p2) { return (p1.X() < p2.X()); })->X();
        ymax=std::max_element(coordinates.begin(), coordinates.end(), [](HPoint & p1,HPoint & p2) { return (p1.Y() < p2.Y()); })->Y();


        os << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << std::endl;

        os << "<svg" << std::endl;
        os << "width=\"" << xmax << "\"" << std::endl;
        os << "height=\"" << ymax << "\"" << std::endl;
        os << "id=\"svg2\"" << std::endl;
        os << "version=\"1.1\">" << std::endl;
        os << "<g>" << std::endl;
        os << "<path" << std::endl;
        os << "style=\"fill:none;stroke:"<< colorName <<";stroke-width:" << stroke_width <<"px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"" << std::endl;
        os << "d=\"M ";
        for(auto instance : coordinates)
            os << instance.X() << "," << instance.Y() << " ";
        os << "\"/>" << std::endl;
        os << "</g>" << std::endl;
        os << "</svg>";

        os.close ();
}
/*!
  \brief Generate a string with the curve representation as SVG.
  Returns a std::string representing the HilbertCurve using SVG format,
  \a colorName and \a stroke_width
  are used to stylish the output.
  \note \a colorName must be a valid HTML color.
  \note \a stroke_width must be in range [0,1]
*/
std::string HilbertCurve::curveToSvg(const char *colorName, float stroke_width)
{
    reflectY ();
    std::vector<HPoint> coordinates = m_curve;
    reflectY ();
    hsize xmax, ymax;

    xmax=std::max_element(coordinates.begin(), coordinates.end(), [](HPoint & p1,HPoint & p2) { return (p1.X() < p2.X()); })->X();
    ymax=std::max_element(coordinates.begin(), coordinates.end(), [](HPoint & p1,HPoint & p2) { return (p1.Y() < p2.Y()); })->Y();

    std::ostringstream os;
    os << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << std::endl;

    os << "<svg" << std::endl;
    os << "width=\"" << xmax << "\"" << std::endl;
    os << "height=\"" << ymax << "\"" << std::endl;
    os << "id=\"svg2\"" << std::endl;
    os << "version=\"1.1\">" << std::endl;
    os << "<g>" << std::endl;
    os << "<path" << std::endl;
    os << "style=\"fill:none;stroke:"<< colorName <<";stroke-width:" << stroke_width <<"px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1\"" << std::endl;
    os << "d=\"M ";
    for(auto instance : coordinates)
        os << instance.X() << "," << instance.Y() << " ";
    os << "\"/>" << std::endl;
    os << "</g>" << std::endl;
    os << "</svg>";

    return os.str ();
}

HilbertCurve HilbertCurve::createCurve(hsize width, hsize height, HilbertCurve::CurveType type, HPoint origen, QuasiSquare::Orientation orientation, bool differenceCurve)
{
    return HilbertCurve(width, height, type, origen, orientation, differenceCurve);
}

void HilbertCurve::Build()
{
    switch (m_type) {
        case H0:
            BuildCurveH0 ();
            break;
        case H1:
            BuildCurve1H ();
            break;
        case H2:
            BuildCurve2H ();
            break;
        case H3:
            BuildCurve3H ();
            break;
        case H4:
            BuildCurve4H ();
            break;
        case H5:
            BuildCurve5H ();
            break;
        case H6:
            BuildCurve6H ();
            break;
        case H7:
            BuildCurve7H ();
            break;
        case H8:
            BuildCurve8H ();
            break;
        case H9:
            BuildCurve9H ();
            break;
        case H10:
            BuildCurve10H ();
            break;
        case H11:
            BuildCurve11H ();
            break;
        case H12:
            BuildCurve12H ();
            break;
        case H13:
            BuildCurve13H ();
            break;
        case H14:
            BuildCurve14H ();
            break;
        case H15:
            BuildCurve15H ();
            break;
        case H16:
            BuildCurve16H ();
            break;
        case H17:
            BuildCurve17H ();
            break;
        case H18:
            BuildCurve18H ();
            break;
        case H19:
            BuildCurve19H ();
            break;
        case H20:
            BuildCurve20H ();
            break;
        case H21:
            BuildCurve21H ();
            break;
        case H22:
            BuildCurve22H ();
            break;
        case H23:
            BuildCurve23H ();
            break;
        case H24:
            BuildCurve24H ();
            break;
        case H25:
            BuildCurve25H ();
            break;
        case H26:
            BuildCurve26H ();
            break;
        case H27:
            BuildCurve27H ();
            break;
        case H28:
            BuildCurve28H ();
            break;
        case H29:
            BuildCurve29H ();
            break;
        case H30:
            BuildCurve30H ();
            break;
        case H31:
            BuildCurve31H ();
            break;
        case H32:
            BuildCurve32H ();
            break;
        case H33:
            BuildCurve33H ();
            break;
        case H34:
            BuildCurve34H ();
            break;
        case H35:
            BuildCurve35H ();
            break;
        case H36:
            BuildCurve36H ();
            break;
        case H37:
            BuildCurve37H ();
            break;
        case H38:
            BuildCurve38H ();
            break;
        case H39:
            BuildCurve39H ();
            break;

    }
}
/*!
  Assign index and difference values to HPoints.
*/
void HilbertCurve::BuildDifference()
{
    uint ind =0;
    for (auto &v : m_curve)
    {
        v.index = ind++;
    }

    // Calculate mean_difference

    // ## Dunno why on Android std::sort dont work as expected
#ifdef Q_OS_ANDROID
    qSort(m_curve.begin (), m_curve.end ());
#else
    std::sort (m_curve.begin (), m_curve.end ());
#endif

    std::vector<HPoint>::size_type i=0;
    std::vector<HPoint>::size_type j=0;
    std::vector<HPoint>::size_type k=0;
    hfloat dif=0;
    hfloat p=0;
    hfloat mean=0;
    hfloat delta=0;
    hfloat val=0;
    hsize width = this->width ();
    hsize height = this->height ();

    for(j=0; j< height; ++j)
    {
        for(i=0; i< width; ++i)
        {
            int count = 0;
            dif = 0;
            p=m_curve[j*width+i].index;
            // Right
            if(i < (width - 1))
            {
                dif += std::fabs(p-m_curve[j*width+i+1].index);
                ++count;
            }
            // Left
            if(i > 0)
            {
                dif += std::fabs(p-m_curve[j*width+i-1].index);
                ++count;
            }
            // Up
            if(j < (height - 1))
            {
                dif += std::fabs(p-m_curve[(j+1)*width+i].index);
                ++count;
            }
            // Down
            if(j > 0)
            {
                dif += std::fabs(p-m_curve[(j-1)*width+i].index);
                ++count;
            }
            // Right Up
            if((i < (width - 1)) && (j < (height - 1)))
            {
                 dif += std::fabs(p-m_curve[(j+1)*width+i+1].index);
                ++count;
            }
            // Left Up
            if((i > 0) && (j < (height-1)))
            {
                dif += std::fabs(p-m_curve[(j+1)*width+i-1].index);
                ++count;
            }
            // Right Down
            if((i < (width - 1)) && (j > 0))
            {
                dif += std::fabs(p-m_curve[(j-1)*width+i+1].index);
                ++count;
            }
            // Left Down
            if((i > 0) && (j > 0))
            {
                dif += std::fabs(p-m_curve[(j-1)*width+i-1].index);
                ++count;
            }
            m_curve[j*width+i].difference=val=dif/count;

            delta=val-mean;
            mean=mean+delta/++k;
        }
    }
    m_mean_difference=mean;
    std::sort(m_curve.begin (), m_curve.end (), indexCmp);
}

void HilbertCurve::BuildCurveH0()
{
    m_curve.assign (n * m, 0);
    BuildCurve(m_curve, 0);
    // Ok esto funciona bien y es bastante optimo, el riesgo esta en que otro thread
    // pushee tareas en thread_pool::instance. En ese caso este proceso va a tomar tiempo
    // en completar esas tareas antes de retornar. Esto no es muy OP. En cualquier caso
    // si el resto del codigo cliente no usa threads para otras tareas esto es seguro.
    // Mas adelante solucionare esto con un thread_pool dentro de la clase y una bandera de
    // compilacion indicando si se usa o no.
    while (thread_pool::instance ().isWorking ())
    {
        thread_pool::instance ().run_task ();
    }

}
void HilbertCurve::BuildCurve1H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    switch (oABCD)
    {
        case A:
        case C:
            o1 = D;
            o2 = D;
            o3 = B;
            o4 = B;
            break;
        case B:
        case D:
            o1 = C;
            o2 = A;
            o3 = A;
            o4 = C;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, H0, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, H0, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, H0, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, H0, coord + HPoint(w1, 0), o4);


    //Reversing curves for orientation
    cuadrant1.reverse ();
    cuadrant2.reverse ();
    cuadrant3.reverse ();
    cuadrant4.reverse ();
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);

}
void HilbertCurve::BuildCurve2H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    switch (oABCD)
    {
        case A:
        case C:
            o1 = C;
            o2 = A;
            o3 = A;
            o4 = C;
            break;
        case B:
        case D:
            o1 = D;
            o2 = D;
            o3 = B;
            o4 = B;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, H0, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, H0, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, H0, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, H0, coord + HPoint(w1, 0), o4);
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve3H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    switch (oABCD)
    {
        case A:
            o1 = C;
            o2 = D;
            o3 = B;
            o4 = C;
            break;
        case B:
            o1 = D;
            o2 = D;
            o3 = A;
            o4 = C;
            break;
        case C:
            o1 = D;
            o2 = A;
            o3 = A;
            o4 = B;
            break;
        case D:
            o1 = C;
            o2 = A;
            o3 = B;
            o4 = B;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, H0, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, H0, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, H0, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, H0, coord + HPoint(w1, 0), o4);

    //Reversing curves for orientation
    cuadrant1.reverse ();
    cuadrant2.reverse ();
    cuadrant3.reverse ();
    cuadrant4.reverse ();
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve4H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    switch (oABCD)
    {
        case A:
            o1 = B;
            o2 = A;
            o3 = A;
            o4 = C;
            break;
        case B:
            o1 = A;
            o2 = D;
            o3 = B;
            o4 = B;
            break;
        case C:
            o1 = C;
            o2 = A;
            o3 = D;
            o4 = C;
            break;
        case D:
            o1 = D;
            o2 = D;
            o3 = C;
            o4 = B;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, H0, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, H0, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, H0, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, H0, coord + HPoint(w1, 0), o4);
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve5H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    switch (oABCD)
    {
        case A:
            o1 = C;
            o2 = D;
            o3 = B;
            o4 = B;
            break;
        case B:
            o1 = D;
            o2 = A;
            o3 = A;
            o4 = C;
            break;
        case C:
            o1 = D;
            o2 = D;
            o3 = A;
            o4 = B;
            break;
        case D:
            o1 = C;
            o2 = A;
            o3 = B;
            o4 = C;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, H0, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, H0, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, H0, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, H0, coord + HPoint(w1, 0), o4);

    //Reversing curves for orientation
    cuadrant1.reverse ();
    cuadrant2.reverse ();
    cuadrant3.reverse ();
    cuadrant4.reverse ();
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve6H ()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    switch (oABCD)
    {
        case A:
        case C:
            o1 = C;
            o2 = A;
            o3 = A;
            o4 = C;
            break;
        case B:
        case D:
            o1 = D;
            o2 = D;
            o3 = B;
            o4 = B;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, H5, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, H5, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, H5, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, H5, coord + HPoint(w1, 0), o4);

    cuadrant2.reflectAndReverse ();
    cuadrant4.reflectAndReverse ();
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve7H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    switch (oABCD)
    {
        case A:
            o1 = C;
            o2 = A;
            o3 = A;
            o4 = D;
            break;
        case B:
            o1 = A;
            o2 = D;
            o3 = B;
            o4 = B;
            break;
        case C:
            o1 = C;
            o2 = B;
            o3 = A;
            o4 = C;
            break;
        case D:
            o1 = D;
            o2 = D;
            o3 = C;
            o4 = B;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, H5, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, H5, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, H5, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, H5, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            cuadrant2.reflectAndReverse ();
            break;
        case B:
            cuadrant1.reflectAndReverse ();
            cuadrant2.reflectAndReverse ();
            cuadrant4.reflectAndReverse ();
            break;
        case C:
            cuadrant4.reflectAndReverse ();
            break;
        case D:
            cuadrant2.reflectAndReverse ();
            cuadrant3.reflectAndReverse ();
            cuadrant4.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve8H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    switch (oABCD)
    {
        case A:
            o1 = B;
            o2 = A;
            o3 = A;
            o4 = D;
            break;
        case B:
            o1 = A;
            o2 = C;
            o3 = B;
            o4 = B;
            break;
        case C:
            o1 = C;
            o2 = B;
            o3 = D;
            o4 = C;
            break;
        case D:
            o1 = D;
            o2 = D;
            o3 = C;
            o4 = A;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, H5, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, H5, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, H5, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, H5, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            cuadrant1.reflectAndReverse ();
            cuadrant2.reflectAndReverse ();
            break;
        case B:
            cuadrant1.reflectAndReverse ();
            cuadrant4.reflectAndReverse ();
            break;
        case C:
            cuadrant3.reflectAndReverse ();
            cuadrant4.reflectAndReverse ();
            break;
        case D:
            cuadrant2.reflectAndReverse ();
            cuadrant3.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve9H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    switch (oABCD)
    {
        case A:
        case C:
            o1 = D;
            o2 = D;
            o3 = B;
            o4 = B;
            break;
        case B:
        case D:
            o1 = C;
            o2 = A;
            o3 = A;
            o4 = C;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, H5, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, H5, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, H5, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, H5, coord + HPoint(w1, 0), o4);

    cuadrant1.reverse ();
    cuadrant3.reverse ();
    cuadrant2.reflect ();
    cuadrant4.reflect ();

    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve10H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    switch (oABCD)
    {
        case A:
            o1 = C;
            o2 = D;
            o3 = B;
            o4 = C;
            break;
        case B:
            o1 = D;
            o2 = D;
            o3 = A;
            o4 = C;
            break;
        case C:
            o1 = D;
            o2 = A;
            o3 = A;
            o4 = B;
            break;
        case D:
            o1 = C;
            o2 = A;
            o3 = B;
            o4 = B;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, H5, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, H5, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, H5, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, H5, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            cuadrant1.reflect ();
            cuadrant2.reflect ();
            cuadrant3.reverse ();
            cuadrant4.reverse ();
            break;
        case B:
            cuadrant1.reflect ();
            cuadrant2.reverse ();
            cuadrant3.reverse ();
            cuadrant4.reflect ();
            break;
        case C:
            cuadrant1.reverse ();
            cuadrant2.reverse ();
            cuadrant3.reflect ();
            cuadrant4.reflect ();
            break;
        case D:
            cuadrant1.reverse ();
            cuadrant2.reflect ();
            cuadrant3.reflect ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve11H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    switch (oABCD)
    {
        case A:
            o1 = C;
            o2 = D;
            o3 = B;
            o4 = B;
            break;
        case B:
            o1 = C;
            o2 = D;
            o3 = A;
            o4 = C;
            break;
        case C:
            o1 = D;
            o2 = D;
            o3 = A;
            o4 = B;
            break;
        case D:
            o1 = C;
            o2 = A;
            o3 = A;
            o4 = B;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, H5, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, H5, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, H5, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, H5, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            cuadrant1.reflect ();
            cuadrant2.reflect ();
            cuadrant3.reverse ();
            cuadrant4.reflect ();
            break;
        case B:
            cuadrant1.reverse ();
            cuadrant2.reverse ();
            cuadrant3.reverse ();
            cuadrant4.reflect ();
            break;
        case C:
            cuadrant1.reverse ();
            cuadrant2.reflect ();
            cuadrant3.reflect ();
            cuadrant4.reflect ();
            break;
        case D:
            cuadrant1.reverse ();
            cuadrant2.reflect ();
            cuadrant3.reverse ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve12H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = B; t1 = H3;
            o2 = A; t2 = H5;
            o3 = A; t3 = H5;
            o4 = D; t4 = H3;
            break;
        case B:
            o1 = A; t1 = H3;
            o2 = C; t2 = H3;
            o3 = B; t3 = H5;
            o4 = B; t4 = H5;
            break;
        case C:
            o1 = C; t1 = H5;
            o2 = B; t2 = H3;
            o3 = D; t3 = H3;
            o4 = C; t4 = H5;
            break;
        case D:
            o1 = D; t1 = H5;
            o2 = D; t2 = H5;
            o3 = C; t3 = H3;
            o4 = A; t4 = H3;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
        case D:
            cuadrant2.reflectAndReverse ();
            break;
        case B:
        case C:
            cuadrant4.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve13H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = D; t1 = H3;
            o2 = D; t2 = H5;
            o3 = B; t3 = H5;
            o4 = B; t4 = H3;
            break;
        case B:
            o1 = C; t1 = H3;
            o2 = A; t2 = H3;
            o3 = A; t3 = H5;
            o4 = C; t4 = H5;
            break;
        case C:
            o1 = D; t1 = H5;
            o2 = D; t2 = H3;
            o3 = B; t3 = H3;
            o4 = B; t4 = H5;
            break;
        case D:
            o1 = C; t1 = H5;
            o2 = A; t2 = H5;
            o3 = A; t3 = H3;
            o4 = C; t4 = H3;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    cuadrant1.reverse ();
    cuadrant3.reverse ();
    switch (oABCD)
    {
        case A:
        case D:
            cuadrant2.reflect ();
            cuadrant4.reverse ();
            break;
        case B:
        case C:
            cuadrant2.reverse ();
            cuadrant4.reflect ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve14H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = B; t1 = H3;
            o2 = A; t2 = H5;
            o3 = A; t3 = H5;
            o4 = D; t4 = H5;
            break;
        case B:
            o1 = A; t1 = H5;
            o2 = C; t2 = H3;
            o3 = B; t3 = H5;
            o4 = B; t4 = H5;
            break;
        case C:
            o1 = C; t1 = H5;
            o2 = B; t2 = H5;
            o3 = D; t3 = H3;
            o4 = C; t4 = H5;
            break;
        case D:
            o1 = D; t1 = H5;
            o2 = D; t2 = H5;
            o3 = C; t3 = H5;
            o4 = A; t4 = H3;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            cuadrant2.reflectAndReverse ();
            break;
        case B:
            cuadrant1.reflectAndReverse ();
            cuadrant4.reflectAndReverse ();
            break;
        case C:
            cuadrant4.reflectAndReverse ();
            break;
        case D:
            cuadrant2.reflectAndReverse ();
            cuadrant3.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve15H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = D; t1 = H3;
            o2 = D; t2 = H5;
            o3 = B; t3 = H5;
            o4 = C; t4 = H5;
            break;
        case B:
            o1 = D; t1 = H5;
            o2 = A; t2 = H3;
            o3 = A; t3 = H5;
            o4 = C; t4 = H5;
            break;
        case C:
            o1 = D; t1 = H5;
            o2 = A; t2 = H5;
            o3 = B; t3 = H3;
            o4 = B; t4 = H5;
            break;
        case D:
            o1 = C; t1 = H5;
            o2 = A; t2 = H5;
            o3 = B; t3 = H5;
            o4 = C; t4 = H3;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            cuadrant1.reverse ();
            cuadrant2.reflect ();
            cuadrant3.reverse ();
            cuadrant4.reverse ();
            break;
        case B:
            cuadrant1.reflect ();
            cuadrant2.reverse ();
            cuadrant3.reverse ();
            cuadrant4.reflect ();
            break;
        case C:
            cuadrant1.reverse ();
            cuadrant2.reverse ();
            cuadrant3.reverse ();
            cuadrant4.reflect ();
            break;
        case D:
            cuadrant1.reverse ();
            cuadrant2.reflect ();
            cuadrant3.reflect ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve16H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = B; t1 = H3;
            o2 = A; t2 = H5;
            o3 = A; t3 = H5;
            o4 = C; t4 = H5;
            break;
        case B:
            o1 = D; t1 = H5;
            o2 = C; t2 = H3;
            o3 = B; t3 = H5;
            o4 = B; t4 = H5;
            break;
        case C:
            o1 = C; t1 = H5;
            o2 = A; t2 = H5;
            o3 = D; t3 = H3;
            o4 = C; t4 = H5;
            break;
        case D:
            o1 = D; t1 = H5;
            o2 = D; t2 = H5;
            o3 = B; t3 = H5;
            o4 = A; t4 = H3;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            cuadrant2.reflectAndReverse ();
            cuadrant4.reflectAndReverse ();
            break;
        case B:
            cuadrant4.reflectAndReverse ();
            break;
        case C:
            cuadrant2.reflectAndReverse ();
            cuadrant4.reflectAndReverse ();
            break;
        case D:
            cuadrant2.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve17H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = D; t1 = H3;
            o2 = D; t2 = H5;
            o3 = B; t3 = H5;
            o4 = B; t4 = H5;
            break;
        case B:
            o1 = C; t1 = H5;
            o2 = A; t2 = H3;
            o3 = A; t3 = H5;
            o4 = C; t4 = H5;
            break;
        case C:
            o1 = D; t1 = H5;
            o2 = D; t2 = H5;
            o3 = B; t3 = H3;
            o4 = B; t4 = H5;
            break;
        case D:
            o1 = C; t1 = H5;
            o2 = A; t2 = H5;
            o3 = A; t3 = H5;
            o4 = C; t4 = H3;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            cuadrant1.reverse ();
            cuadrant2.reflect ();
            cuadrant3.reverse ();
            cuadrant4.reflect ();
            break;
        case B:
            cuadrant1.reverse ();
            cuadrant2.reverse ();
            cuadrant3.reverse ();
            cuadrant4.reflect ();
            break;
        case C:
            cuadrant1.reverse ();
            cuadrant2.reflect ();
            cuadrant3.reverse ();
            cuadrant4.reflect ();
            break;
        case D:
            cuadrant1.reverse ();
            cuadrant2.reflect ();
            cuadrant3.reverse ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);

}
void HilbertCurve::BuildCurve18H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = B; t1 = H4;
            o2 = A; t2 = H0;
            o3 = A; t3 = H0;
            o4 = D; t4 = H4;
            break;
        case B:
            o1 = A; t1 = H4;
            o2 = C; t2 = H4;
            o3 = B; t3 = H0;
            o4 = B; t4 = H0;
            break;
        case C:
            o1 = C; t1 = H0;
            o2 = B; t2 = H4;
            o3 = D; t3 = H4;
            o4 = C; t4 = H0;
            break;
        case D:
            o1 = D; t1 = H0;
            o2 = D; t2 = H0;
            o3 = C; t3 = H4;
            o4 = A; t4 = H4;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
        case B:
            cuadrant1.reflectAndReverse ();
            break;
        case C:
        case D:
            cuadrant3.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve19H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = C; t1 = H4;
            o2 = A; t2 = H0;
            o3 = A; t3 = H0;
            o4 = C; t4 = H4;
            break;
        case B:
            o1 = D; t1 = H4;
            o2 = D; t2 = H4;
            o3 = B; t3 = H0;
            o4 = B; t4 = H0;
            break;
        case C:
            o1 = C; t1 = H0;
            o2 = A; t2 = H4;
            o3 = A; t3 = H4;
            o4 = C; t4 = H0;
            break;
        case D:
            o1 = D; t1 = H0;
            o2 = D; t2 = H0;
            o3 = B; t3 = H4;
            o4 = B; t4 = H4;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
        case B:
            cuadrant1.reflectAndReverse ();
            break;
        case C:
        case D:
            cuadrant3.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve20H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = B; t1 = H4;
            o2 = A; t2 = H0;
            o3 = A; t3 = H0;
            o4 = C; t4 = H4;
            break;
        case B:
            o1 = D; t1 = H4;
            o2 = C; t2 = H4;
            o3 = B; t3 = H0;
            o4 = B; t4 = H0;
            break;
        case C:
            o1 = C; t1 = H0;
            o2 = A; t2 = H4;
            o3 = D; t3 = H4;
            o4 = C; t4 = H0;
            break;
        case D:
            o1 = D; t1 = H0;
            o2 = D; t2 = H0;
            o3 = B; t3 = H4;
            o4 = A; t4 = H4;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
        case B:
            cuadrant1.reflectAndReverse ();
            break;
        case C:
        case D:
            cuadrant3.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve21H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = C; t1 = H4;
            o2 = D; t2 = H0;
            o3 = B; t3 = H0;
            o4 = C; t4 = H4;
            break;
        case B:
            o1 = D; t1 = H4;
            o2 = D; t2 = H4;
            o3 = A; t3 = H0;
            o4 = C; t4 = H0;
            break;
        case C:
            o1 = D; t1 = H0;
            o2 = A; t2 = H4;
            o3 = A; t3 = H4;
            o4 = B; t4 = H0;
            break;
        case D:
            o1 = C; t1 = H0;
            o2 = A; t2 = H0;
            o3 = B; t3 = H4;
            o4 = B; t4 = H4;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    cuadrant1.reverse ();
    cuadrant3.reverse ();
    switch (oABCD)
    {
        case A:
        case D:
            cuadrant2.reverse ();
            cuadrant4.reflect ();
            break;
        case C:
        case B:
            cuadrant2.reflect ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve22H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = D; t1 = H4;
            o2 = D; t2 = H0;
            o3 = B; t3 = H0;
            o4 = B; t4 = H4;
            break;
        case B:
            o1 = C; t1 = H4;
            o2 = A; t2 = H4;
            o3 = A; t3 = H0;
            o4 = C; t4 = H0;
            break;
        case C:
            o1 = D; t1 = H0;
            o2 = D; t2 = H4;
            o3 = B; t3 = H4;
            o4 = B; t4 = H0;
            break;
        case D:
            o1 = C; t1 = H0;
            o2 = A; t2 = H0;
            o3 = A; t3 = H4;
            o4 = C; t4 = H4;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    cuadrant1.reverse ();
    cuadrant3.reverse ();
    switch (oABCD)
    {
        case A:
        case D:
            cuadrant2.reverse ();
            cuadrant4.reflect ();
            break;
        case C:
        case B:
            cuadrant2.reflect ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve23H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = D; t1 = H4;
            o2 = D; t2 = H0;
            o3 = B; t3 = H0;
            o4 = C; t4 = H4;
            break;
        case B:
            o1 = D; t1 = H4;
            o2 = A; t2 = H4;
            o3 = A; t3 = H0;
            o4 = C; t4 = H0;
            break;
        case C:
            o1 = D; t1 = H0;
            o2 = A; t2 = H4;
            o3 = B; t3 = H4;
            o4 = B; t4 = H0;
            break;
        case D:
            o1 = C; t1 = H0;
            o2 = A; t2 = H0;
            o3 = B; t3 = H4;
            o4 = C; t4 = H4;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);
    cuadrant1.reverse ();
    cuadrant3.reverse ();
    switch (oABCD)
    {
        case A:
        case D:
            cuadrant2.reverse ();
            cuadrant4.reflect ();
            break;
        case B:
        case C:
            cuadrant2.reflect ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve24H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = C; t1 = H0;
            o2 = D; t2 = H0;
            o3 = B; t3 = H0;
            o4 = B; t4 = H4;
            break;
        case B:
            o1 = C; t1 = H4;
            o2 = D; t2 = H0;
            o3 = A; t3 = H0;
            o4 = C; t4 = H0;
            break;
        case C:
            o1 = D; t1 = H0;
            o2 = D; t2 = H4;
            o3 = A; t3 = H0;
            o4 = B; t4 = H0;
            break;
        case D:
            o1 = C; t1 = H0;
            o2 = A; t2 = H0;
            o3 = A; t3 = H4;
            o4 = B; t4 = H0;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    cuadrant1.reverse ();
    cuadrant3.reverse ();
    switch (oABCD)
    {
        case A:
            cuadrant2.reverse ();
            cuadrant4.reflect ();
            break;
        case B:
        case D:
            cuadrant2.reverse ();
            cuadrant4.reverse ();
            break;
        case C:
            cuadrant2.reflect ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve25H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = D; t1 = H0;
            o2 = D; t2 = H0;
            o3 = B; t3 = H0;
            o4 = C; t4 = H4;
            break;
        case B:
            o1 = D; t1 = H4;
            o2 = A; t2 = H0;
            o3 = A; t3 = H0;
            o4 = C; t4 = H0;
            break;
        case C:
            o1 = D; t1 = H0;
            o2 = A; t2 = H4;
            o3 = B; t3 = H0;
            o4 = B; t4 = H0;
            break;
        case D:
            o1 = C; t1 = H0;
            o2 = A; t2 = H0;
            o3 = B; t3 = H4;
            o4 = C; t4 = H0;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    cuadrant1.reverse ();
    cuadrant3.reverse ();
    switch (oABCD)
    {
        case A:
            cuadrant2.reverse ();
            cuadrant4.reflect ();
            break;
        case B:
        case D:
            cuadrant2.reverse ();
            cuadrant4.reverse ();
            break;
        case C:
            cuadrant2.reflect ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve26H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = D; t1 = H0;
            o2 = D; t2 = H0;
            o3 = B; t3 = H0;
            o4 = B; t4 = H4;
            break;
        case B:
            o1 = C; t1 = H4;
            o2 = A; t2 = H0;
            o3 = A; t3 = H0;
            o4 = C; t4 = H0;
            break;
        case C:
            o1 = D; t1 = H0;
            o2 = D; t2 = H4;
            o3 = B; t3 = H0;
            o4 = B; t4 = H0;
            break;
        case D:
            o1 = C; t1 = H0;
            o2 = A; t2 = H0;
            o3 = A; t3 = H4;
            o4 = C; t4 = H0;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    cuadrant1.reverse ();
    cuadrant3.reverse ();
    switch (oABCD)
    {
        case A:
            cuadrant2.reverse ();
            cuadrant4.reflect ();
            break;
        case B:
        case D:
            cuadrant2.reverse ();
            cuadrant4.reverse ();
            break;
        case C:
            cuadrant2.reflect ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve27H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = C; t1 = H0;
            o2 = A; t2 = H0;
            o3 = A; t3 = H0;
            o4 = C; t4 = H4;
            break;
        case B:
            o1 = D; t1 = H4;
            o2 = D; t2 = H0;
            o3 = B; t3 = H0;
            o4 = B; t4 = H0;
            break;
        case C:
            o1 = C; t1 = H0;
            o2 = A; t2 = H4;
            o3 = A; t3 = H0;
            o4 = C; t4 = H0;
            break;
        case D:
            o1 = D; t1 = H0;
            o2 = D; t2 = H0;
            o3 = B; t3 = H4;
            o4 = B; t4 = H0;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            break;
        case B:
            cuadrant1.reflectAndReverse ();
            break;
        case C:
            break;
        case D:
            cuadrant3.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve28H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = C; t1 = H0;
            o2 = A; t2 = H0;
            o3 = A; t3 = H0;
            o4 = D; t4 = H4;
            break;
        case B:
            o1 = A; t1 = H4;
            o2 = D; t2 = H0;
            o3 = B; t3 = H0;
            o4 = B; t4 = H0;
            break;
        case C:
            o1 = C; t1 = H0;
            o2 = B; t2 = H4;
            o3 = A; t3 = H0;
            o4 = C; t4 = H0;
            break;
        case D:
            o1 = D; t1 = H0;
            o2 = D; t2 = H0;
            o3 = C; t3 = H4;
            o4 = B; t4 = H0;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            break;
        case B:
            cuadrant1.reflectAndReverse ();
            break;
        case C:
            break;
        case D:
            cuadrant3.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve29H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = B; t1 = H0;
            o2 = A; t2 = H0;
            o3 = A; t3 = H0;
            o4 = C; t4 = H4;
            break;
        case B:
            o1 = D; t1 = H4;
            o2 = D; t2 = H0;
            o3 = B; t3 = H0;
            o4 = B; t4 = H0;
            break;
        case C:
            o1 = C; t1 = H0;
            o2 = A; t2 = H4;
            o3 = D; t3 = H0;
            o4 = C; t4 = H0;
            break;
        case D:
            o1 = D; t1 = H0;
            o2 = D; t2 = H0;
            o3 = B; t3 = H4;
            o4 = A; t4 = H0;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            break;
        case B:
            cuadrant1.reflectAndReverse ();
            break;
        case C:
            break;
        case D:
            cuadrant3.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve30H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = B; t1 = H0;
            o2 = A; t2 = H0;
            o3 = A; t3 = H0;
            o4 = D; t4 = H4;
            break;
        case B:
            o1 = A; t1 = H4;
            o2 = C; t2 = H0;
            o3 = B; t3 = H0;
            o4 = B; t4 = H0;
            break;
        case C:
            o1 = C; t1 = H0;
            o2 = B; t2 = H4;
            o3 = D; t3 = H0;
            o4 = C; t4 = H0;
            break;
        case D:
            o1 = D; t1 = H0;
            o2 = D; t2 = H0;
            o3 = C; t3 = H4;
            o4 = A; t4 = H0;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            break;
        case B:
            cuadrant1.reflectAndReverse ();
            break;
        case C:
            break;
        case D:
            cuadrant3.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve31H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = C; t1 = H0;
            o2 = D; t2 = H0;
            o3 = B; t3 = H0;
            o4 = C; t4 = H4;
            break;
        case B:
            o1 = D; t1 = H4;
            o2 = D; t2 = H0;
            o3 = A; t3 = H0;
            o4 = C; t4 = H0;
            break;
        case C:
            o1 = D; t1 = H0;
            o2 = A; t2 = H4;
            o3 = A; t3 = H0;
            o4 = B; t4 = H0;
            break;
        case D:
            o1 = C; t1 = H0;
            o2 = A; t2 = H0;
            o3 = B; t3 = H4;
            o4 = B; t4 = H0;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    cuadrant1.reverse ();
    cuadrant3.reverse ();
    switch (oABCD)
    {
        case A:
            cuadrant2.reverse ();
            cuadrant4.reflect ();
            break;
        case B:
        case D:
            cuadrant2.reverse ();
            cuadrant4.reverse ();
            break;
        case C:
            cuadrant2.reflect ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve32H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = C; t1 = H1;
            o2 = A; t2 = H5;
            o3 = A; t3 = H5;
            o4 = C; t4 = H1;
            break;
        case B:
            o1 = D; t1 = H1;
            o2 = D; t2 = H1;
            o3 = B; t3 = H5;
            o4 = B; t4 = H5;
            break;
        case C:
            o1 = C; t1 = H5;
            o2 = A; t2 = H1;
            o3 = A; t3 = H1;
            o4 = C; t4 = H5;
            break;
        case D:
            o1 = D; t1 = H5;
            o2 = D; t2 = H5;
            o3 = B; t3 = H1;
            o4 = B; t4 = H1;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
        case D:
            cuadrant2.reflectAndReverse ();
            break;
        case B:
        case C:
            cuadrant4.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve33H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = C; t1 = H1;
            o2 = D; t2 = H5;
            o3 = B; t3 = H5;
            o4 = C; t4 = H1;
            break;
        case B:
            o1 = D; t1 = H1;
            o2 = D; t2 = H1;
            o3 = A; t3 = H5;
            o4 = C; t4 = H5;
            break;
        case C:
            o1 = D; t1 = H5;
            o2 = A; t2 = H1;
            o3 = A; t3 = H1;
            o4 = B; t4 = H5;
            break;
        case D:
            o1 = C; t1 = H5;
            o2 = A; t2 = H5;
            o3 = B; t3 = H1;
            o4 = B; t4 = H1;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    cuadrant1.reverse ();
    cuadrant3.reverse ();
    switch (oABCD)
    {
        case A:
        case D:
            cuadrant2.reflect ();
            cuadrant4.reverse ();
            break;
        case B:
        case C:
            cuadrant2.reverse ();
            cuadrant4.reflect ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve34H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = C; t1 = H5;
            o2 = A; t2 = H5;
            o3 = A; t3 = H5;
            o4 = C; t4 = H1;
            break;
        case B:
            o1 = D; t1 = H1;
            o2 = D; t2 = H5;
            o3 = B; t3 = H5;
            o4 = B; t4 = H5;
            break;
        case C:
            o1 = C; t1 = H5;
            o2 = A; t2 = H1;
            o3 = A; t3 = H5;
            o4 = C; t4 = H5;
            break;
        case D:
            o1 = D; t1 = H5;
            o2 = D; t2 = H5;
            o3 = B; t3 = H1;
            o4 = B; t4 = H5;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            cuadrant2.reflectAndReverse ();
            break;
        case B:
            cuadrant2.reflectAndReverse ();
            cuadrant4.reflectAndReverse ();
            break;
        case C:
            cuadrant4.reflectAndReverse ();
            break;
        case D:
            cuadrant2.reflectAndReverse ();
            cuadrant4.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve35H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = D; t1 = H5;
            o2 = D; t2 = H5;
            o3 = B; t3 = H5;
            o4 = C; t4 = H1;
            break;
        case B:
            o1 = D; t1 = H1;
            o2 = A; t2 = H5;
            o3 = A; t3 = H5;
            o4 = C; t4 = H5;
            break;
        case C:
            o1 = D; t1 = H5;
            o2 = A; t2 = H1;
            o3 = B; t3 = H5;
            o4 = B; t4 = H5;
            break;
        case D:
            o1 = C; t1 = H5;
            o2 = A; t2 = H5;
            o3 = B; t3 = H1;
            o4 = C; t4 = H5;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    cuadrant1.reverse ();
    cuadrant3.reverse ();
    switch (oABCD)
    {
        case A:
            cuadrant2.reflect ();
            cuadrant4.reverse ();
            break;
        case B:
        case D:
            cuadrant2.reflect ();
            cuadrant4.reflect ();
            break;
        case C:
            cuadrant2.reverse ();
            cuadrant4.reflect ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve36H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = B; t1 = H5;
            o2 = A; t2 = H5;
            o3 = A; t3 = H5;
            o4 = C; t4 = H1;
            break;
        case B:
            o1 = D; t1 = H1;
            o2 = C; t2 = H5;
            o3 = B; t3 = H5;
            o4 = B; t4 = H5;
            break;
        case C:
            o1 = C; t1 = H5;
            o2 = A; t2 = H1;
            o3 = D; t3 = H5;
            o4 = C; t4 = H5;
            break;
        case D:
            o1 = D; t1 = H5;
            o2 = D; t2 = H5;
            o3 = B; t3 = H1;
            o4 = A; t4 = H5;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            cuadrant1.reflectAndReverse ();
            cuadrant2.reflectAndReverse ();
            break;
        case B:
            cuadrant4.reflectAndReverse ();
            break;
        case C:
            cuadrant3.reflectAndReverse ();
            cuadrant4.reflectAndReverse ();
            break;
        case D:
            cuadrant2.reflectAndReverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve37H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = C; t1 = H5;
            o2 = D; t2 = H5;
            o3 = B; t3 = H5;
            o4 = C; t4 = H1;
            break;
        case B:
            o1 = D; t1 = H1;
            o2 = D; t2 = H5;
            o3 = A; t3 = H5;
            o4 = C; t4 = H5;
            break;
        case C:
            o1 = D; t1 = H5;
            o2 = A; t2 = H1;
            o3 = A; t3 = H5;
            o4 = B; t4 = H5;
            break;
        case D:
            o1 = C; t1 = H5;
            o2 = A; t2 = H5;
            o3 = B; t3 = H1;
            o4 = B; t4 = H5;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
            cuadrant1.reflect ();
            cuadrant2.reflect ();
            cuadrant3.reverse ();
            cuadrant4.reverse ();
            break;
        case B:
            cuadrant1.reverse ();
            cuadrant2.reverse ();
            cuadrant3.reverse ();
            cuadrant4.reflect ();
            break;
        case C:
            cuadrant1.reverse ();
            cuadrant2.reverse ();
            cuadrant3.reflect ();
            cuadrant4.reflect ();
            break;
        case D:
            cuadrant1.reverse ();
            cuadrant2.reflect ();
            cuadrant3.reverse ();
            cuadrant4.reverse ();
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve38H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = B; t1 = H3;
            o2 = A; t2 = H5;
            o3 = A; t3 = H5;
            o4 = C; t4 = H1;
            break;
        case B:
            o1 = D; t1 = H1;
            o2 = C; t2 = H3;
            o3 = B; t3 = H5;
            o4 = B; t4 = H5;
            break;
        case C:
            o1 = C; t1 = H5;
            o2 = A; t2 = H1;
            o3 = D; t3 = H3;
            o4 = C; t4 = H5;
            break;
        case D:
            o1 = D; t1 = H5;
            o2 = D; t2 = H5;
            o3 = B; t3 = H1;
            o4 = A; t4 = H3;
            break;
        default:
            throw HilbertBadOrientation();
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    switch (oABCD)
    {
        case A:
        case D:
            cuadrant2.reflectAndReverse ();
            break;
        case B:
        case C:
            cuadrant4.reflectAndReverse ();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}
void HilbertCurve::BuildCurve39H()
{
    hsize w1, w2, h1, h2;
    //Even partition of the quasiquare
    w2 = width ()/2;
    w1 = width () - w2;
    h2 = height ()/2;
    h1 = height ()- h2;

    Orientation o1, o2, o3, o4;
    CurveType t1, t2, t3, t4;
    switch (oABCD)
    {
        case A:
            o1 = D; t1 = H3;
            o2 = D; t2 = H5;
            o3 = B; t3 = H5;
            o4 = C; t4 = H1;
            break;
        case B:
            o1 = D; t1 = H1;
            o2 = A; t2 = H3;
            o3 = A; t3 = H5;
            o4 = C; t4 = H5;
            break;
        case C:
            o1 = D; t1 = H5;
            o2 = A; t2 = H1;
            o3 = B; t3 = H3;
            o4 = B; t4 = H5;
            break;
        case D:
            o1 = C; t1 = H5;
            o2 = A; t2 = H5;
            o3 = B; t3 = H1;
            o4 = C; t4 = H3;
            break;
    }
    //Principal 4 quasisquares
    HilbertCurve cuadrant1(w1, h1, t1, coord + HPoint(0, 0), o1);
    HilbertCurve cuadrant2(w1, h2, t2, coord + HPoint(0, h1), o2);
    HilbertCurve cuadrant3(w2, h2, t3, coord + HPoint(w1, h1), o3);
    HilbertCurve cuadrant4(w2, h1, t4, coord + HPoint(w1, 0), o4);

    cuadrant1.reverse ();
    cuadrant3.reverse ();
    switch (oABCD)
    {
        case A:
        case D:
            cuadrant2.reflect ();
            cuadrant4.reverse ();
            break;
        case B:
        case C:
            cuadrant2.reverse ();
            cuadrant4.reflect ();
            break;
    }
    joinCurve (cuadrant1, cuadrant2, cuadrant3, cuadrant4);
}

/*!
  \class HilbertCurve
  \inmodule hilbertlib
  \ingroup hcurve
  \brief The HilbertCurve class represent a Hilbert curve.

  The HilbertCurve class represent a Hilbert curve with a given dimension, orientations
  and start coordinate. Internally it inherits from QuasiSquare in a protected way so
  it hides the QuasiSquare interface. HilbertCurve class can compute the 40 different set
  of homogeneous and no homogeneus curves given by E. Estevez-Rams et al. in the article
  \e {"Hilbert Curves in two dimensions"}
 */
/*!
  Default constructors.
  Create a single point curve.
*/
HilbertCurve::HilbertCurve():
    HilbertCurve(1,1,H0,0,A)
{}
/*!
  \brief General constructor

  Construct the Hilbert Curve defined by \a type with dimensions \a width, \a height; with the
  orientation given by \a orientation and starting at \a coord. The flag \a differenceCurve is needed
  in order to plot the Difference Map, by default is true.
*/
HilbertCurve::HilbertCurve(hsize width, hsize height, CurveType type, HPoint coord, Orientation orientation, bool differenceCurve):
    QuasiSquare(height, width,  coord, orientation),
    m_type(type)
{
    Build ();
    if(differenceCurve)
    {
        BuildDifference ();
        reflectY ();
    }

}
/*!
  \fn HilbertCurve::HilbertCurve(const HilbertCurve &curve)
  \brief Copy Constructor
  Construct a copy of the given \a curve
*/
HilbertCurve::HilbertCurve(const HilbertCurve &curve):
    QuasiSquare(curve.n, curve.m, curve.coord, curve.oABCD),
    m_type(curve.m_type),
    m_curve(curve.m_curve),
    m_mean_difference(curve.meanDifference ())
{
}
/*!
  Returns the mean differnece of the curve. BuildDifference() must be called first.
*/
hfloat HilbertCurve::meanDifference() const
{
     return m_mean_difference;
}

/*!
  \fn void HilbertCurve::reflectX()
  \brief Perform an horizontal reflection of the curve
*/
void HilbertCurve::reflectX()
{
    for_each_parallel(m_curve.begin (), m_curve.end (),
                  [&](HPoint &p)
                  {
                      p.X (m  - 1 - p.X() + 2 * coord.X ());
                  });
//    for(auto iter = m_curve.begin(); iter != m_curve.end(); ++iter)
//    {
//        iter->X(width ()-1-iter->X() + 2*coord.X ());

//    }
}
/*!
  \fn void HilbertCurve::reflectY()
  \brief Perform a vertical reflection of the curve
*/
void HilbertCurve::reflectY ()
{
    for_each_parallel(m_curve.begin (), m_curve.end (),
                  [&](HPoint &p)
                  {
                      p.Y (n  - 1 - p.Y() + 2 * coord.Y ());
                  });
//    for(auto iter = m_curve.begin(); iter != m_curve.end(); ++iter)
//    {
//        iter->Y(height ()-1-iter->Y() + 2*coord.Y());
//    }
}
/*!
  \fn void HilbertCurve::joinCurve(HilbertCurve &cuadrant1, HilbertCurve &cuadrant2, HilbertCurve &cuadrant3, HilbertCurve &cuadrant4)
  \brief Join the four cuadrants according to orientation

  This function is used for HilbertCurve for joining the 4 cuadrants according to the
  orientation. The cuadrants order is as follow:
  \list
  \li \a cuadrant1  Lower Left Cuadrant
  \li \a cuadrant2  Upper Left Cuadrant
  \li \a cuadrant3  Upper Right Cuadrant
  \li \a cuadrant3  Lower Right Cuadrant
  \endlist
*/
void HilbertCurve::joinCurve(HilbertCurve &cuadrant1, HilbertCurve &cuadrant2, HilbertCurve &cuadrant3, HilbertCurve &cuadrant4)
{
    m_curve.clear ();
    m_curve.reserve (width ()*height ());
    try
    {
        switch (oABCD)
        {
            case A:
                //Building the curve
                m_curve.insert(m_curve.end(), cuadrant1.m_curve.begin(), cuadrant1.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant2.m_curve.begin(), cuadrant2.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant3.m_curve.begin(), cuadrant3.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant4.m_curve.begin(), cuadrant4.m_curve.end());
                break;
            case B:
                //Building the curve
                m_curve.insert(m_curve.end(), cuadrant1.m_curve.begin(), cuadrant1.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant4.m_curve.begin(), cuadrant4.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant3.m_curve.begin(), cuadrant3.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant2.m_curve.begin(), cuadrant2.m_curve.end());
                break;
            case C:
                //Building the curve
                m_curve.insert(m_curve.end(), cuadrant3.m_curve.begin(), cuadrant3.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant4.m_curve.begin(), cuadrant4.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant1.m_curve.begin(), cuadrant1.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant2.m_curve.begin(), cuadrant2.m_curve.end());
                break;
            case D:
                //Building the curve
                m_curve.insert(m_curve.end(), cuadrant3.m_curve.begin(), cuadrant3.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant2.m_curve.begin(), cuadrant2.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant1.m_curve.begin(), cuadrant1.m_curve.end());
                m_curve.insert(m_curve.end(), cuadrant4.m_curve.begin(), cuadrant4.m_curve.end());
                break;
            default:
                throw HilbertBadOrientation();
                break;
        }
    }
    catch(HilbertBadOrientation &o)
    {
        throw o;
    }
    catch(std::bad_alloc &ba)
    {
        throw HilbertBadAlloc();
    }
}
/*!
  \fn void HilbertCurve::reflectAndReverse()
  \brief Perform reflect() and reverse() continusly
  \sa reverse(), reflect()
*/
void HilbertCurve::reflectAndReverse()
{
    reflect ();
    reverse ();
}
/*!
  \fn void HilbertCurve::reverse()
  \brief Reverse curve
*/
void HilbertCurve::reverse()
{
    reverse_parallel(m_curve.begin (), m_curve.end ());
}
/*!
  \fn void HilbertCurve::reflect()
  \brief Reflect curve according to orientation

  Reflect the curve based on the orientation. If orientations
  is \c A or \c C reflectX() will be called, else if orientation
  is \c B or \c D reflectY() will be called instead.
*/
void HilbertCurve::reflect()
{
    if(oABCD == A || oABCD == C)
        reflectX ();
    else if(oABCD == B || oABCD == D)
        reflectY ();
    else
        throw HilbertBadOrientation();

}
