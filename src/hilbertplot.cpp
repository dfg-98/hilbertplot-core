/*!
   \headerfile "hilbertplot.h"

   \title Hilbert Plot Declaration

   \brief The "hilbertplot.h" header define HilbertPlot class
 */
#include "hilbertplot.h"
#include <cmath>
#include <fftw3.h>
#include <limits>
#include <iostream>


/*!
   \class HilbertPlot
   \since 0.1
   \inmodule hilbertlib
   \ingroup hcurve
   \ingroup hdata
   \brief The \c HilbertPlot class define a linear data vector
   with an associated HilbertCurve.

   The \c HilbertPlot class implements a read-only class representing a HilbertCurve
    with a DataSequence. It inheriths public from HilbertCurve, implementing its whole interfaces.
  It also have functionallity for accessing internal DataSequence values.

*/

/*!

  \brief Default Construnctor

       Construct the HilbertPlot as an empty data vector.
 */
HilbertPlot::HilbertPlot():
    HilbertPlot(DataSequence(), 0, 0, H0)
{}
/*!
   \brief General Constructor

     Constructs the \c HilbertPlot using \a data. If \a width and \a height are given greather
     than zero data will be shrinked or expanding according to given values. If zero is given
     best dimension will be computed. The \c HilbertCurve used is given by \a type.
 */
HilbertPlot::HilbertPlot(const DataSequence &data, hsize width, hsize height, CurveType type):
    HilbertCurve (constructCurve (data.size (), width, height, type)),
    m_data(data)
{
    while(m_data.size () < width * height)
    {
        m_data.push_back(0);
    }
    while(m_data.size () > width * height)
    {
        m_data.pop_back();
    }
    m_plotToCurve = std::vector<std::vector<hint>>(width, std::vector<hint>(height, 0));
    for(auto point = HilbertCurve::begin (); point != HilbertCurve::end (); ++point)
    {
        m_plotToCurve[point->X ()][point->Y()] = point->index;
    }
    if(m_data.size () > 0)
    {
        m_min = m_data.min ();
        m_max = m_data.max ();
    }
    else
    {
        m_min = 0;
        m_max = 0;
    }
}
/*!
  \brief Copy Constructor
   Construct a copy of \a hilbertplot
 */
HilbertPlot::HilbertPlot(const HilbertPlot &hilbertplot):
     HilbertCurve(hilbertplot),
     m_data(hilbertplot.m_data),
    m_min(hilbertplot.m_min),
    m_max(hilbertplot.m_max),
    m_plotToCurve(hilbertplot.m_plotToCurve)
{
}
/*!
  \brief Constant acces operator.
  Returns a read-only (constant) HPoint reference at \a index.
*/
std::vector<HPoint>::const_reference HilbertPlot::operator [](std::vector<HPoint>::size_type index) const
{
    if(index >= HilbertCurve::lenght ())
        throw HilbertIndexOutOfRange();
    return HilbertCurve::operator [] (index);
}
/*!
  \brief Constant access operator.
  Returns a read only (constant) HPoint reference at \a index.
  \sa operator[]
*/
std::vector<HPoint>::const_reference HilbertPlot::at(std::vector<HPoint>::size_type index)
{
    if(index >= HilbertCurve::lenght ())
        throw HilbertIndexOutOfRange();
    return HilbertCurve::operator [] (index);
}
/*!
  \overload at()
  Returns the HPoint at positions \a x and \a y on the plot.
  \sa valueAt()
*/
std::vector<HPoint>::const_reference HilbertPlot::at(std::vector<HPoint>::size_type x, std::vector<HPoint>::size_type y)
{
    if(x >= width () || y >= height ())
        throw HilbertIndexOutOfRange();
    int index = m_plotToCurve[x][y];
    return HilbertCurve::operator [] (index);
}

/*!
  Returns the value stored at given \a index of the curve.
*/
hfloat HilbertPlot::valueAt(std::vector<hfloat>::size_type index) const
{
    if(index >= m_data.size ())
        throw HilbertIndexOutOfRange();
    return m_data.at (index);
}
/*!
  \overload valueAt()
  Returns the value stored in the plot at coordinates \a x and \a y.
  \note HilbertIndexOutOfRange() exception is thrown if the given indexes aren't valid.
*/
hfloat HilbertPlot::valueAt(std::vector<hfloat>::size_type x, std::vector<hfloat>::size_type y) const
{
    if(x >= width () || y >= height ())
        throw HilbertIndexOutOfRange();
    hsize index = m_plotToCurve[x][y];
    if(index >= m_data.size ())
        throw HilbertIndexOutOfRange();
    return m_data.at (index);
}

hfloat HilbertPlot::valueNormalizedAt(std::vector<hfloat>::size_type index) const
{
    if(index >= m_data.size ())
        throw HilbertIndexOutOfRange();
    return (m_data.at (index) - m_min) * (m_min == m_max? 0.0 : 1.0/(m_max - m_min));
}

hfloat HilbertPlot::valueNormalizedAt(std::vector<hfloat>::size_type x, std::vector<hfloat>::size_type y) const
{
    if(x >= width () || y >= height ())
        throw HilbertIndexOutOfRange();
    hsize index = m_plotToCurve[x][y];
    return valueNormalizedAt (index);
}

void HilbertPlot::replaceValueAt(std::vector<hfloat>::size_type index, hfloat value)
{
    if(index >= m_data.size ())
        throw HilbertIndexOutOfRange();
    m_data[index] = value;
    m_max = m_data.max ();
    m_min = m_data.min ();
}

void HilbertPlot::replaceValueAt(std::vector<hfloat>::size_type x, std::vector<hfloat>::size_type y, hfloat value)
{
    if(x >= width () || y >= height ())
        throw HilbertIndexOutOfRange();
    hsize index = m_plotToCurve[x][y];
    replaceValueAt (index, value);
}
/*!
  \brief Return the index of the curve corresponding to the \a x , \a y coordinate.
*/
hint HilbertPlot::indexOf(hint x, hint y) const
{
    if(x >= width () || y >= height ())
        throw HilbertIndexOutOfRange();
    return m_plotToCurve[x][y];
}

/*!

  An HImage is a matrix with dimension
  \c width() x \c height() represented as a vector
  of vector of \c hfloat. Each value is a normalized
  value in range [0-1] representing the pixel intensity.
  It's represented this way so it can be used for externals
  tools to generates images with a given color gradient.

  \note If \a threshold is given greather than 0 difference map will be computed
  assigning a value of 2 to the points with difference value greather than \a threshold
*/
HImage HilbertPlot::generateImage(hfloat threshold)
{
    std::vector<std::vector<hfloat>> image(width (), std::vector<hfloat>(height (), 0));
    hfloat minmax = m_max == m_min ? 0.0 : 1.0/(m_max - m_min);
    if(threshold > 0)
    {
        for(auto point = HilbertCurve::cbegin (); point != HilbertCurve::cend (); ++point)
        {
            hfloat value = (valueAt (point->X(), point->Y()) - m_min) * minmax;
            if(point->DifferenceValue () / meanDifference ()> threshold)
                value = 2;
            image[point->X()][point->Y ()] = value;
        }
    }
    else
    {
        for(auto point = HilbertCurve::cbegin (); point != HilbertCurve::cend (); ++point)
        {
            hfloat value = valueAt (point->X(), point->Y()) * minmax;
            image[point->X()][point->Y ()] = value;
        }
    }
    return image;
}

/*!
  Returns a copy of the parent DataSequence
*/
DataSequence HilbertPlot::dataCopy() const
{
    return m_data;
}
/*!
  Replace the current data with the given \a data.
  \note An exception is thrown if data sizes aren't equals.
*/
void HilbertPlot::replaceData(const DataSequence &data)
{
    if(m_data.size () != data.size ())
    {
        throw HilbertBadSize();
    }

    hfloat max = data.max ();
    hfloat min = data.min ();
    hfloat minmax = 1.0/(max - min);
    m_data.clear();
    m_data.reserve (data.size ());
    for(const hfloat &val : data)
    {
        m_data.push_back ((val - min) * minmax);
    }
}
/*!
  \brief Compute the Fourier Transform of the 2D HilbertPlot.
  \note If \a logflag is set to true. Plot be computed in logarithm base.
*/
DataSequence HilbertPlot::hpFourierTransform(bool logflag) const
{
    if(m_data.size () == 0) throw HilbertBadOperation();
    double *datainput;
    fftw_complex *dataoutput;
    fftw_plan p;
    uint width = this->width ();
    uint height = this->height ();
    uint data_size = width*height;
    uint y=0;
    uint x=0;
    uint w2=width/2;
    uint i=0;
    std::vector<double> output;

    double max, max2;
    max = max2 =std::numeric_limits<double>::min();
    double min =std::numeric_limits<double>::max();

    uint data_size2=height*(w2+1);

    try{
        datainput = (double*) fftw_malloc(sizeof(double) * data_size);
        dataoutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (data_size2+2));
        p = fftw_plan_dft_r2c_2d(height, width, datainput, dataoutput, FFTW_ESTIMATE);

        for(y = 0; y < height; ++y)
        {
            for(x = 0; x < width; ++x)
            {
                double val = valueAt (x, y);
                datainput[i++]= val;
            }
        }
        fftw_execute(p); // the transform

        output.clear();
        output.assign (data_size, 0);
    }
    catch (std::bad_alloc& ba)
    {
            throw HilbertBadAlloc();
    }

    double v;
    for(i=0; i < data_size2; ++i)
    {
        v=dataoutput[i][0]*dataoutput[i][0]+dataoutput[i][1]*dataoutput[i][1];
        if(max < v)
            max=v;
        if  (min > v)
            min=v;
        if(max2 < v && v < max)
            max2=v;
    }
    double maxmin=max2-min;
    double wdf=0;
    if(logflag)
        maxmin=log(maxmin);

    width++;
    // arranging transform data into redundant centered output
    for(y=0; y < height; ++y)
    {
        for(x=0; x <= w2; ++x)
        {
            uint index1 = this->indexOf (x, y);
            uint index2 = this->indexOf (this->width ()-x-1, y);
            
            if(!logflag)
            {
                wdf=dataoutput[y*(w2+1)+x][0]*dataoutput[y*(w2+1)+x][0]
                    +dataoutput[y*(w2+1)+x][1]*dataoutput[y*(w2+1)+x][1];
                if(wdf >= max)
                    wdf=max2;
                output[index1]=output[index2]= (wdf-min)/maxmin;
            }
            else
            {
                wdf=dataoutput[y*(w2+1)+x][0]*dataoutput[y*(w2+1)+x][0]+dataoutput[y*(w2+1)+x][1]*dataoutput[y*(w2+1)+x][1];
                if(wdf > max2)
                    wdf=max2;
                output[index1]=output[index2]= log(wdf-min+1)/maxmin;
            }
        }
        if(!logflag)
        {
            wdf=dataoutput[y*(w2+1)][0]*dataoutput[y*(w2+1)][0]+dataoutput[y*(w2+1)][1]*dataoutput[y*(w2+1)][1];
            if(wdf >= max)
                wdf=max2;
            output[this->indexOf (w2, y)]= (wdf-min)/maxmin;
        }
        else
        {
            wdf=dataoutput[y*(w2+1)][0]*dataoutput[y*(w2+1)][0]+dataoutput[y*(w2+1)][1]*dataoutput[y*(w2+1)][1];
            if(wdf > max2)
                wdf=max2;
            if(wdf-min >0)
                output[this->indexOf (w2, y)]= log(wdf-min)/maxmin;
        }
    }
    fftw_destroy_plan(p);
    fftw_free(datainput);
    fftw_free(dataoutput);
    DataSequence dataOutput(output.size ());
    i = 0;
    for(auto val : output)
    {
        dataOutput[i++] = val;
    }
    return dataOutput;
}


/*!
  \brief Compute best dimensions to adjust the given \a lenght

       For HilbertPlot is needed that data lenght be adjusted to fit a square or quasisquere.
       The best dimension is given by the values of width and height that guarantee the less
       data loss.
 */
std::pair<hsize, hsize> HilbertPlot::bestDimensions(hsize lenght)
{
    std::pair<hsize, hsize> dim;
    hfloat sq = std::sqrt(lenght);
    if(hsize(sq)*hsize(sq) == lenght)//If it's a perfect square
    {
        dim.first = dim.second = sq;
        return dim;
    }
    // Find the minimum for data loss
    hsize f = std::floor(sq);
    hsize c = std::ceil(sq);
    int min1 = std::abs(int(lenght - f*f));
    int min2 = std::abs (int(lenght - c*c));
    int min3 = std::abs (int(lenght - f*c));
    if(min1 < min2)
    {
        if(min1 < min3)
        {
            dim.first = dim.second = f;
            return dim;
        }
        dim.first = c;
        dim.second = f;
        return dim;
    }
    if(min2 < min3)
    {
        dim.first = dim.second = c;
        return dim;
    }
    dim.first = c;
    dim.second = f;
    return dim;
}
/*!
  Generate the HilbertCurve for the constructor
*/
const HilbertCurve HilbertPlot::constructCurve(hsize lenght, hsize &width, hsize &height, CurveType type)
{
    if(width == 0 || height == 0)
    {
        auto dim = bestDimensions (lenght);
        width = dim.first;
        height = dim.second;

    }
    return HilbertCurve(width, height, type, 0, A, true);
}
