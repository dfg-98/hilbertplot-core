/*!
  * \headerfile "datasequence.h"
  *
  *  \title DataSequence
  *
  *  \brief The "datasequence.h" header define DataSequence class.
  *
  */
#include "datasequence.h"

#include <fftw3.h>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <exception>
#include <numeric>
#include <sstream>
#include <fstream>
#include <iostream>

hfloat zlog(hfloat val);

/*!
  \class DataSequence

  \brief Represents a linear data vector
  \inherits std::vector<hfloat>
  \ingroup hilbert
  \inmodule hilbertlib

  \c DataSequence class inheriths from std::vector<hfloat> in a public way.
  Implements the whole interface of std::vector, with
  overloading operators for perfoming data operations as
  member by member arithmetics and logical operations, and accessors operators.
*/
/*!
  Default constructor. Constructs an empty data.
*/
DataSequence::DataSequence(): std::vector<hfloat>()
{}
/*!
  General constructor.
  Constructs the data with \a elements copy of \a defaultValue.
*/
DataSequence::DataSequence(int elements, hfloat defaultValue):
    std::vector<hfloat>(elements, defaultValue)
{}
/*!
  Constructs the DataSequence with values provided by \a data.
*/
DataSequence::DataSequence(const std::vector<hfloat> &data):
    std::vector<hfloat>(data)
{}

/*!
  Copy constructor. Create a copy of \a data
*/
DataSequence::DataSequence(const DataSequence &data):
    std::vector<hfloat>(data)
{}
/*!
  Move constructor. Transfers \a data
*/
DataSequence::DataSequence(DataSequence &&data):
    std::vector<hfloat>(data)
{}
DataSequence::DataSequence(std::initializer_list<hfloat> data):
    std::vector<hfloat>(data)
{}

DataSequence &DataSequence::operator=(const DataSequence &other)
{
    this->clear ();
    this->reserve (other.size ());
    for(auto &val: other)
        this->push_back (val);
    return *this;
}
/*!
  Returns a new \c DataSequence as a sum of elements one by one
  between \c this and \a d.

  \sa DataSequence::operator+()
*/
DataSequence DataSequence::operator+(const DataSequence &d) const
{
    DataSequence newData;
    for(size_t i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(this->at(i));
        }
        else
        {
            newData.push_back(this->at (i) + d[i]);
        }
    }
    return newData;
}
/*!
  \overload DataSequence::operator+()

  Returns a new \c DataSequence resulted from the sum of elements
with \a value given.
  \note Returned data size is determined by first data operand
 \sa DataSequence::operator+()
*/
DataSequence DataSequence::operator + (const hfloat &value)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back(this->at (i) + value);
    }
    return newData;
}
/*!
  Returns a new \c DataSequence as a substractions of elements one by one
  between \c this and \a d.
  \note Returned data size is determined by first data operand
  \sa DataSequence::operator+()
*/
DataSequence DataSequence::operator-(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(this->at (i));
        }
        else
        {
            newData.push_back(this->at (i) - d[i]);
        }
    }
    return newData;
}
/*!
  \overload DataSequence::operator-()

  Returns a new \c DataSequence resulted from the diference of elements
with \a value given.
  \note Returned data size is determined by first data operand
 \sa DataSequence::operator+()
*/
DataSequence DataSequence::operator-(const hfloat &value)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back(this->at (i) - value);
    }
    return newData;
}
/*!
  \fn DataSequence DataSequence::operator*(const DataSequence &d) const

  Returns a new \c DataSequence as a multiplication of elements one by one between
  \c this and \a d.

  \note Returned data size is determined by first data operand
  \sa DataSequence::operator+()
*/
DataSequence DataSequence::operator*(const DataSequence &d) const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(this->at(i));
        }
        else
        {
            newData.push_back(this->at(i) * d[i]);
        }
    }
    return newData;
}
/*!
  \overload DataSequence::operator*()

  Returns a new \c DataSequence resulted from a multiplication of elements
with \a value given.
  \note Returned data size is determined by first data operand
 \sa DataSequence::operator*()
*/
DataSequence DataSequence::operator * (const hfloat &value) const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back(this->at(i) * value);
    }
    return newData;
}
/*!
  Returns a DataSequence resulting from perfom division between data and
  \a d elements. If zero division detected HilbertZeroDivision exception will be throw.
*/
DataSequence DataSequence::operator/(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(this->at(i));
        }
        else
        {
            if(d[i] == 0)
            {
                throw HilbertZeroDivision();
            }
            newData.push_back(this->at(i) * d[i]);
        }
    }
    return newData;
}
/*!
  \overload operator/()

*/
DataSequence DataSequence::operator / (const hfloat &val)  const
{
    if(val == 0.0)
    {
        throw HilbertZeroDivision();
    }
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back(this->at(i) * val);
    }
    return newData;
}
/*!
  Returns a DataSequence resulting from perfom potentiation between data and
  \a d elements.
*/
DataSequence DataSequence::operator^(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(this->at(i));
        }
        else
        {
            newData.push_back(std::pow(this->at(i),  d[i]));
        }
    }
    return newData;
}
/*!
  \overload operator^()
*/
DataSequence DataSequence::operator ^ (const hfloat &val)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back(std::pow(this->at(i), val));
    }
    return newData;
}
/*!
  \overload DataSequence::operator^()
*/
DataSequence operator^(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        newData.push_back(std::pow(val, d[i]));
    }
    return newData;
}
/*!
  Returns a DataSequence with \c 1 where elements of \c this and \a d are equals,
  \c 0 otherwise.
  \sa hammingDistance()
*/
DataSequence DataSequence::operator==(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(0);
        }
        else
        {
            newData.push_back(this->at(i) == d[i]);
        }
    }
    return newData;
}
/*!
  \overload operator==()
*/
DataSequence DataSequence::operator ==(const hfloat &val)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back(this->at(i) == val);
    }
    return newData;
}
/*!
  \relates DataSequence
  Returns a DataSequence with the same size as \a d where elemnts
  will be \c 1 if the corresponding element in \a d and \a val are equal,
    \c 0 otherwise.
*/
DataSequence operator==(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        newData.push_back(val == d[i]);
    }
    return newData;
}
/*!
  Returns a DataSequence with the same size as \c this where elemnts
  will be \c 1 if the corresponding element in \c this and \a d are different,

  \sa hammingDistance(), operator==()
*/
DataSequence DataSequence::operator!=(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(1);
        }
        else
        {
            newData.push_back(this->at(i) != d[i]);
        }
    }
    return newData;
}
/*!
  Returns a DataSequence with the same size as \c this where elemnts
  will be \c 1 if the corresponding element in data is different of \a val,
    \c 0 otherwise.
*/
DataSequence DataSequence::operator !=(const hfloat &val)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back(this->at(i) != val);
    }

    return newData;
}
/*!
  Returns a DataSequence with the same size as \a d where elemnts
  will be \c 1 if the corresponding element in \a d is different of \a val,
    \c 0 otherwise.
*/
DataSequence operator!=(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        newData.push_back(val != d[i]);
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing greater than comparsion
    operation between data and \a d values.
*/
DataSequence DataSequence::operator>(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(1);
        }
        else
        {
            newData.push_back(this->at(i) > d[i]);
        }
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing greater than comparsion
    operation between data values and \a val.
*/
DataSequence DataSequence::operator >(const hfloat &val)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back(this->at(i) > val);
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing greater than comparsion
    operation between data and \a d values.
*/
DataSequence operator>(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        newData.push_back(val > d[i]);
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing less than comparsion
    operation between data and \a d values.
*/
DataSequence DataSequence::operator<(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(0);
        }
        else
        {
            newData.push_back(this->at(i) < d[i]);
        }
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing less than comparsion
    operation between data values and \a val.
*/
DataSequence DataSequence::operator <(const hfloat &val)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back(this->at(i) < val);
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing less comparsion
    operation between data and \a d values.
*/
DataSequence operator<(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        newData.push_back(val < d[i]);
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing greater-equal than comparsion
    operation between data and \a d values.
*/
DataSequence DataSequence::operator>=(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(1);
        }
        else
        {
            newData.push_back(this->at(i) >= d[i]);
        }
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing greater-equal than comparsion
    operation between data values and \a val.
*/
DataSequence DataSequence::operator >=(const hfloat &val)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back(this->at(i) >= val);
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing greater-equal than comparsion
    operation between \a val and \a d values.
*/
DataSequence operator>=(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        newData.push_back(val >= d[i]);
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing less-equal than comparsion
    operation between data and \a d values.
*/
DataSequence DataSequence::operator<=(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(0);
        }
        else
        {
            newData.push_back(this->at(i) <= d[i]);
        }
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing less-equal comparsion
    operation between data values and \a val.
*/
DataSequence DataSequence::operator <=(const hfloat &val)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back(this->at(i) <= val);
    }

    return newData;
}
/*!
    Returns a DataSequence resulting for performing less-equal comparsion
    operation between \a val and \a d values.
*/
DataSequence operator<=(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        newData.push_back(val <= d[i]);
    }

    return newData;
}
/*!
  Perform and AND operation between data and \a d values.
  Returns the resulting DataSequence.
  \note A value is \c true if is greather than 0
*/
DataSequence DataSequence::operator&&(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(0);
        }
        else
        {
            newData.push_back((this->at(i) > 0) && (d[i] > 0));
        }
    }

    return newData;
}
/*!
  Perform and AND operation between data values and \a val.
  Returns the resulting DataSequence
  \note A value is \c true if is greather than 0
*/
DataSequence DataSequence::operator &&(const hfloat &val)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back((this->at(i) > 0) && (val > 0));
    }

    return newData;
}
/*!
  Perform and AND operation between \a val and \a d values.
  Returns the resulting DataSequence
  \note A value is \c true if is greather than 0
*/
DataSequence operator&&(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        newData.push_back((val > 0) && (d[i] > 0));
    }

    return newData;
}
/*!
  Perform and OR operation between data and \a d values.
  Returns the resulting DataSequence
  \note A value is \c true if is greather than 0
*/
DataSequence DataSequence::operator|(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(1);
        }
        else
        {
            newData.push_back((this->at(i) > 0) || (d[i] > 0));
        }
    }

    return newData;
}
/*!
  \overload operator|()
    Perform and OR operation between data values and \a val.
  Returns the resulting DataSequence
  \note A value is \c true if is greather than 0
*/
DataSequence DataSequence::operator |(const hfloat &val)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back((this->at(i) > 0) || (val > 0));
    }

    return newData;
}
/*!
  \overload operator|()
  Perform and OR operation between \a val and \a d values.
  Returns the resulting DataSequence
  \note A value is \c true if is greather than 0
*/
DataSequence operator|(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        newData.push_back((val > 0) || (d[i] > 0));
    }

    return newData;
}
/*!
  Perform and XOR operation between data and \a d values.
  Returns the resulting DataSequence
  \note A value is \c true if is greather than 0
*/
DataSequence DataSequence::operator||(const DataSequence &d)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(i >= d.size ())
        {
            newData.push_back(1);
        }
        else
        {
            newData.push_back((this->at(i) > 0) ^ (d[i] > 0));
        }
    }

    return newData;
}
/*!
  Returns the fourier transform of the given data.
  If \a logflag is set to \c true the values will be normalized
  using logarithm.
*/
DataSequence DataSequence::fourierTransform(bool logflag) const
{
    if(size () == 0) throw HilbertBadOperation();

    double *datainput=NULL;
    fftw_complex *dataoutput=NULL;
    fftw_plan p;
    int data_size = size ();
    int i=0;
    int data_size2=data_size/2;
    DataSequence output;

    try
    {
        datainput = (double*) fftw_malloc(sizeof(double) * data_size);
        dataoutput = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (data_size2+2));
        if(datainput == NULL || dataoutput == NULL)
          throw std::bad_alloc(); // memory exhausted

        assert(datainput != NULL || dataoutput != NULL);

        p = fftw_plan_dft_r2c_1d(data_size, datainput, dataoutput, FFTW_ESTIMATE); // the transform

        for(auto val : (*this))
        { // filling the data
            assert(i < data_size);
            datainput[i++]=val;
        }

        fftw_execute(p); // the transform

        output.assign (data_size+1, 0);

    }catch (std::bad_alloc& ba)
    {
        throw HilbertBadAlloc();
    }

    // arranging transform data into redundant centered output
    for(i=1; i <= data_size2; ++i)
    {
        assert(data_size2+i <= data_size && data_size2+i >= 0);
        assert(data_size2-i <= data_size && data_size2-i >= 0);

        output[data_size2+i]=output[data_size2-i]=dataoutput[i][0]*dataoutput[i][0]+dataoutput[i][1]*dataoutput[i][1];
        if(logflag)
        {
               if(output[data_size2-i]>0)
                   output[data_size2-i]=output[data_size2+i]=log(sqrt(output[data_size2-i]));
        }
    }
    output[data_size2]=dataoutput[0][0]*dataoutput[0][0]+dataoutput[0][1]*dataoutput[0][1];

    if(logflag)
    {
           if(output[data_size2]>0)
               output[data_size2]=log(sqrt(output[data_size2]));
    }

    fftw_destroy_plan(p);
    fftw_free(datainput);
    fftw_free(dataoutput);

    output.pop_back();

    return output;
}
/*!
  Returns a new data wich elements are the Hamming distance between
  \c this and \a d elements.
*/
DataSequence DataSequence::hammingDistance(const DataSequence &d) const
{
    DataSequence hammingData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        hfloat val = 0;
        if(i < d.size ())
        {
            if(this->at(i) == d[i]) val = 1;
        }
        hammingData.push_back (val);
    }
    return hammingData;
}
/*!
  Returns a new data wich elements are the Manhattan distance between
  \c this and \a d elements.
*/
DataSequence DataSequence::manhattanDistance(const DataSequence &d) const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        hfloat val = 0;
        if(i < d.size ())
        {
            val = std::abs(this->at(i) - d[i]);
        }
        newData.push_back (val);
    }
    return newData;
}
/*!
  Returns a new data resulting from the elements that
  evaluates \c true in \a filterFunction.
*/
DataSequence DataSequence::filter(bool (*filterFunction)(hfloat))
{
    DataSequence newData;
    for(hfloat val: (*this))
    {
        if(filterFunction(val))
            newData.push_back (val);
    }
    return newData;
}
/*!
  Filter the data using \a filterFunction by performing comparsion
  element to element with \a other data. Returns a new data with the filtered elements
*/
DataSequence DataSequence::filterByComparsion(DataSequence &other, bool (*filterFunction)(const hfloat&,const hfloat&)) const
{
    DataSequence newData;
    unsigned int lenght = std::min(size (), other.size ());
    for(unsigned int i = 0; i < lenght; ++i)
    {
        if(filterFunction(this->at(i), other[i]))
        {
            newData.push_back (this->at(i));
        }
    }
    return newData;
}
/*!
  \overload filterByComparsion()
  Filter the data using \a filterFunction by performing comparsion
  between data elements and \a val. Returns a new data with the filtered elements
 */
DataSequence DataSequence::filterByComparsion(hfloat val, bool (*filterFunction)(const hfloat&,const hfloat&)) const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        if(filterFunction(this->at(i), val))
        {
            newData.push_back (this->at(i));
        }
    }
    return newData;
}
/*!
  Thresholds the values on the data according to the given \a thresholdFunction.
  Returns a new data resulting from thresholding.
*/
DataSequence DataSequence::thresholdData(hfloat (*thresholdFunction)(const hfloat &))
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back (thresholdFunction(this->at (i)));
    }
    return newData;
}
/*!
  \brief Compute granularity \a n on data.
*/
DataSequence &DataSequence::granularity(uint n)
{
    std::vector<hfloat>::const_iterator it = begin();
    std::vector<hfloat> result;
    hfloat sum=0;

    if(size() == 0)
    {
        std::cerr << "Data size is zero." << std::endl;
        return *this;
    }

    if(n > size()/2 || n == 0)
    {
        std::cerr << "Bad operation." << std::endl;
        return *this;
    }

    try
    {
        while(it < end()-(n-1))
        {
            sum=0;
            for(unsigned int i=0; i < n; ++i)
            {
                sum=sum+(*it);
                ++it;
            }
            sum=sum/hfloat(n);
            for(unsigned int i=0; i < n; ++i)
            {
                result.push_back(sum);
            }
        }
        while(it < end()) // last values (beyond size-gr) are not averaged
        {
            result.push_back(*it++);

        }
    }
    catch (std::bad_alloc& ba)
    {
        std::cerr << "Error applying granularity" << std::endl;
        throw HilbertBadAlloc();
    }
    (*this) = result;
    return *this;
}
/*!
  \brief Returns the minimum value in the data.
*/
hfloat DataSequence::max() const
{
    return *std::max_element(begin (), end ());
}
/*!
  \brief Returns the maximum value in the data.
*/
hfloat DataSequence::min() const
{
    return *std::min_element(begin (), end ());
}
/*!
 * \fn DataSequence::mean() const
 *  Mean value of the data
 */
hfloat DataSequence::mean() const
{
    if(size () == 0) return 0;
    hfloat sum = std::accumulate(begin (), end (), 0);
    return sum/hfloat(size ());
}
/*!
 * \brief DataSequence::stdDeviation
 *
 */
hfloat DataSequence::stdDeviation() const
{
    if(size () < 2) return 0;
    hfloat meanValue = mean ();
    hfloat sum = 0;
    for(const hfloat &val : (*this))
    {
        sum += (val - meanValue)*(val - meanValue);
    }
    return std::sqrt(1.0/hfloat(size ()-1)*sum);
}
/*!
 * \fn DataSequence::Entropy() const
 * \brief Compute the Shannon entropy of the data
 *  Shannon information entropy
 */
hfloat DataSequence::Entropy() const
{
    std::vector<unsigned long> freq;
    hfloat max=0;
    hfloat min=0;
    hfloat minmax=0;
    int index=0;

    if(size() == 0)
      throw HilbertBadSize();

    try
    {
        freq.assign (ENTROPY_LEVELS+3, 0);
    } catch  (std::bad_alloc& ba)
    {
        throw HilbertBadAlloc();
    }

    max = this->max ();
    min = this->min ();

    minmax = ENTROPY_LEVELS/(max-min);

    for(auto instance : (*this))
    {
        index=static_cast<unsigned int>(floor((instance-min)*minmax));
        assert(index < ENTROPY_LEVELS+1);
        freq[index]++;
    }

    hfloat val=0;
    int nbins = 0;
    for(auto instance : freq)
    {
        nbins += instance != 0;
        val=val+instance*zlog(instance);
    }
    nbins += nbins == 1;
    return (-val/size()+zlog(size()))/std::log(nbins);
}


/*!
  \brief Load a data in plain text format from \a input stream.
*/
DataSequence DataSequence::fromPlainText(std::istream &input)
{
    std::vector<hfloat> data;
    // Get the lenght of the stream
    input.seekg(0, input.end);
    int lenght = input.tellg();
    input.seekg(0, input.beg);
    // Buffer to store the file
    std::string buffer(lenght, '\0');
    input.read(&buffer[0], lenght); // Reading the whole file into buffer
    onlyNumbers (buffer);
    std::stringstream stream(buffer);
    hfloat val;
    while (stream >> val)
    {
        data.push_back (val);
    }
    return DataSequence(data);
}
/*!
  \overload fromPlainText()
  \brief Load a data in plain text format from \a input string.
*/
DataSequence DataSequence::fromPlainText(std::string &input)
{
    std::stringstream ss(input);
    return fromPlainText (ss);
}
/*!
  \brief Get rid of non numerical characters.
  Iterates throw \a input_string replacing non numerical character for space.
  \note As this function is used for extrac floating point values, characters using
  to represent floats like '.' , '-', '+' and 'e' are keeped.
*/
std::string DataSequence::onlyNumbers(std::string &input_string)
{
    for(std::string::size_type i = 0; i < input_string.size (); ++i)
    {
        if(!isNumeric(input_string.at (i)))
        {
            if(input_string.at (i) == 'e' && i > 0 && isNumeric (input_string.at (i-1)))
                continue;
            input_string[i] = ' ';
            continue;
        }
    }
    return  input_string;
}
/*!
  \brief Check if a character \a ch is a numeric value.
  Returns true if \a ch is a numeric value.
*/
bool DataSequence::isNumeric(char ch)
{
    return std::isdigit (ch)|| ch == '.' || ch == '-' || ch == '+';
}
/*!
  \overload operator||()
  Returns a DataSequence as result of XOR operation between data values
  and \a val.
  \note Positive values are taking as \c true
  \note Output values are \c 1 and \c 0.
*/
DataSequence DataSequence::operator ||(const hfloat &val)  const
{
    DataSequence newData;
    for(unsigned int i = 0; i < size (); ++i)
    {
        newData.push_back((this->at(i) > 0) ^ (val > 0));
    }

    return newData;
}
/*!
  \relates DataSequence
  Returns a DataSequence as result of applying XOR operator between
  \a val and \a d values.
*/
DataSequence operator||(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        newData.push_back((val > 0) ^ (d[i] > 0));
    }
    return newData;
}
/*!
  \relates DataSequence
  Returns a DataSequence as result of sum \a val with \a d elements.
*/
DataSequence operator+(const hfloat &val, const DataSequence &d)
{
    return d + val;
}
/*!
  \relates DataSequence
  Returns a DataSequence as result of substract \a val with \a d elements.
*/
DataSequence operator-(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        newData.push_back(val - d[i]);
    }
    return newData;
}
/*!
  \relates DataSequence
  Returns a DataSequence as result of multiplicate \a val with \a d elements.
*/
DataSequence operator*(const hfloat &val, const DataSequence &d)
{
    return d * val;
}
/*!
  \relates DataSequence
  Returns a DataSequence as result of divide \a val with \a d elements.
  If zero division detected HilbertZeroDivision exception will be raised.
*/
DataSequence operator/(const hfloat &val, const DataSequence &d)
{
    DataSequence newData;
    for(unsigned int i = 0; i < d.size (); ++i)
    {
        if(d[i] == 0)
        {
            throw HilbertZeroDivision();
        }
        newData.push_back(val / d[i]);
    }

    return newData;
}
/*!
  Takes logarithm of \a val giving zero if \a val is non positive.
*/
inline hfloat zlog(hfloat val)
{
    hfloat r=0;
    if (val > 0)
        r=log(val);
    return r;
}
/*!
  \brief Print the data.
*/
std::ostream &operator<<(std::ostream &out, DataSequence &data)
{
    for(auto val:data)
        out << val;
    return  out;
}
