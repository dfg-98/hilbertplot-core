#ifndef DATASEQUENCE_H
#define DATASEQUENCE_H

#include <vector>
#include <istream>
#include "hilbertdefines.h"

static const int ENTROPY_LEVELS = 65535;
static const hfloat LOG2_1 = 1.44269504088896340735992468100;
static const hfloat LOG2 = 0.693147180559945309417232121458;

//std::vector<hfloat> GiveMeDoubleVec(std::ifstream &in);

class DataSequence : public std::vector<hfloat>
{
    public:
        //Constructors
        DataSequence(void);
        DataSequence(int elements, hfloat defaultValue = 0.0);
        DataSequence(const std::vector<hfloat> &data);
        DataSequence(const DataSequence &data);
        DataSequence(DataSequence &&data);
        DataSequence(std::initializer_list<hfloat> data);

        DataSequence& operator= (const DataSequence& other);

        ///-----------------------------------------------------------
        ///  Operators overload
        /// ----------------------------------------------------------
        ///
        ///  Arithmetics:
        DataSequence operator+(const DataSequence & d) const;
        DataSequence operator+(const hfloat &value) const;
        friend DataSequence operator+(const hfloat &val,const  DataSequence &d);
        DataSequence operator-(const DataSequence & d) const;
        DataSequence operator-(const hfloat &val) const;
        friend DataSequence operator-(const hfloat &val,const  DataSequence &d);
        DataSequence operator*(const DataSequence & d) const;
        DataSequence operator*(const hfloat &value) const;
        friend DataSequence operator*(const hfloat &val,const  DataSequence &d);
        DataSequence operator/(const DataSequence & d) const;
        DataSequence operator/(const hfloat &val) const;
        friend DataSequence operator/(const hfloat &val,const  DataSequence &d);
        DataSequence operator^(const DataSequence & d) const;
        DataSequence operator^(const hfloat &val) const;
        friend DataSequence operator^(const hfloat &val,const  DataSequence &d);
        DataSequence operator==(const DataSequence & d) const;
        DataSequence operator==(const hfloat &val) const;
        friend DataSequence operator==(const hfloat &val,const  DataSequence &d);
        DataSequence operator!=(const DataSequence & d) const;
        DataSequence operator!=(const hfloat &val) const;
        friend DataSequence operator!=(const hfloat &val,const  DataSequence &d);
        DataSequence operator>(const DataSequence & d) const;
        DataSequence operator>(const hfloat &val) const;
        friend DataSequence operator>(const hfloat &val,const  DataSequence &d);
        DataSequence operator<(const DataSequence & d) const;
        DataSequence operator<(const hfloat &val) const;
        friend DataSequence operator<(const hfloat &val,const  DataSequence &d);
        DataSequence operator>=(const DataSequence & d) const;
        DataSequence operator>=(const hfloat &val) const;
        friend DataSequence operator>=(const hfloat &val,const  DataSequence &d);
        DataSequence operator<=(const DataSequence & d) const;
        DataSequence operator<=(const hfloat &val) const;
        friend DataSequence operator<=(const hfloat &val,const  DataSequence &d);
        DataSequence operator&&(const DataSequence & d) const;
        DataSequence operator&&(const hfloat &val) const;
        friend DataSequence operator&&(const hfloat &val,const  DataSequence &d);
        DataSequence operator|(const DataSequence & d) const;
        DataSequence operator|(const hfloat &val) const;
        friend DataSequence operator|(const hfloat &val,const  DataSequence &d);
        DataSequence operator||(const DataSequence & d) const;
        DataSequence operator||(const hfloat &val) const;
        friend DataSequence operator||(const hfloat &val,const  DataSequence &d);

        // Fourier
        DataSequence fourierTransform(bool logflag) const;

        // Distances
        DataSequence hammingDistance(const DataSequence &d) const;
        DataSequence manhattanDistance(const DataSequence &d) const;

        /// FILTERING
        DataSequence filter(bool (*filterFunction)(hfloat));
        DataSequence filterByComparsion(DataSequence &other, bool (*filterFunction)(const hfloat &, const hfloat &)) const;
        DataSequence filterByComparsion(hfloat val, bool (*filterFunction)(const hfloat &, const hfloat &)) const;

        /// Thresholidng
        DataSequence thresholdData(hfloat (*thresholdFunction)(const hfloat&));

        DataSequence &granularity(uint n);

        /// Data information:
        hfloat max() const;
        hfloat min() const;
        hfloat mean() const;
        hfloat stdDeviation() const;
        hfloat Entropy() const;

        friend std::ostream & operator<<(std::ostream & out, DataSequence & data);

        static DataSequence fromPlainText(std::istream &input);
        static DataSequence fromPlainText(std::string &input);
        static std::string onlyNumbers(std::string &input_string);
        static bool isNumeric(char ch);
};
DataSequence operator+(const hfloat &val,const  DataSequence &d);
DataSequence operator-(const hfloat &val,const  DataSequence &d);
DataSequence operator+(const hfloat &val,const  DataSequence &d);
DataSequence operator*(const hfloat &val,const  DataSequence &d);
DataSequence operator^(const hfloat &val,const  DataSequence &d);
DataSequence operator==(const hfloat &val,const  DataSequence &d);
DataSequence operator!=(const hfloat &val,const  DataSequence &d);
DataSequence operator>(const hfloat &val,const  DataSequence &d);
DataSequence operator<(const hfloat &val,const  DataSequence &d);
DataSequence operator>=(const hfloat &val,const  DataSequence &d);
DataSequence operator<=(const hfloat &val,const  DataSequence &d);
DataSequence operator&&(const hfloat &val,const  DataSequence &d);
DataSequence operator|(const hfloat &val,const  DataSequence &d);
DataSequence operator||(const hfloat &val,const  DataSequence &d);
std::ostream & operator<<(std::ostream & out, DataSequence & data);
#endif // DATA_H

