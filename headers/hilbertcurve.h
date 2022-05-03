#ifndef HILBERTCURVE_H
#define HILBERTCURVE_H

#include <vector>

#include "hpoint.h"
#include "hilbertdefines.h"

class QuasiSquare
{
    public:

        enum Orientation {A, B, C, D};

        //class constructors
        QuasiSquare(const QuasiSquare & q);
        QuasiSquare();
        QuasiSquare(hsize nn, hsize mm, HPoint point, Orientation o);

        //operator overloading
        QuasiSquare & operator=(const QuasiSquare & qs);
        //Build curve functions
        std::vector<HPoint> &BuildCurve(std::vector<HPoint> & coordinates_list, hsize index);
        friend class HilbertCurve;
    protected:
        hsize n;
        hsize m;
        HPoint coord; // Origen
        Orientation oABCD;
        //Partition function
        std::vector<QuasiSquare> & Partition(std::vector<QuasiSquare> & partition_vec);

};

class HilbertCurve : protected QuasiSquare
{
    public:
        typedef typename std::vector<HPoint>::iterator iterator;
        typedef typename std::vector<HPoint>::const_iterator const_iterator;
        enum CurveType{H0,  H1,  H2,  H3,  H4,  H5,  H6,  H7,  H8,  H9,
                       H10, H11, H12, H13, H14, H15, H16, H17, H18, H19,
                       H20, H21, H22, H23, H24, H25, H26, H27, H28, H29,
                       H30, H31, H32, H33, H34, H35, H36, H37, H38, H39};

        HilbertCurve(void);
        HilbertCurve(hsize width, hsize height, CurveType type = H0, HPoint origen = 0, Orientation orientation = A, bool differenceCurve = false);
        HilbertCurve(const HilbertCurve &curve);

        hfloat meanDifference() const;
        hsize lenght() const;
        hsize width() const;
        hsize height() const;
        CurveType type() const;

        std::vector<HPoint>::reference operator[](std::vector<HPoint>::size_type index);
        std::vector<HPoint>::const_reference operator[](std::vector<HPoint>::size_type index) const;
        const_iterator begin() const;
        const_iterator cbegin() const;
        const_iterator end() const;
        const_iterator cend() const;
        void SaveSVG(const char* filename, const char* colorName, float stroke_width = 0.2);
        std::string curveToSvg(const char* colorName = "red", float stroke_width = 0.2);

        static HilbertCurve createCurve(hsize width, hsize height, CurveType type = H0, HPoint origen = 0, Orientation orientation = A, bool differenceCurve = false);
        friend class HilbertPlotForm;

    private:
        CurveType m_type;
        std::vector<HPoint> m_curve;
        hfloat m_mean_difference;

        void Build();
        void BuildDifference();
        void BuildCurveH0();
        void BuildCurve1H();
        void BuildCurve2H();
        void BuildCurve3H();
        void BuildCurve4H();
        void BuildCurve5H();
        void BuildCurve6H();
        void BuildCurve7H();
        void BuildCurve8H();
        void BuildCurve9H();
        void BuildCurve10H();
        void BuildCurve11H();
        void BuildCurve12H();
        void BuildCurve13H();
        void BuildCurve14H();
        void BuildCurve15H();
        void BuildCurve16H();
        void BuildCurve17H();
        void BuildCurve18H();
        void BuildCurve19H();
        void BuildCurve20H();
        void BuildCurve21H();
        void BuildCurve22H();
        void BuildCurve23H();
        void BuildCurve24H();
        void BuildCurve25H();
        void BuildCurve26H();
        void BuildCurve27H();
        void BuildCurve28H();
        void BuildCurve29H();
        void BuildCurve30H();
        void BuildCurve31H();
        void BuildCurve32H();
        void BuildCurve33H();
        void BuildCurve34H();
        void BuildCurve35H();
        void BuildCurve36H();
        void BuildCurve37H();
        void BuildCurve38H();
        void BuildCurve39H();

        void joinCurve(HilbertCurve &cuadrant1, HilbertCurve &cuadrant2, HilbertCurve &cuadrant3, HilbertCurve &cuadrant4);
        void reflectAndReverse();
        void reverse();
        void reflect();
        void reflectX();
        void reflectY ();
};

#endif // HILBERTCURVE_H
