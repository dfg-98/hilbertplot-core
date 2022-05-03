#ifndef HILBERTPLOT_H
#define HILBERTPLOT_H

#include "hilbertcurve.h"
#include "datasequence.h"
#include <utility>

#ifdef QT_CORE
    #include <QDataStream>
#endif

class HilbertPlot : public HilbertCurve
{
    public:


        HilbertPlot();
        HilbertPlot(const DataSequence &data, hsize width = 0, hsize height = 0, CurveType type = H0);
        HilbertPlot(const HilbertPlot &hilbertplot);

        std::vector<HPoint>::const_reference operator [] (std::vector<HPoint>::size_type index) const;
        std::vector<HPoint>::const_reference at (std::vector<HPoint>::size_type index);
        std::vector<HPoint>::const_reference at (std::vector<HPoint>::size_type x, std::vector<HPoint>::size_type y);
        hfloat valueAt(std::vector<hfloat>::size_type index) const;
        hfloat valueAt(std::vector<hfloat>::size_type x, std::vector<hfloat>::size_type y) const;
        hfloat valueNormalizedAt(std::vector<hfloat>::size_type index) const;
        hfloat valueNormalizedAt(std::vector<hfloat>::size_type x, std::vector<hfloat>::size_type y) const;
        void replaceValueAt(std::vector<hfloat>::size_type index, hfloat value);
        void replaceValueAt(std::vector<hfloat>::size_type x, std::vector<hfloat>::size_type y, hfloat value);
        hint indexOf(hint x, hint y) const;

        hfloat min() const {return  m_min;}
        hfloat max() const {return  m_max;}

        HImage generateImage(hfloat threshold = 0);

        DataSequence dataCopy() const;
        void replaceData(const DataSequence &data);

        DataSequence hpFourierTransform(bool logflag) const;
        static std::pair<hsize, hsize> bestDimensions(hsize lenght);

    private:
        DataSequence m_data;
        hfloat m_min;
        hfloat m_max;
        std::vector<std::vector<hint>> m_plotToCurve;
        static const HilbertCurve constructCurve(hsize lenght, hsize &width, hsize &height, CurveType type);
};
#endif // HILBERTPLOT_H
