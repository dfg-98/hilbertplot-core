#ifndef HILBERTDEFINES_H
#define HILBERTDEFINES_H

#include <exception>
#include <vector>

typedef unsigned int hint;
typedef unsigned int hsize;
typedef double hfloat;

#define vecvecfloat std::vector<std::vector<hfloat>>

typedef vecvecfloat HImage;

class HilbertBadAlloc : public std::bad_alloc{};
class HilbertBadOrientation : public std::exception{};
class HilbertBadOperation : public std::exception{};
class HilbertIndexOutOfRange : public std::exception{};
class HilbertZeroDivision : public std::exception{};
class HilbertBadSize : public std::exception{};

#endif // HILBERTDEFINES_H
