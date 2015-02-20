#ifndef AOUTOUT_HPP
#define AOUTOUT_HPP

#include <vector>
#include <string>

template< typename XType, typename YType >
class AOut {
public:
    AOut(std::size_t coordCount) : _x(coordCount), _y(coordCount) { }

    void addCoordinate(XType x, YType y) {
        _x.push_back(x);
        _y.push_back(y);
    }

    virtual void out(const std::string& name) = 0;

protected:
    std::vector< XType > _x;
    std::vector< YType > _y;
};

typedef AOut< float, float > f_AOut;

#endif // AOUTOUT_HPP
