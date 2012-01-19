
#pragma once

namespace MathLib
{
template<typename T, unsigned N>
class TemplateVector
{
public:
    T* getRaw() {return _data;};
private:
    T _data[N];
};

typedef TemplateVector<double, 2> Vector2D;
typedef TemplateVector<double, 3> Vector3D;
}
