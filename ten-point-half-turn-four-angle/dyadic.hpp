#pragma once

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <stdexcept>

using Integer = std::int64_t;
using Wide = __int128_t;
constexpr Integer unit = 1LL << 43;
constexpr Integer limit = 1LL << 60;

inline Integer checked(Wide v) {
    if (v < -limit || v > limit) throw std::overflow_error("dyadic range");
    return static_cast<Integer>(v);
}

inline Integer floor_ratio(Wide a, Wide b) {
    if (b <= 0) throw std::domain_error("nonpositive divisor");
    return checked(a/b-(a%b<0));
}

inline Integer ceil_ratio(Wide a, Wide b) { return -floor_ratio(-a,b); }

struct Interval { Integer lo=0, hi=0; };
inline Interval rat(Integer a, Integer b=1) {
    return {floor_ratio(Wide(a)*unit,b),ceil_ratio(Wide(a)*unit,b)};
}

inline Interval operator+(Interval a, Interval b) {
    return {checked(Wide(a.lo)+b.lo),checked(Wide(a.hi)+b.hi)};
}
inline Interval operator-(Interval a) { return {-a.hi,-a.lo}; }
inline Interval operator-(Interval a, Interval b) { return a+(-b); }
inline Interval operator*(Interval a, Interval b) {
    auto [l,h]=std::minmax({Wide(a.lo)*b.lo,Wide(a.lo)*b.hi,Wide(a.hi)*b.lo,Wide(a.hi)*b.hi});
    return {floor_ratio(l,unit),ceil_ratio(h,unit)};
}
inline Interval operator*(Interval a, int k) {
    if (k<0) return (-a)*(-k);
    return {checked(Wide(a.lo)*k),checked(Wide(a.hi)*k)};
}
inline Interval operator/(Interval a, Integer k) {
    return {floor_ratio(a.lo,k),ceil_ratio(a.hi,k)};
}
inline Interval operator/(Interval a, Interval b) {
    if (b.lo<=0) throw std::domain_error("nonpositive interval divisor");
    return {floor_ratio(Wide(a.lo)*unit,a.lo<0?b.lo:b.hi),
            ceil_ratio(Wide(a.hi)*unit,a.hi>0?b.lo:b.hi)};
}
inline Integer mag(Interval a) { return std::max(std::abs(a.lo),std::abs(a.hi)); }
inline Interval square(Interval a) {
    Integer low=a.lo<=0 && a.hi>=0 ? 0 : std::min(std::abs(a.lo),std::abs(a.hi));
    return {floor_ratio(Wide(low)*low,unit),ceil_ratio(Wide(mag(a))*mag(a),unit)};
}
inline Integer sqrt_upper(Integer x) {
    if (x<0) throw std::domain_error("negative radicand");
    Wide n=Wide(x)*unit;
    auto root=static_cast<Integer>(std::sqrt(static_cast<double>(n)));
    while (Wide(root)*root<n) ++root;
    while (root>0 && Wide(root-1)*(root-1)>=n) --root;
    return root;
}
inline Interval point(Integer a) { return {a,a}; }
inline Interval center(Interval a) { return point((a.lo+a.hi)/2); }
inline Interval outer(Integer a, Integer b, Integer c, Integer d) { return {rat(a,b).lo,rat(c,d).hi}; }
inline bool nonzero(Interval a) { return a.lo>0 || a.hi<0; }
inline bool meet(Interval &a, Interval b) {
    a={std::max(a.lo,b.lo),std::min(a.hi,b.hi)};
    return a.lo<=a.hi;
}
inline double midpoint(Interval a) { return double(a.lo+a.hi)/(2*unit); }
