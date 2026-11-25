// Four-angle cover. All conclusions use the integer interval operations.
#include "dyadic.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <utility>
#include <vector>

using Angles = std::array<Interval,4>;
using Radii = std::array<Interval,5>;
using ForceAllowance = std::array<Integer,5>;
using Vec = std::array<Interval,9>;
using Matrix = std::array<Vec,9>;
using FloatMatrix = std::array<std::array<double,9>,9>;

const Interval one=rat(1), tau=outer(62831853071795LL,10000000000000LL,
                                     62831853071797LL,10000000000000LL);
const Interval lower_gap=rat(43,200), upper_gap=rat(3,10);
const Interval gap_cos=outer(-1519,5000,271,1250);
const Interval inner_band=outer(43,400,53,400), outer_band=outer(89,100,1271,1000);
const Interval mass=rat(9,2);

struct Direction { Interval cos, sin; };
struct Edge {
    int i, j;
    Direction center, range;
};
using Edges = std::array<Edge,10>;
using Turn = std::array<Direction,4>;
std::array<Turn,2> anchors;
std::array<Radii,2> anchor_radii;
Integer nodes=0, geometry=0, energy=0, force=0, krawczyk=0, basin=0, saddle=0, unknown=0;
Integer radial_inclusions=0, split=0, zero_midpoints=0;
Integer energy_radial_localizations=0, local_hessian_refinements=0, sharper_radial_inclusions=0;
Integer radial_subcovers=0, radial_subcover_nodes=0;
int deepest=0, max_depth=0;
bool small_midpoints=false;
Interval force_slack, energy_slack, projection_slack, midpoint_sq, individual_midpoint_sq, midpoint_test;

bool outside_force_slack(Interval x) {
    return small_midpoints ? x.lo>force_slack.hi || x.hi<-force_slack.hi : nonzero(x);
}

Direction direction(Interval angle) {
    // Reduce the exact midpoint to within one eighth turn of a quadrant.
    Interval t=center(angle);
    Integer q=floor_ratio(Wide(t.lo)*4+unit/2,unit);
    Interval x=(t-point(q*(unit/4)))*tau, y=square(x);
    Interval cosine=one, sine=one, ec=one, es=one;
    for (int k=1;k<=7;++k) {
        ec=(-ec*y)/(2*k*(2*k-1));
        es=(-es*y)/(2*k*(2*k+1));
        cosine=cosine+ec;
        sine=sine+es;
    }
    if (y.hi>=unit) throw std::logic_error("unreduced trigonometric argument");
    cosine=cosine+Interval{0,mag(ec)/240+1};
    sine=x*(sine+Interval{0,mag(es)/272+1});
    switch ((q%4+4)%4) {
      case 1: {auto old=cosine; cosine=-sine; sine=old; break;}
      case 2: cosine=-cosine; sine=-sine; break;
      case 3: {auto old=cosine; cosine=sine; sine=-old; break;}
    }
    const Interval delta=(angle-t)*tau;
    const Interval small_cos{(one-square(delta)/2).lo,unit};
    const Interval small_sin{-mag(delta),mag(delta)};
    Direction result{cosine*small_cos-sine*small_sin, sine*small_cos+cosine*small_sin};
    if (!meet(result.cos,{-unit,unit}) || !meet(result.sin,{-unit,unit}))
        throw std::logic_error("unit circle enclosure");
    return result;
}

bool inside_bands(const Radii &r) {
    if (r[0].lo<rat(43,400).hi || r[0].hi>rat(53,400).lo) return false;
    for (int i=1; i<5; ++i)
        if (r[i].lo<rat(89,100).hi || r[i].hi>rat(1271,1000).lo) return false;
    return true;
}

bool inside_log_domain(const Radii &r) {
    if (!inside_bands(r)) return false;
    Interval product=one, min_product=one;
    for (int i=1; i<5; ++i) {
        product=product*r[i];
        min_product=min_product*rat(43,40);
    }
    return product.lo>=min_product.hi;
}

bool gap_constraints(Angles &a) {
    Interval allowed{lower_gap.lo,upper_gap.hi};
    for (int repeat=0;repeat<2;++repeat) {
        for (int i=0;i<4;++i) {
            int j=(i+1)%4;
            Interval wrap=j==0?one:Interval{};
            if (!meet(a[j],a[i]+allowed-wrap) || !meet(a[i],a[j]-allowed+wrap)) return false;
        }
    }
    return true;
}

Edges make_edges(const Angles &a) {
    std::array<Interval,5> t{}, t0{};
    for (int i=0;i<4;++i) {
        t[i+1]=a[i];
        t0[i+1]=center(a[i]);
    }
    Edges edges;
    int k=0;
    for (int i=0;i<5;++i)
        for (int j=i+1;j<5;++j)
            edges[k++]={i,j,direction(t0[i]-t0[j]),direction(t[i]-t[j])};
    return edges;
}

bool geometry_ok(const Edges &edges) {
    for (auto e:edges) if (e.i!=0) {
        const bool adjacent=e.j-e.i!=2;
        if (!meet(e.range.cos,adjacent?gap_cos:outer(-1,1,-81,100))) return false;
    }
    return true;
}

int orbit(const Angles &angle) {
    Turn d;
    for (int j=0;j<4;++j) d[j]=direction(angle[j]);
    for (int k=0;k<2;++k) {
        Interval e=k==0?rat(11,500):rat(3,125);
        Interval chord_lower=square(e)-square(square(e))/12;
        for (int shift=0;shift<4;++shift) {
            bool inside=true;
            for (int j=0;j<4;++j) {
                auto a=anchors[k][(j+shift)%4];
                auto chord=square(a.cos-d[j].cos)+square(a.sin-d[j].sin);
                inside &= chord.hi<chord_lower.lo;
            }
            if (inside) return k+1;
        }
    }
    return 0;
}

Turn angle_directions(const Angles &angle) {
    Turn d;
    for (int j=0;j<4;++j) d[j]=direction(angle[j]);
    return d;
}

int paired_orbit(const Turn &d,const Radii &r,Interval repair={}) {
    for (int k=0;k<2;++k) {
        auto diff=square((r[0]-anchor_radii[k][0])/rat(13,20));
        for (int shift=0;shift<4;++shift) {
            auto total=diff;
            for (int j=0;j<4;++j) {
                int u=(j+shift)%4;
                auto a=anchors[k][u];
                auto chord=square(a.cos-d[j].cos)+square(a.sin-d[j].sin);
                total=total+square((r[j+1]-anchor_radii[k][u+1])/rat(9,5))
                      + chord*rat(1271,2000);
            }
            const auto amplitude=point(sqrt_upper(total.hi))+repair;
            if (k==0 && amplitude.hi<rat(7,250).lo) return k+1;
            if (k==1 && ((square(amplitude)+midpoint_sq)*2).hi<square(rat(21,500)).lo)
                return k+1;
        }
    }
    return 0;
}

struct ComplexInterval {Interval re,im;};
ComplexInterval operator*(ComplexInterval a, ComplexInterval b) {
    return {a.re*b.re-a.im*b.im,a.re*b.im+a.im*b.re};
}
ComplexInterval operator-(ComplexInterval a, ComplexInterval b) {
    return {a.re-b.re,a.im-b.im};
}
ComplexInterval scale(ComplexInterval a,Interval b) { return {a.re*b,a.im*b}; }
Integer magnitude(ComplexInterval a) { return sqrt_upper((square(a.re)+square(a.im)).hi); }

Angles unscaled_angle_box(Angles a) {
    const auto shift=(individual_midpoint_sq*6/5)/tau;
    for (auto &v:a) v=v+Interval{-shift.hi,shift.hi};
    return a;
}

Radii unscaled_radius_box(Radii r) {
    for (int i=1;i<5;++i) r[i]=r[i]-Interval{0,individual_midpoint_sq.hi};
    return r;
}

using RealQuadratic=std::array<std::array<Integer,10>,10>;
template<std::size_t N>
bool positive_definite(std::array<std::array<Integer,N>,N> a) {
    for (std::size_t i=0;i<N;++i) {
        if (a[i][i]<=0) return false;
        for (std::size_t j=i+1;j<N;++j) {
            if (a[j][j]<=0 || Wide(a[i][j])*a[i][j]>=Wide(a[i][i])*a[j][j]) return false;
        }
        for (std::size_t j=i+1;j<N;++j) {
            for (std::size_t k=j;k<N;++k)
                a[j][k]=checked(Wide(a[j][k])-ceil_ratio(Wide(a[i][j])*a[i][k],a[i][i]));
            a[j][j]-=static_cast<Integer>(N-1-i);
        }
    }
    return true;
}

ComplexInterval midpoint_weight(Interval ri,Interval rj,Direction u,Direction relative,Interval d) {
    auto p=rj*relative.cos, q=rj*relative.sin;
    ComplexInterval numerator{ri+p,q}, conjugate_difference{ri-p,q};
    auto weight=numerator*conjugate_difference*conjugate_difference*ComplexInterval{u.cos,-u.sin};
    const auto den=square(d);
    return {weight.re*2/den,weight.im*2/den};
}

Integer rectangular_error(ComplexInterval weight) {
    return magnitude({weight.re-center(weight.re),weight.im-center(weight.im)});
}

struct MidpointAngles {
    std::array<Direction,5> box{}, center{};
    std::array<Direction,10> relative{}, relative_center{};
    Radii movement{};
};

MidpointAngles midpoint_angles(const Angles &raw) {
    std::array<Interval,5> a{};
    for (int i=1;i<5;++i) a[i]=raw[i-1];
    MidpointAngles prepared;
    for (int i=0;i<5;++i) {
        prepared.box[i]=direction(a[i]);
        prepared.center[i]=direction(center(a[i]));
        prepared.movement[i]=point(mag((a[i]-center(a[i]))*tau));
    }
    int edge=0;
    for (int i=0;i<5;++i) for (int j=i+1;j<5;++j,++edge) {
        prepared.relative[edge]=direction(a[j]-a[i]);
        prepared.relative_center[edge]=direction(center(a[j])-center(a[i]));
    }
    return prepared;
}

bool requires_zero_midpoints(const MidpointAngles &a,const Radii &r,Interval threshold) {
    Radii radius_center{}, movement{};
    for (int i=0;i<5;++i) {
        radius_center[i]=center(r[i]);
        movement[i]=point(mag(r[i]-radius_center[i]))+r[i]*a.movement[i];
    }
    RealQuadratic lower{};
    for (int i=0;i<10;++i) lower[i][i]=10*unit-threshold.hi;
    int edge=0;
    for (int i=0;i<5;++i) for (int j=i+1;j<5;++j,++edge) {
        const auto relative=a.relative[edge];
        auto d=square(r[i]-r[j])+r[i]*r[j]*(one-relative.cos)*2;
        if (d.lo<=0) return false;
        auto weight=midpoint_weight(r[i],r[j],a.box[i],relative,d);
        auto error=rectangular_error(weight);
        const auto sqrt_d=point(sqrt_upper(d.lo)-1);
        if (sqrt_d.lo>0) {
            const auto cross=r[i]*r[j]*relative.cos*6;
            const auto first=square(r[i])+square(r[j])*9+cross;
            const auto second=square(r[i])*9+square(r[j])+cross;
            const auto slope=point(sqrt_upper(first.hi))*movement[i]+point(sqrt_upper(second.hi))*movement[j];
            const auto change=slope*2/(point(d.lo)*sqrt_d);
            if (change.hi<error) {
                const auto rc=radius_center[i], sc=radius_center[j];
                const auto rel=a.relative_center[edge];
                const auto dist=square(rc-sc)+rc*sc*(one-rel.cos)*2;
                if (dist.lo>0) {
                    const auto value=midpoint_weight(rc,sc,a.center[i],rel,dist);
                    const auto centered_error=change.hi+rectangular_error(value);
                    if (centered_error<error) { weight=value; error=centered_error; }
                }
            }
        }
        auto re=center(weight.re), im=center(weight.im);
        std::array<std::array<Integer,2>,2> block{{
            {2*(re.lo-error),-2*im.lo},
            {-2*im.lo,2*(-re.lo-error)}
        }};
        for (int x=0;x<2;++x) for (int y=0;y<2;++y) {
            lower[2*i+x][2*i+y]+=block[x][y];
            lower[2*j+x][2*j+y]+=block[x][y];
            lower[2*i+x][2*j+y]-=block[x][y];
            lower[2*j+x][2*i+y]-=block[x][y];
        }
    }
    return positive_definite(lower);
}

struct Derivatives {
    Vec g{};
    std::array<Interval,10> b{},c{},distance{};
};

bool derivatives(const Radii &r, const Edges &e, bool at_center, Derivatives &out, bool hessian=true) {
    for (int i=0;i<5;++i) out.g[i]=rat(9)-r[i]*10;
    for (int k=0;k<10;++k) {
        const auto [i,j,mid,full]=e[k];
        const auto [cos,sin]=at_center?mid:full;
        Interval p=r[i]*r[j]*2, base=square(r[i])+square(r[j]);
        Interval difference=(r[i]-r[j])*(r[i]+r[j]);
        Interval d=square(r[i]-r[j])+p*(one-cos);
        if (d.lo<=0) return false;
        Interval radial=difference*2/d, angular=p*sin*2/d;
        out.g[i]=out.g[i]+radial;
        out.g[j]=out.g[j]-radial;
        if (i) out.g[i+4]=out.g[i+4]+angular;
        out.g[j+4]=out.g[j+4]-angular;
        out.distance[k]=d;
        if (hessian) {
            auto den=square(d);
            out.b[k]=p*(p-base*cos)*2/den;
            out.c[k]=-(p*difference*sin*2/den);
        }
    }
    return true;
}

// Left-multiply the derivative by P while preserving the common edge weights.
Matrix transformed(const Matrix &p, int rows, const Derivatives &d, const Radii &r,
                   const Edges &edges, bool radial=true, bool angular=true) {
    Matrix out{};
    for (int a=0;a<rows;++a) {
        if (radial) for (int i=0;i<5;++i) out[a][i]=-(p[a][i]*(r[i]*10));
        for (int k=0; k<10; ++k) {
            const int i=edges[k].i,j=edges[k].j;
            const auto u=p[a][i]-p[a][j];
            const auto v=i?p[a][i+4]-p[a][j+4]:-p[a][j+4];
            if (radial) {
                Interval value=u*d.b[k]+v*d.c[k];
                out[a][i]=out[a][i]+value;
                out[a][j]=out[a][j]-value;
            }
            if (angular) {
                Interval value=u*d.c[k]-v*d.b[k];
                if (i) out[a][i+4]=out[a][i+4]+value;
                out[a][j+4]=out[a][j+4]-value;
            }
        }
        if (radial) for (int i=0;i<5;++i) out[a][i]=out[a][i]/r[i];
        if (angular) for (int i=5;i<9;++i) out[a][i]=out[a][i]*tau;
    }
    return out;
}

bool inverse(FloatMatrix a, int n, FloatMatrix &out) {
    out={};
    for (int i=0;i<n;++i) out[i][i]=1;
    for (int j=0;j<n;++j) {
        int b=j;
        for (int k=j+1;k<n;++k) if (std::abs(a[k][j])>std::abs(a[b][j])) b=k;
        std::swap(a[j],a[b]); std::swap(out[j],out[b]);
        double p=a[j][j];
        if (std::abs(p)<1e-8) return false;
        for (int k=0;k<n;++k) { a[j][k]/=p; out[j][k]/=p; }
        for (int i=0;i<n;++i) if (i!=j) {
            const double scale=a[i][j];
            for (int k=0;k<n;++k) { a[i][k]-=scale*a[j][k]; out[i][k]-=scale*out[j][k]; }
        }
    }
    return true;
}

FloatMatrix float_jacobian(const Derivatives &d, const Radii &r, const Edges &e) {
    FloatMatrix a{};
    for (int i=0;i<5;++i) a[i][i]=-10*midpoint(r[i]);
    for (int k=0;k<10;++k) {
        int i=e[k].i,j=e[k].j;
        double b=midpoint(d.b[k]), c=midpoint(d.c[k]);
        a[i][i]+=b; a[j][j]+=b; a[i][j]-=b; a[j][i]-=b;
        a[j][j+4]+=c; a[i][j+4]-=c; a[j+4][j]+=c; a[j+4][i]-=c;
        a[j+4][j+4]-=b;
        if (i) {
            a[i][i+4]+=c; a[j][i+4]-=c; a[i+4][i]+=c; a[i+4][j]-=c;
            a[i+4][i+4]-=b; a[i+4][j+4]+=b; a[j+4][i+4]+=b;
        }
    }
    for (auto &row:a) {
        for (int i=0;i<5;++i) row[i]/=midpoint(r[i]);
        for (int i=5;i<9;++i) row[i]*=midpoint(tau);
    }
    return a;
}

bool propose_matrix(const FloatMatrix &input,int n,Matrix &out, bool invert=true) {
    FloatMatrix a=input;
    if (invert && !inverse(input,n,a)) return false;
    for (int i=0;i<n;++i) for (int j=0;j<n;++j) {
        if (!std::isfinite(a[i][j]) || std::abs(a[i][j])>128) return false;
        out[i][j]=point(std::llround(a[i][j]*unit));
    }
    return true;
}

Radii reference(const Edges &e, const Radii &seed) {
    // Floating Newton iterations propose a dyadic point; feasibility is checked exactly.
    std::array<double,5> r;
    for (int i=0;i<5;++i) r[i]=midpoint(seed[i]);
    std::array<double,10> c;
    for (int k=0;k<10;++k) c[k]=midpoint(e[k].center.cos);
    for (int step=0;step<8;++step) {
        FloatMatrix a{}, inv;
        std::array<double,5> g;
        for (int i=0;i<5;++i) { g[i]=9-10*r[i]; a[i][i]=-10*r[i]; }
        bool safe=true;
        for (int k=0;k<10;++k) {
            int i=e[k].i,j=e[k].j;
            double p=2*r[i]*r[j], base=r[i]*r[i]+r[j]*r[j], d=base-p*c[k];
            if (d<=.01) { safe=false; break; }
            double f=2*(r[i]*r[i]-r[j]*r[j])/d, b=2*p*(p-base*c[k])/(d*d);
            g[i]+=f; g[j]-=f; a[i][i]+=b; a[j][j]+=b; a[i][j]-=b; a[j][i]-=b;
        }
        if (!safe || !inverse(a,5,inv)) break;
        std::array<double,5> next=r;
        double change=0;
        for (int i=0;i<5;++i) {
            double ds=0;
            for (int j=0;j<5;++j) ds-=inv[i][j]*g[j];
            change=std::max(change,std::abs(ds));
            next[i]=r[i]*(1+std::clamp(ds,-.2,.2));
            safe &= next[i]>.06 && next[i]<1.7;
        }
        if (!safe) break;
        r=next;
        if (change<1e-11) break;
    }
    double inner=std::clamp(r[0],.1075000001,.1324999999);
    double redistribute=(4.5-inner-r[1]-r[2]-r[3]-r[4])/4;
    for (int i=1;i<5;++i) r[i]+=redistribute;
    r[0]=inner;
    Radii result;
    Integer sum=0;
    for (int i=0;i<4;++i) {
        result[i]=point(std::llround(r[i]*unit));
        sum+=result[i].lo;
    }
    result[4]=point(mass.lo-sum);
    return inside_log_domain(result)?result:seed;
}

Interval exp_upper(Interval x) {
    Interval sum=one, term=one;
    for (int k=1;k<=6;++k) { term=(term*x)/k; sum=sum+term; }
    return sum+(term*x/7)/(one-x/8);
}

const Interval product_witness=rat(469765409,1LL<<26);
struct EnergyData { Interval product,angular,radial; };

EnergyData energy_data(const Derivatives &mid,const Derivatives &fixed,const Edges &edges,
                       const Radii &ref,const Angles &a,Interval slack) {
    Vec disp{};
    for (int i=0;i<4;++i) disp[i+5]=(a[i]-center(a[i]))*tau;
    Interval increase=slack;
    Interval radial;
    for (int i=0;i<5;++i) radial=radial+square(fixed.g[i])*2;
    for (int i=5;i<9;++i) increase=increase+mid.g[i]*disp[i];
    Interval product=one;
    for (auto r:ref) product=product*r;
    for (int k=0;k<10;++k) {
        const auto [i,j,theta0,theta]=edges[k];
        Interval p=ref[i]*ref[j]*2;
        if ((p-(square(ref[i])+square(ref[j]))*theta.cos).lo<0) {
            Interval err=disp[j+4]-(i?disp[i+4]:Interval{});
            increase=increase+p*square(err)/mid.distance[k];
        }
        if (k==4) product=product/(1<<16);
        product=product*square(mid.distance[k]);
    }
    return {product,increase,radial};
}

bool energy_excluded(const EnergyData &data) {
    Interval bound=point(std::max<Integer>(0,(data.angular+data.radial).hi));
    return bound.hi<=unit && (data.product*exp_upper(bound)).hi<product_witness.lo;
}

// Complete the radial concavity square at the feasible reference. A large
// objective confines the repaired radii even without an exact radial root.
bool energy_radii(const EnergyData &data,const Radii &ref,const Derivatives &fixed,Radii &r) {
    const auto budget=data.product/product_witness-one+data.angular+data.radial;
    if (budget.hi>unit/16) return false;
    const auto radius=sqrt_upper(8*std::max<Integer>(0,budget.hi));
    Integer total=0,max_total=0;
    for (int i=0;i<5;++i) {
        const Interval band=i==0?inner_band:outer_band;
        const auto delta=(fixed.g[i]*4+Interval{-radius,radius})*band;
        r[i]=ref[i]+delta;
        if (!meet(r[i],band)) return false;
        total+=r[i].hi-r[i].lo;
        max_total+=band.hi-band.lo;
    }
    return total*5<max_total*4;
}

Integer curvature_upper(Interval r,Interval s,Interval c,Integer fallback) {
    const auto ratio=r/s;
    auto cosh_ratio=[](Interval q){return (q+one/q)/2;};
    const auto left=cosh_ratio(point(ratio.lo)), right=cosh_ratio(point(ratio.hi));
    Interval q{ratio.lo<=unit && unit<=ratio.hi?unit:std::min(left.lo,right.lo),
               std::max(left.hi,right.hi)};
    const auto slope=rat(2)-square(q)-c*q;
    if (q.lo<=c.hi || (slope.lo<0 && slope.hi>0)) return fallback;
    const auto endpoint=point(slope.lo>=0?c.hi:c.lo);
    auto curvature=[&](Integer z) {
        auto t=point(z);
        return ((one-t*endpoint)*2/square(t-endpoint)).hi;
    };
    return std::min(fallback,std::max(curvature(q.lo),curvature(q.hi)));
}

bool radial_hessian_floor(const Radii &r,const Edges &e,const Derivatives &data,
                          Interval inner,Interval outer) {
    std::array<std::array<Integer,5>,5> lower{};
    for (int i=0;i<5;++i) lower[i][i]=10*r[i].lo-(i==0?inner:outer).hi;
    for (int k=0;k<10;++k) {
        const auto [i,j,mid,range]=e[k];
        const auto b=curvature_upper(r[i],r[j],range.cos,data.b[k].hi);
        lower[i][i]-=b;
        lower[j][j]-=b;
        lower[i][j]=b;
    }
    return positive_definite(lower);
}

struct RadialEllipsoid {
    Radii stiffness{}, center{};
    Interval budget{};
    bool available=false;
};

bool local_energy_contraction(Radii &r,const EnergyData &data,const Edges &e,
                              const Radii &ref,const Derivatives &fixed,RadialEllipsoid &ellipsoid) {
    Radii hull;
    for (int i=0;i<5;++i) hull[i]={std::min(r[i].lo,ref[i].lo),std::max(r[i].hi,ref[i].hi)};
    Derivatives ranges;
    if (!derivatives(hull,e,false,ranges)) return true;
    for (const auto &[inner,outer]:std::array<std::pair<Interval,Interval>,3>{{
        {rat(3,4),rat(3,2)},{rat(1,2),one},{rat(1,3),rat(2,3)}
    }}) {
        if (!radial_hessian_floor(hull,e,ranges,inner,outer)) continue;
        ++local_hessian_refinements;
        auto budget=(data.product/product_witness-one+data.angular)*2;
        for (int i=0;i<5;++i)
            budget=budget+square(fixed.g[i])/(i==0?inner:outer);
        if (budget.hi<0) return false;
        for (int i=0;i<5;++i) {
            const auto stiffness=i==0?inner:outer;
            ellipsoid.stiffness[i]=stiffness;
            ellipsoid.center[i]=fixed.g[i]/stiffness;
            const auto radius=sqrt_upper((budget/stiffness).hi);
            if (!meet(r[i],ref[i]+(ellipsoid.center[i]+Interval{-radius,radius})*hull[i]))
                return false;
        }
        ellipsoid.budget=budget;
        ellipsoid.available=true;
        return true;
    }
    return true;
}

Interval log_relative_point(Interval r,Interval ref) {
    const auto t=(r-ref)/(r+ref);
    if (mag(t)>unit/4) throw std::logic_error("radius ratio for logarithm");
    const auto t2=square(t);
    Interval series=t, term=t;
    for (int k=1;k<=4;++k) {
        term=term*t2;
        series=series+term/(2*k+1);
    }
    const auto error=mag((term*t2*2)/(one-t2)/11);
    return series*2+Interval{-error,error};
}

Interval log_relative(Interval r,Interval ref) {
    return {log_relative_point(point(r.lo),ref).lo,log_relative_point(point(r.hi),ref).hi};
}

bool tighten_ellipsoid(Radii &r,const Radii &ref,const RadialEllipsoid &ellipsoid,Interval target) {
    for (int i=0;i<5;++i) {
        auto remaining=target;
        for (int j=0;j<5;++j) if (j!=i) remaining=remaining-r[j];
        if (!meet(r[i],remaining)) return false;
    }
    Radii contribution;
    Interval cost;
    for (int i=0;i<5;++i) {
        contribution[i]=ellipsoid.stiffness[i]*square(log_relative(r[i],ref[i])-ellipsoid.center[i]);
        cost=cost+contribution[i];
    }
    if (cost.lo>ellipsoid.budget.hi) return false;
    for (int i=0;i<5;++i) {
        const auto budget=point(ellipsoid.budget.hi-cost.lo+contribution[i].lo);
        const auto radius=sqrt_upper((budget/ellipsoid.stiffness[i]).hi);
        const Interval hull{std::min(r[i].lo,ref[i].lo),std::max(r[i].hi,ref[i].hi)};
        if (!meet(r[i],ref[i]+(ellipsoid.center[i]+Interval{-radius,radius})*hull))
            return false;
    }
    return true;
}

bool radial_subcover(Radii r,const MidpointAngles &angle,const Turn &directions,const Radii &ref,
                     const RadialEllipsoid &ellipsoid,Interval target,Interval threshold,int depth) {
    ++radial_subcover_nodes;
    if (!tighten_ellipsoid(r,ref,ellipsoid,target)) return true;
    if (requires_zero_midpoints(angle,r,threshold) || paired_orbit(directions,r)!=0) return true;
    if (!depth) return false;
    int axis=0;
    for (int i=1;i<5;++i)
        if (r[i].hi-r[i].lo>r[axis].hi-r[axis].lo) axis=i;
    const auto midpoint=center(r[axis]).lo, upper=r[axis].hi;
    if (midpoint==r[axis].lo || midpoint==upper) return false;
    r[axis].hi=midpoint;
    if (!radial_subcover(r,angle,directions,ref,ellipsoid,target,threshold,depth-1)) return false;
    r[axis]={midpoint,upper};
    return radial_subcover(r,angle,directions,ref,ellipsoid,target,threshold,depth-1);
}

struct EdgeBounds {
    Integer largest=0;
    Radii incident{};

    void include(int i,int j,Interval bound) {
        largest=std::max(largest,bound.hi);
        incident[i]=incident[i]+bound;
        incident[j]=incident[j]+bound;
    }

    Interval quartic_bound() const {
        const auto peak=std::max_element(incident.begin(),incident.end(),
                                        [](Interval a,Interval b){return a.hi<b.hi;})->hi;
        return point(std::min(ceil_ratio(Wide(25)*largest,4),4*peak))*square(midpoint_sq);
    }
};

struct ProjectionErrors { Interval energy,stationarity; };

ProjectionErrors projection_errors(const Radii &r,const Edges &edges) {
    EdgeBounds quartic, sixth;
    const auto eta_sq=midpoint_sq*2;
    for (const auto &edge:edges) {
        const auto i=edge.i, j=edge.j;
        const auto p=r[i]*r[j]*2;
        const auto twice_dot_sq=p*(one+edge.range.cos);
        const auto twice_dot=Interval{std::max<Integer>(0,sqrt_upper(twice_dot_sq.lo)-1),
                                       sqrt_upper(twice_dot_sq.hi)};
        const auto long_dist=r[i]+r[j]+twice_dot;
        const auto product=square(r[i]-r[j])+p*(one-edge.range.cos);
        const auto short_dist=product/long_dist;
        Interval energy,stationarity;
        for (const auto floor:{point(std::max(short_dist.lo,square(rat(4,7)).lo)),
                               point(std::max(long_dist.lo,square(rat(99,100)).lo))}) {
            const auto rho=eta_sq/floor;
            energy=energy+(one+rho)/(square(floor)*(one-rho));
            stationarity=stationarity+rat(2)/(floor*square(floor-eta_sq));
        }
        quartic.include(i,j,energy);
        sixth.include(i,j,stationarity);
    }
    return {point(std::min(projection_slack.hi,quartic.quartic_bound().hi)),
            point(std::min(midpoint_test.hi,sixth.quartic_bound().hi*2))};
}

bool projection_energy_contraction(Radii &r,const Angles &a,const Edges &edges,const Radii &ref,
                                   Interval slack,RadialEllipsoid &ellipsoid) {
    // The raw projection need not belong to K. Certify its curvature only
    // along the segment from the reference to its current radial box.
    Derivatives mid,fixed;
    if (!derivatives(ref,edges,true,mid) || !derivatives(ref,edges,false,fixed)) return true;
    const auto data=energy_data(mid,fixed,edges,ref,a,slack);
    return local_energy_contraction(r,data,edges,ref,fixed,ellipsoid);
}

Vec multiply(const Matrix &m,const Vec &x,int rows) {
    Vec out{};
    for (int i=0;i<rows;++i) for (int j=0;j<9;++j) out[i]=out[i]+m[i][j]*x[j];
    return out;
}

struct RadialModel { Matrix preconditioner{}; Vec residual{}; };

Integer preconditioned_force_slack(const Vec &row) {
    Interval norm_sq, centered_norm_sq, sum;
    for (int j=0;j<5;++j) sum=sum+row[j];
    const auto average=center(sum/5);
    for (int j=0;j<9;++j) {
        norm_sq=norm_sq+square(row[j]);
        centered_norm_sq=centered_norm_sq+square(row[j]-(j<5?average:Interval{}));
    }
    const auto direct=point(sqrt_upper(norm_sq.hi))*force_slack;
    const auto centered=point(sqrt_upper(centered_norm_sq.hi))*force_slack
                       +point(mag(average))*individual_midpoint_sq*10;
    return std::min(direct.hi,centered.hi);
}

bool radial_model(const Derivatives &mid,const Derivatives &fixed,const Edges &e,
                  const Radii &ref,const Angles &a,const FloatMatrix &jac,RadialModel &model) {
    auto &p=model.preconditioner;
    if (!propose_matrix(jac,5,p)) return false;
    auto &error=model.residual;
    error=multiply(p,mid.g,5);
    auto mixed=transformed(p,5,fixed,ref,e,false,true);
    for (int i=0;i<5;++i)
        for (int j=0;j<4;++j) error[i]=error[i]+mixed[i][j+5]*(a[j]-center(a[j]));
    return true;
}

bool radial_inclusion(const Edges &e,const Radii &ref,const RadialModel &model,
                      bool approximate,Radii &out,const ForceAllowance *known_noise=nullptr) {
    auto error=model.residual;
    std::array<Integer,5> widths{};
    const int padding=approximate?64:8;
    for (int i=0;i<5;++i) {
        if (approximate) {
            const auto noise=known_noise?(*known_noise)[i]:preconditioned_force_slack(model.preconditioner[i]);
            error[i]=error[i]+Interval{-noise,noise};
        }
        widths[i]=mag(error[i])+std::max<Integer>(128,mag(error[i])/padding);
    }
    for (int k=0;k<(approximate?12:6);++k) {
        for (int i=0;i<5;++i)
            out[i]=ref[i]+Interval{-widths[i],widths[i]};
        if (!inside_log_domain(out)) return false;
        Derivatives ranges;
        if (!derivatives(out,e,false,ranges)) return false;
        auto rad=transformed(model.preconditioner,5,ranges,out,e,true,false);
        bool enclosed=true;
        Vec image{};
        for (int i=0;i<5;++i) {
            image[i]=-error[i];
            for (int j=0;j<5;++j)
                image[i]=image[i]+((i==j?one:Interval{})-rad[i][j])*Interval{-widths[j],widths[j]};
            enclosed &= mag(image[i])<widths[i];
        }
        if (enclosed) {
            // P is nonsingular when I-P Dg is a contraction in the weighted sup norm.
            // The inclusion bounds used above also bound its weighted norm strictly below 1.
            for (int i=0;i<5;++i) out[i]=ref[i]+image[i];
            return true;
        }
        for (int i=0;i<5;++i)
            widths[i]=std::max(widths[i],mag(image[i])+std::max<Integer>(128,mag(image[i])/padding));
    }
    return false;
}

bool tighten_approximate_radii(Radii &r,const Radii &ref,const Edges &e,const RadialModel &model) {
    const Interval repaired_mass{mass.lo-individual_midpoint_sq.hi,mass.hi};
    for (int pass=0;pass<3;++pass) {
        for (int i=0;i<5;++i) {
            Interval total=repaired_mass;
            for (int j=0;j<5;++j) if (i!=j) total=total-r[j];
            if (!meet(r[i],total)) return false;
        }
        Radii hull;
        for (int i=0;i<5;++i) hull[i]={std::min(r[i].lo,ref[i].lo),std::max(r[i].hi,ref[i].hi)};
        Derivatives ranges;
        if (!derivatives(hull,e,false,ranges)) return true;
        const auto &p=model.preconditioner;
        const auto derivative=transformed(p,5,ranges,hull,e,true,false);
        auto next=r;
        for (int i=0;i<5;++i) {
            const auto noise=preconditioned_force_slack(p[i]);
            Interval image=ref[i]-model.residual[i]+Interval{-noise,noise};
            for (int j=0;j<5;++j)
                image=image+((i==j?one:Interval{})-derivative[i][j])*(r[j]-ref[j]);
            if (!meet(next[i],image)) return false;
        }
        r=next;
    }
    return true;
}

// Bound each row of P times the repaired radial force using its own coefficients.
ForceAllowance radial_force_noise(const Angles &angle,const Radii &r,const Matrix &p) {
    ForceAllowance result;
    for (int i=0;i<5;++i) result[i]=preconditioned_force_slack(p[i]);
    const auto extended_angles=unscaled_angle_box(angle);
    const auto extended=unscaled_radius_box(r);
    const auto edges=make_edges(extended_angles);
    Derivatives ranges;
    if (!derivatives(extended,edges,false,ranges)) return result;
    const auto derivative=transformed(p,5,ranges,extended,edges);
    const Interval remainder=midpoint_sq*2/(square(rat(4,7))-midpoint_sq*2)*15;
    std::array<Interval,10> inverse_cube;
    for (int k=0;k<10;++k) {
        auto d=ranges.distance[k];
        const auto lower=point(sqrt_upper(d.lo)-1);
        if (lower.lo<=0) return result;
        inverse_cube[k]=one/(point(d.lo)*lower);
    }
    for (int row=0;row<5;++row) {
        std::array<Interval,5> incident{};
        Interval max_edge;
        for (int k=0;k<10;++k) {
            int i=edges[k].i,j=edges[k].j;
            auto ri=extended[i], rj=extended[j];
            auto q=scale(ComplexInterval{edges[k].range.cos,-edges[k].range.sin},rj);
            auto left=scale({ri+q.re*3,q.im*3},p[row][i]*ri);
            auto right=scale(q*ComplexInterval{ri*3+q.re,q.im},p[row][j]);
            Interval bound=point(magnitude(left-right))*inverse_cube[k]*4
                           +point(mag(p[row][i])+mag(p[row][j]))*remainder;
            max_edge.hi=std::max(max_edge.hi,bound.hi);
            incident[i]=incident[i]+bound;
            incident[j]=incident[j]+bound;
        }
        const auto max_incident=std::max_element(incident.begin(),incident.end(),
                                                 [](auto a,auto b){return a.hi<b.hi;});
        Interval repair;
        for (int i=1;i<5;++i) {
            Interval bound=point(mag(derivative[row][i]))+point(mag(derivative[row][i+4]))*rat(6,5)/tau;
            repair.hi=std::max(repair.hi,bound.hi);
        }
        Integer total=std::min(max_edge.hi*5,max_incident->hi*2)+repair.hi;
        result[row]=std::min(result[row],(point(total)*midpoint_sq).hi);
    }
    return result;
}

bool tighten_repaired_radii(Radii &r,const Angles &angle,const Radii &ref,const Edges &e,
                            const RadialModel &model) {
    const auto noise=radial_force_noise(angle,r,model.preconditioner);
    // The old box bounds the actual force error. For that fixed error, radial
    // uniqueness also applies to a newly certified inclusion box inside K.
    Radii smaller;
    if (radial_inclusion(e,ref,model,true,smaller,&noise)) {
        ++sharper_radial_inclusions;
        for (int i=0;i<5;++i) if (!meet(r[i],smaller[i])) return false;
    }
    for (int pass=0;pass<2;++pass) {
        Radii hull;
        for (int i=0;i<5;++i) hull[i]={std::min(r[i].lo,ref[i].lo),std::max(r[i].hi,ref[i].hi)};
        Derivatives ranges;
        if (!derivatives(hull,e,false,ranges)) return true;
        auto derivative=transformed(model.preconditioner,5,ranges,hull,e,true,false);
        auto next=r;
        for (int i=0;i<5;++i) {
            Interval image=ref[i]-model.residual[i]+Interval{-noise[i],noise[i]};
            for (int j=0;j<5;++j)
                image=image+((i==j?one:Interval{})-derivative[i][j])*(r[j]-ref[j]);
            if (!meet(next[i],image)) return false;
        }
        r=next;
    }
    return true;
}

enum class Result {none,force,krawczyk};
Result full_test(Angles &a,const Edges &e,const Radii &ref,const Radii &r,const Derivatives &mid,
                 const FloatMatrix &jac) {
    Derivatives ranges;
    if (!derivatives(r,e,false,ranges)) return Result::none;
    for (auto f:ranges.g) if (outside_force_slack(f)) return Result::force;
    Vec delta{}, x0{};
    for (int i=0;i<5;++i) { delta[i]=r[i]-ref[i]; x0[i]=ref[i]; }
    for (int i=0;i<4;++i) { x0[i+5]=center(a[i]); delta[i+5]=a[i]-x0[i+5]; }
    Matrix identity{};
    for (int i=0;i<9;++i) identity[i][i]=one;
    // The centered line can include ref even when the final Krawczyk image
    // (whose intersection contains every root) is offset from ref.
    Radii hull=r;
    for (int i=0;i<5;++i) hull[i]={std::min(r[i].lo,ref[i].lo),std::max(r[i].hi,ref[i].hi)};
    if (!derivatives(hull,e,false,ranges)) return Result::none;
    auto jwhole=transformed(identity,9,ranges,hull,e);
    auto mean=multiply(jwhole,delta,9);
    for (int i=0;i<9;++i) if (outside_force_slack(mid.g[i]+mean[i])) return Result::force;
    Matrix p{};
    if (!propose_matrix(jac,9,p)) return Result::none;
    jwhole=transformed(p,9,ranges,hull,e);
    auto rhs=multiply(p,mid.g,9);
    std::array<Interval,9> image;
    Angles new_box=a;
    for (int i=0;i<9;++i) {
        image[i]=-rhs[i];
        if (small_midpoints) {
            const auto noise=preconditioned_force_slack(p[i]);
            image[i]=image[i]+Interval{-noise,noise};
        }
        for (int j=0;j<9;++j)
            image[i]=image[i]+((i==j?one:Interval{})-jwhole[i][j])*delta[j];
        if (!meet(image[i],delta[i])) return Result::krawczyk;
        if (i>=5) new_box[i-5]=x0[i]+image[i];
    }
    a=new_box;
    return Result::none;
}

void visit(Angles a,int depth,const Radii &seed) {
    ++nodes;
    deepest=std::max(deepest,depth);
    if (!gap_constraints(a)) { ++geometry; return; }
    auto e=make_edges(a);
    if (!geometry_ok(e)) { ++geometry; return; }
    const int terminal=small_midpoints?0:orbit(a);
    if (terminal) { ++(terminal==1?basin:saddle); return; }
    Radii ref=reference(e,seed);
    Derivatives mid,fixed;
    if (derivatives(ref,e,true,mid) && derivatives(ref,e,false,fixed)) {
        const auto energy_bound=energy_data(mid,fixed,e,ref,a,small_midpoints?energy_slack:Interval{});
        if (energy_excluded(energy_bound)) { ++energy; return; }
        auto jac=float_jacobian(mid,ref,e);
        Radii r,energy_box;
        const bool energy_localized=small_midpoints && energy_radii(energy_bound,ref,fixed,energy_box);
        if (energy_localized) ++energy_radial_localizations;
        RadialModel model;
        if (radial_model(mid,fixed,e,ref,a,jac,model)) {
            const bool approximate=small_midpoints && radial_inclusion(e,ref,model,true,r);
            const bool included=approximate || radial_inclusion(e,ref,model,false,r);
            if (included || energy_localized) {
                if (included) ++radial_inclusions;
                if (small_midpoints) {
                    if (!included) r=energy_box;
                    if (included && !approximate) for (int i=0;i<5;++i) {
                        auto bound=i==0?rat(53,400):rat(1271,1000);
                        Integer inflation=(force_slack*4*bound).hi;
                        r[i]=r[i]+Interval{-inflation,inflation};
                        if (!meet(r[i],i==0?inner_band:outer_band)) throw std::logic_error("radial bounds");
                    }
                    if (included && energy_localized) for (int i=0;i<5;++i) {
                        if (!meet(r[i],energy_box[i])) { ++energy; return; }
                    }
                    if (!tighten_approximate_radii(r,ref,e,model) ||
                        !tighten_repaired_radii(r,a,ref,e,model)) { ++force; return; }
                    RadialEllipsoid ellipsoid;
                    if (!local_energy_contraction(r,energy_bound,e,ref,fixed,ellipsoid)) { ++energy; return; }
                    const Interval repaired_mass{mass.lo-individual_midpoint_sq.hi,mass.hi};
                    if (ellipsoid.available && !tighten_ellipsoid(r,ref,ellipsoid,repaired_mass)) {
                        ++energy; return;
                    }
                    const auto projected_angles=unscaled_angle_box(a);
                    const auto midpoint_directions=midpoint_angles(projected_angles);
                    auto projected=unscaled_radius_box(r);
                    if (requires_zero_midpoints(midpoint_directions,projected,midpoint_test)) {
                        ++zero_midpoints; return;
                    }
                    const int paired=paired_orbit(angle_directions(a),r,midpoint_sq);
                    if (paired) { ++(paired==1?basin:saddle); return; }
                    const bool small_angles=std::all_of(a.begin(),a.end(),[](Interval x){
                        return x.hi-x.lo<=unit/1024;
                    });
                    if (small_angles) {
                        RadialEllipsoid projection;
                        const auto projection_edges=make_edges(projected_angles);
                        const auto errors=projection_errors(projected,projection_edges);
                        const Interval projected_mass{mass.lo-midpoint_sq.hi,mass.hi};
                        if (!projection_energy_contraction(projected,projected_angles,projection_edges,ref,
                                                           errors.energy,projection) ||
                            (projection.available && !tighten_ellipsoid(projected,ref,projection,projected_mass))) {
                            ++energy; return;
                        }
                        if (requires_zero_midpoints(midpoint_directions,projected,errors.stationarity)) {
                            ++zero_midpoints; return;
                        }
                        const Turn directions=angle_directions(projected_angles);
                        const int projected_orbit=paired_orbit(directions,projected);
                        if (projected_orbit) { ++(projected_orbit==1?basin:saddle); return; }
                        if (projection.available &&
                            radial_subcover(projected,midpoint_directions,directions,ref,projection,
                                            projected_mass,errors.stationarity,9)) {
                            ++radial_subcovers; return;
                        }
                    }
                }
                auto result=full_test(a,e,ref,r,mid,jac);
                if (result==Result::force) { ++force; return; }
                if (result==Result::krawczyk) { ++krawczyk; return; }
            }
        }
    }
    if (depth>=max_depth) { ++unknown; return; }
    int axis=0;
    for (int j=1;j<4;++j)
        if (a[j].hi-a[j].lo>a[axis].hi-a[axis].lo) axis=j;
    Integer m=center(a[axis]).lo, hi=a[axis].hi;
    if (m==a[axis].lo || m==hi) { ++unknown; return; }
    ++split;
    a[axis].hi=m;
    visit(a,depth+1,ref);
    a[axis]={m,hi};
    visit(a,depth+1,ref);
}

int main(int argc,char **argv) {
    try {
        if (argc!=2 && argc!=3) return 2;
        max_depth=std::stoi(argv[1]);
        if (argc==3) {
            if (std::string(argv[2])!="midpoints") return 2;
            small_midpoints=true;
        }
        if (max_depth<0 || max_depth>30) return 2;
        if (direction(lower_gap).cos.lo<=gap_cos.hi ||
            direction(upper_gap).cos.hi>=gap_cos.lo)
            throw std::logic_error("gap bounds in turns");
        Integer denominator;
        int sixth_order, energy_coefficient, projection_coefficient;
        if (!(std::cin>>denominator>>sixth_order>>energy_coefficient>>projection_coefficient) ||
            denominator<=0 || sixth_order<=0 || energy_coefficient<=0 || projection_coefficient<=0) return 2;
        const auto midpoint_limit=rat(1,denominator);
        midpoint_sq=square(midpoint_limit);
        individual_midpoint_sq=midpoint_sq*4/5;
        force_slack=midpoint_sq*123;
        energy_slack=square(midpoint_sq)*energy_coefficient;
        projection_slack=square(midpoint_sq)*projection_coefficient;
        midpoint_test=square(midpoint_sq)*sixth_order;
        for (auto &row:anchors) for (auto &v:row) {
            if (!(std::cin>>v.cos.lo>>v.cos.hi>>v.sin.lo>>v.sin.hi)) return 2;
        }
        for (auto &row:anchor_radii) for (auto &r:row) {
            if (!(std::cin>>r.lo>>r.hi)) return 2;
        }
        Radii r;
        r[0]=center(rat(12,100));
        for (int i=1;i<4;++i) r[i]=center(rat(1095,1000));
        r[4]=point(mass.lo-r[0].lo-r[1].lo-r[2].lo-r[3].lo);
        if (!inside_log_domain(r)) throw std::logic_error("fallback reference");
        int roots=0;
        constexpr Integer bins=32, width=unit/bins;
        for (int i=0;i<bins;++i) for (int j=i+1;j<bins;++j)
            for (int k=j+1;k<bins;++k) for (int l=k+1;l<bins;++l) {
                ++roots;
                visit({Interval{i*width,(i+1)*width},Interval{j*width,(j+1)*width},
                       Interval{k*width,(k+1)*width},Interval{l*width,(l+1)*width}},0,r);
            }
        if (nodes!=roots+2*split ||
            nodes!=split+geometry+energy+force+krawczyk+basin+saddle+unknown+zero_midpoints+radial_subcovers)
            throw std::logic_error("subdivision counts");
        std::cout<<roots<<" "<<nodes<<" "<<split<<" "<<geometry<<" "<<energy<<" "<<force<<" "<<krawczyk
                 <<" "<<basin<<" "<<saddle<<" "<<unknown<<" "<<zero_midpoints<<" "<<radial_subcovers
                 <<" "<<radial_inclusions<<" "<<energy_radial_localizations<<" "<<local_hessian_refinements
                 <<" "<<sharper_radial_inclusions<<" "<<radial_subcover_nodes<<" "<<deepest<<"\n";
    } catch (const std::exception &e) { std::cerr<<e.what()<<"\n"; return 1; }
}
