// Four-angle cover. All conclusions use the integer interval operations.
#include "dyadic.hpp"

#include <array>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

using Angles = std::array<Interval,4>;
using Radii = std::array<Interval,5>;
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
int deepest=0, max_depth=0;
bool small_midpoints=false;
Interval force_slack, energy_slack, midpoint_sq;

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

int paired_orbit(const Angles &angle,const Radii &r) {
    Turn d;
    for (int j=0;j<4;++j) d[j]=direction(angle[j]);
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
            Interval bound=k==0?rat(27,1000):rat(26,1000);
            if (total.hi<square(bound).lo) return k+1;
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
    const auto shift=(midpoint_sq*6/5)/tau;
    for (auto &v:a) v=v+Interval{-shift.hi,shift.hi};
    return a;
}

Radii unscaled_radius_box(Radii r) {
    for (int i=1;i<5;++i) r[i]=r[i]-Interval{0,midpoint_sq.hi};
    return r;
}

using RealQuadratic=std::array<std::array<Integer,10>,10>;
bool positive_definite(RealQuadratic a) {
    for (int i=0;i<10;++i) {
        if (a[i][i]<=0) return false;
        for (int j=i+1;j<10;++j) {
            if (a[j][j]<=0 || Wide(a[i][j])*a[i][j]>=Wide(a[i][i])*a[j][j]) return false;
        }
        for (int j=i+1;j<10;++j) {
            for (int k=j;k<10;++k)
                a[j][k]=checked(Wide(a[j][k])-ceil_ratio(Wide(a[i][j])*a[i][k],a[i][i]));
            a[j][j]-=(9-i);
        }
    }
    return true;
}

bool requires_zero_midpoints(const Angles &angle,const Radii &repaired) {
    // Return to the unscaled half differences using the componentwise repair bounds.
    std::array<Interval,5> a{};
    const auto raw=unscaled_angle_box(angle);
    const auto r=unscaled_radius_box(repaired);
    for (int i=1;i<5;++i) a[i]=raw[i-1];
    RealQuadratic lower{};
    for (int i=0;i<10;++i) lower[i][i]=10*unit-(midpoint_sq*147).hi;
    for (int i=0;i<5;++i) for (int j=i+1;j<5;++j) {
        auto u=direction(a[i]), relative=direction(a[j]-a[i]);
        auto d=square(r[i]-r[j])+r[i]*r[j]*(one-relative.cos)*2;
        if (d.lo<=0) return false;
        auto p=r[j]*relative.cos, q=r[j]*relative.sin;
        ComplexInterval numerator{r[i]+p,q}, conjugate_difference{r[i]-p,q};
        auto weight=numerator*conjugate_difference*conjugate_difference*ComplexInterval{u.cos,-u.sin};
        auto den=square(d);
        weight={weight.re*2/den,weight.im*2/den};
        auto re=center(weight.re), im=center(weight.im);
        Integer error=sqrt_upper((square(weight.re-re)+square(weight.im-im)).hi);
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

bool energy_excluded(const Derivatives &mid, const Derivatives &fixed, const Edges &edges,
                     const Radii &ref, const Angles &a) {
    Vec disp{};
    for (int i=0;i<4;++i) disp[i+5]=(a[i]-center(a[i]))*tau;
    Interval increase;
    if (small_midpoints) increase=energy_slack;
    for (int i=0;i<5;++i) increase=increase+square(fixed.g[i])*2;
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
    Interval bound=point(std::max<Integer>(0,increase.hi));
    return bound.hi<=unit && (product*exp_upper(bound)).hi<rat(469765409,1LL<<26).lo;
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
                       +point(mag(average))*midpoint_sq*10;
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
                      bool approximate,Radii &out) {
    auto error=model.residual;
    std::array<Integer,5> widths{};
    const int padding=approximate?64:8;
    for (int i=0;i<5;++i) {
        if (approximate) {
            const auto noise=preconditioned_force_slack(model.preconditioner[i]);
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
    const Interval repaired_mass{mass.lo-midpoint_sq.hi,mass.hi};
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
std::array<Integer,5> radial_force_noise(const Angles &angle,const Radii &r,const Matrix &p) {
    std::array<Integer,5> result;
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
        if (energy_excluded(mid,fixed,e,ref,a)) { ++energy; return; }
        auto jac=float_jacobian(mid,ref,e);
        Radii r;
        RadialModel model;
        if (radial_model(mid,fixed,e,ref,a,jac,model)) {
            const bool approximate=small_midpoints && radial_inclusion(e,ref,model,true,r);
            if (approximate || radial_inclusion(e,ref,model,false,r)) {
                ++radial_inclusions;
                if (small_midpoints) {
                    if (!approximate) for (int i=0;i<5;++i) {
                        auto bound=i==0?rat(53,400):rat(1271,1000);
                        Integer inflation=(force_slack*4*bound).hi;
                        r[i]=r[i]+Interval{-inflation,inflation};
                        if (!meet(r[i],i==0?inner_band:outer_band)) throw std::logic_error("radial bounds");
                    }
                    if (!tighten_approximate_radii(r,ref,e,model) ||
                        !tighten_repaired_radii(r,a,ref,e,model)) { ++force; return; }
                    if (requires_zero_midpoints(a,r)) { ++zero_midpoints; return; }
                    const int paired=paired_orbit(a,r);
                    if (paired) { ++(paired==1?basin:saddle); return; }
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
        if (!(std::cin>>denominator) || denominator<=0) return 2;
        const auto midpoint_limit=rat(1,denominator);
        midpoint_sq=square(midpoint_limit);
        force_slack=midpoint_sq*132;
        energy_slack=square(midpoint_sq)*242;
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
        if (nodes!=roots+2*split || nodes!=split+geometry+energy+force+krawczyk+basin+saddle+unknown+zero_midpoints)
            throw std::logic_error("subdivision counts");
        std::cout<<roots<<" "<<nodes<<" "<<split<<" "<<geometry<<" "<<energy<<" "<<force<<" "<<krawczyk
                 <<" "<<basin<<" "<<saddle<<" "<<unknown<<" "<<zero_midpoints
                 <<" "<<radial_inclusions<<" "<<deepest<<"\n";
    } catch (const std::exception &e) { std::cerr<<e.what()<<"\n"; return 1; }
}
