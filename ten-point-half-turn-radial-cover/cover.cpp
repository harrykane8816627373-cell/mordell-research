// Exact dyadic interval pruning in the angle-based half-turn chart.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>

[[noreturn]] void fail(int line) {
    std::cerr << "Arithmetic range or subdivision failure at source line " << line << std::endl;
    std::abort();
}

using I64 = int64_t;
using Z = __int128_t;
constexpr I64 S = 1LL << 43;
constexpr I64 LIMIT = (1LL << 60);

I64 floorq(Z a, Z b) {
    if (b <= 0) fail(__LINE__);
    Z q=a/b;
    if (a%b < 0) --q;
    if (q < -LIMIT || q > LIMIT) fail(__LINE__);
    return static_cast<I64>(q);
}
I64 ceilq(Z a, Z b) {return -floorq(-a,b);}

struct Iv {I64 l,h;};
Iv Q(I64 a, I64 b=1) {return {floorq(Z(a)*S,b),ceilq(Z(a)*S,b)};}
I64 safe(Z v) {
    if (v < -LIMIT || v > LIMIT) fail(__LINE__);
    return static_cast<I64>(v);
}
Iv operator+(Iv a,Iv b) {return {safe(Z(a.l)+b.l),safe(Z(a.h)+b.h)};}
Iv operator-(Iv a) {return {-a.h,-a.l};}
Iv operator-(Iv a,Iv b) {return a+-b;}
Iv operator*(Iv a,Iv b) {
    auto [l,h]=std::minmax({Z(a.l)*b.l,Z(a.l)*b.h,Z(a.h)*b.l,Z(a.h)*b.h});
    return {floorq(l,S),ceilq(h,S)};
}
Iv times(Iv a,int m) {return m>=0?Iv{safe(Z(a.l)*m),safe(Z(a.h)*m)}:Iv{safe(Z(a.h)*m),safe(Z(a.l)*m)};}
Iv operator/(Iv a,Iv b) {
    if(b.l<=0) fail(__LINE__);
    if(a.l>=0) return {floorq(Z(a.l)*S,b.h),ceilq(Z(a.h)*S,b.l)};
    if(a.h<=0) return {floorq(Z(a.l)*S,b.l),ceilq(Z(a.h)*S,b.h)};
    return {floorq(Z(a.l)*S,b.l),ceilq(Z(a.h)*S,b.l)};
}
Iv sq(Iv a) {
    I64 low=(a.l<=0 && a.h>=0)?0:std::min(std::abs(a.l),std::abs(a.h));
    I64 high=std::max(std::abs(a.l),std::abs(a.h));
    return {floorq(Z(low)*low,S),ceilq(Z(high)*high,S)};
}
I64 isqrt(Z a) {
    if(a<0) fail(__LINE__);
    // Double is a speed hint only. The two integer comparisons fix the exact result.
    I64 x=static_cast<I64>(std::sqrt(static_cast<double>(a)));
    while(Z(x)*x>a) --x;
    while(Z(x+1)*(x+1)<=a) ++x;
    return x;
}
Iv sqrt_iv(Iv a) {return {isqrt(Z(std::max<I64>(a.l,0))*S),isqrt(Z(a.h)*S)+1};}
bool cut(Iv &x,Iv y) {x.l=std::max(x.l,y.l); x.h=std::min(x.h,y.h); return x.l<=x.h;}
bool excludes0(Iv a) {return a.l>0 || a.h<0;}

constexpr Iv ZERO{0,0}, ONE{S,S}, UNIT{-S,S};
struct C {Iv x,y;};
C operator*(C a,C b) {return {a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x};}
bool clipunit(C &a) {return cut(a.x,UNIT) && cut(a.y,UNIT);}

using Box=std::array<Iv,6>;
enum Outcome {Split, Geometry, Force, Value, Undecided, Newton};
const Iv OUTER{Q(89,100).l,Q(1271,1000).h};
const Iv RING_COS{Q(-1519,5000).l,Q(271,1250).h};
const Iv TWOSTEP_COS{-S,Q(-81,100).h};
const Iv RING_SIN{Q(19,20).l,S};
const Iv DISPERSE=Q(47,2000);
const Iv WITNESS=Q(469765409,1LL<<50);

Iv ward_weight(Iv r) {return Q(289,40)-times(sq(r-Q(17,20)),10);}

constexpr int DIM=7;
using Gradient=std::array<Iv,DIM>;
Gradient parameter_delta;
struct Jet {
    Iv value;
    Iv center;
    Gradient derivative{};
    Jet():value(ZERO),center(ZERO) {}
    Jet(Iv v):value(v),center(v) {}
};
void centered_enclosure(Jet &r) {
    Iv value=r.center;
    for(int k=0;k<DIM;++k) value=value+r.derivative[k]*parameter_delta[k];
    if(!cut(r.value,value)) fail(__LINE__);
}
Jet operator+(const Jet &a,const Jet &b) {
    Jet r(a.value+b.value);
    r.center=a.center+b.center;
    for(int k=0;k<DIM;++k) r.derivative[k]=a.derivative[k]+b.derivative[k];
    centered_enclosure(r);
    return r;
}
Jet operator-(const Jet &a) {
    Jet r(-a.value);
    r.center=-a.center;
    for(int k=0;k<DIM;++k) r.derivative[k]=-a.derivative[k];
    return r;
}
Jet operator-(const Jet &a,const Jet &b) {return a+-b;}
Jet operator*(const Jet &a,const Jet &b) {
    Jet r(a.value*b.value);
    r.center=a.center*b.center;
    for(int k=0;k<DIM;++k) r.derivative[k]=a.derivative[k]*b.value+a.value*b.derivative[k];
    centered_enclosure(r);
    return r;
}
Jet operator/(const Jet &a,const Jet &b) {
    Jet r(a.value/b.value);
    r.center=a.center/b.center;
    for(int k=0;k<DIM;++k) r.derivative[k]=(a.derivative[k]-r.value*b.derivative[k])/b.value;
    centered_enclosure(r);
    return r;
}
Jet times(const Jet &a,int m) {
    Jet r(times(a.value,m));
    r.center=times(a.center,m);
    for(int k=0;k<DIM;++k) r.derivative[k]=times(a.derivative[k],m);
    centered_enclosure(r);
    return r;
}
Jet sq(const Jet &a) {
    Jet r(sq(a.value));
    r.center=sq(a.center);
    for(int k=0;k<DIM;++k) r.derivative[k]=times(a.value*a.derivative[k],2);
    centered_enclosure(r);
    return r;
}
Jet sqrt_jet(const Jet &a) {
    Jet r(sqrt_iv(a.value));
    r.center=sqrt_iv(a.center);
    for(int k=0;k<DIM;++k) r.derivative[k]=a.derivative[k]/times(r.value,2);
    centered_enclosure(r);
    return r;
}
struct JC {Jet x,y;};
JC operator*(const JC &a,const JC &b) {return {a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x};}
Jet weight_jet(const Jet &a) {return Jet(Q(289,40))-times(sq(a-Jet(Q(17,20))),10);}
Jet inverse_jet(const Jet &a) {return (Jet(Q(17))+sqrt_jet(Jet(Q(289))-times(a,40)))/Jet(Q(20));}

bool correlated_force(const std::array<Iv,DIM> &box,int side,std::array<Jet,11> &output) {
    std::array<Jet,DIM> x;
    for(int k=0;k<DIM;++k) {
        x[k]=Jet(box[k]);
        I64 midpoint=(box[k].l+box[k].h)/2;
        x[k].center={midpoint,midpoint};
        parameter_delta[k]=box[k]-x[k].center;
        x[k].derivative[k]=ONE;
    }
    const Jet one(ONE);
    auto direction=[&](const Jet &c) {return JC{c,sqrt_jet(one-sq(c))};};
    JC g1=direction(x[3]),g2=direction(x[4]),g3=direction(x[5]);
    JC h12=g1*g2,h23=g2*g3,last=h12*g3;
    Jet t2=sq(x[2]),div=one+t2;
    std::array<JC,5> u;
    u[0]={one,Jet(ZERO)};
    u[2]={times((one-t2)/div,side),times(x[2],2)/div};
    u[1]=u[2]*JC{g1.x,-g1.y};
    u[3]=u[2]*g2;
    u[4]=u[2]*h23;
    Jet b0=(Jet(Q(17))-times(x[0],10))*x[0],b1=weight_jet(x[1]);
    Jet a3=(-b0*u[4].y-b1*last.y)/g3.y,a4=(b0*u[3].y+b1*h12.y)/g3.y;
    Jet s3=-h23.y/g3.y,s4=g2.y/g3.y,b2=weight_jet(x[6]);
    Jet b3=a3+s3*b2,b4=a4+s4*b2;
    if((Q(289)-times(b3.value,40)).l<=Q(1,4).h ||
       (Q(289)-times(b4.value,40)).l<=Q(1,4).h) return false;
    std::array<Jet,5> r{x[0],x[1],x[6],inverse_jet(b3),inverse_jet(b4)},r2,real,imag;
    std::array<std::array<JC,5>,5> angles;
    for(int i=1;i<5;++i) angles[0][i]=u[i];
    angles[1][2]=g1; angles[2][3]=g2; angles[3][4]=g3;
    angles[1][3]=h12;angles[2][4]=h23;angles[1][4]=last;
    Jet mass=Jet(Q(-9,2));
    for(int i=0;i<5;++i) {
        r2[i]=sq(r[i]);real[i]=Jet(Q(9,2))-times(r[i],5);imag[i]=Jet(ZERO);
        mass=mass+r[i];
    }
    for(int i=0;i<5;++i) for(int j=i+1;j<5;++j) {
        Jet rr=r[i]*r[j], dd=sq(r[i]-r[j])+times(rr*(one-angles[i][j].x),2);
        if(dd.value.l<=Q(1,5).h) return false;
        Jet dr=(r2[i]-r2[j])/dd,di=times(rr*angles[i][j].y,2)/dd;
        real[i]=real[i]+dr;real[j]=real[j]-dr;
        imag[i]=imag[i]+di;imag[j]=imag[j]-di;
    }
    for(int i=0;i<5;++i) {
        output[2*i]=real[i];output[2*i+1]=imag[i];
    }
    output[10]=mass;
    return true;
}

using Preconditioner=std::array<std::array<Iv,DIM>,DIM>;
bool choose_inverse(const std::array<Jet,11> &center,std::array<int,DIM> &rows,Preconditioner &inverse) {
    std::array<std::array<double,11>,DIM> transpose{};
    std::array<int,11> perm;
    for(int j=0;j<11;++j) {
        perm[j]=j;
        for(int i=0;i<DIM;++i) transpose[i][j]=double(center[j].derivative[i].l+center[j].derivative[i].h)/(2*S);
    }
    // Complete pivoting selects a nonsingular seven-row proposal.
    for(int k=0;k<DIM;++k) {
        int a=k,b=k;
        for(int i=k;i<DIM;++i) for(int j=k;j<11;++j)
            if(std::abs(transpose[i][j])>std::abs(transpose[a][b])) {a=i;b=j;}
        if(std::abs(transpose[a][b])<1e-8) return false;
        std::swap(transpose[a],transpose[k]);
        for(int i=0;i<DIM;++i) std::swap(transpose[i][b],transpose[i][k]);
        std::swap(perm[b],perm[k]);
        for(int i=k+1;i<DIM;++i) {
            double v=transpose[i][k]/transpose[k][k];
            for(int j=k+1;j<11;++j) transpose[i][j]-=v*transpose[k][j];
        }
        rows[k]=perm[k];
    }
    std::array<std::array<double,2*DIM>,DIM> m{};
    for(int i=0;i<DIM;++i) {
        for(int j=0;j<DIM;++j) m[i][j]=double(center[rows[i]].derivative[j].l+center[rows[i]].derivative[j].h)/(2*S);
        m[i][DIM+i]=1;
    }
    for(int k=0;k<DIM;++k) {
        int pivot=k;
        for(int i=k+1;i<DIM;++i) if(std::abs(m[i][k])>std::abs(m[pivot][k])) pivot=i;
        if(std::abs(m[pivot][k])<1e-8) return false;
        std::swap(m[k],m[pivot]);
        double d=m[k][k];
        for(int j=k;j<2*DIM;++j) m[k][j]/=d;
        for(int i=0;i<DIM;++i) if(i!=k) {
            double v=m[i][k];
            for(int j=k;j<2*DIM;++j) m[i][j]-=v*m[k][j];
        }
    }
    for(int i=0;i<DIM;++i) for(int j=0;j<DIM;++j) {
        if(!std::isfinite(m[i][DIM+j]) || std::abs(m[i][DIM+j])>32) return false;
        I64 v=static_cast<I64>(std::round(m[i][DIM+j]*S));
        inverse[i][j]={v,v};
    }
    return true;
}

bool correlated_refine(Box &box,Iv radius,int side) {
    std::array<Iv,DIM> input,center,delta;
    for(int k=0;k<DIM;++k) {
        input[k]=k==6?radius:box[k];
        I64 mid=(input[k].l+input[k].h)/2;
        center[k]={mid,mid};delta[k]=input[k]-center[k];
    }
    std::array<Jet,11> whole,point;
    if(!correlated_force(input,side,whole)||!correlated_force(center,side,point)) return true;
    for(int i=0;i<11;++i) {
        Iv mean=point[i].value;
        for(int j=0;j<DIM;++j) {
            if(whole[i].derivative[j].l<-128*S || whole[i].derivative[j].h>128*S) return true;
            mean=mean+whole[i].derivative[j]*delta[j];
        }
        if(excludes0(whole[i].value)||excludes0(mean)) return false;
    }
    std::array<int,DIM> rows;
    Preconditioner inverse;
    if(!choose_inverse(point,rows,inverse)) return true;
    Preconditioner jacobian;
    Gradient rhs;
    for(int i=0;i<DIM;++i) {
        rhs[i]=ZERO;
        for(int k=0;k<DIM;++k) rhs[i]=rhs[i]-inverse[i][k]*point[rows[k]].value;
        for(int j=0;j<DIM;++j) {
            jacobian[i][j]=ZERO;
            for(int k=0;k<DIM;++k) jacobian[i][j]=jacobian[i][j]+inverse[i][k]*whole[rows[k]].derivative[j];
        }
    }
    for(int i=0;i<DIM;++i) {
        Iv proposed=center[i]+rhs[i];
        for(int j=0;j<DIM;++j) {
            proposed=proposed+(Q(i==j?1:0)-jacobian[i][j])*delta[j];
        }
        if(!cut(input[i],proposed)) return false;
    }
    for(int i=0;i<DIM;++i) delta[i]=input[i]-center[i];
    for(int iteration=0;iteration<2;++iteration) for(int i=0;i<DIM;++i) {
        Iv diagonal=jacobian[i][i];
        if(diagonal.l<=Q(1,8).h && diagonal.h>=-Q(1,8).h) continue;
        Iv value=rhs[i];
        for(int j=0;j<DIM;++j) if(j!=i) value=value-jacobian[i][j]*delta[j];
        if(value.l<-1024*S || value.h>1024*S) continue;
        Iv proposal=diagonal.l>0?value/diagonal:(-value)/(-diagonal);
        if(!cut(delta[i],proposal)) return false;
        input[i]=center[i]+delta[i];
    }
    for(int k=0;k<6;++k) box[k]=input[k];
    return true;
}

using Angles=std::array<std::array<C,5>,5>;

Outcome physical_test(const std::array<Iv,5> &r,const Angles &angles) {
    I64 lo=r[1].l,hi=r[1].h;
    for(int i=2;i<5;++i) {lo=std::min(lo,r[i].l);hi=std::max(hi,r[i].h);}
    if(hi-lo<=DISPERSE.l) return Geometry;
    std::array<Iv,5> real{}, imag{}, r2{};
    for(int i=0;i<5;++i) {r2[i]=sq(r[i]); real[i]=Q(9,2)-times(r[i],5); imag[i]=ZERO;}
    Iv product=ONE;
    for(int i=0;i<5;++i) product=product*r[i];
    for(int i=0;i<5;++i) for(int j=i+1;j<5;++j) {
        Iv co=angles[i][j].x, si=angles[i][j].y;
        Iv rr=r[i]*r[j];
        Iv radial_cut=sq(r[i]+r[j]-Q(2,5))-times(rr*(ONE+co),2);
        if(radial_cut.h<=0) return Geometry;
        Iv dd=sq(r[i]-r[j])+times(rr*(ONE-co),2);
        if(dd.l<=0) {
            if(i!=0) fail(__LINE__);
            if(!cut(dd,Iv{sq(Q(89,100)-Q(53,400)).l,LIMIT})) return Geometry;
        }
        Iv radd=(r2[i]-r2[j])/dd, iadd=times(rr*si,2)/dd;
        real[i]=real[i]+radd; real[j]=real[j]-radd;
        imag[i]=imag[i]+iadd; imag[j]=imag[j]-iadd;
        product=product*sq(dd/Q(4));
    }
    for(int i=0;i<5;++i) if(excludes0(real[i])||excludes0(imag[i])) return Force;
    if(product.h<WITNESS.l) return Value;
    return Split;
}

bool relation(Iv &variable,Iv target,Iv constant,Iv slope) {
    if(!cut(target,constant+slope*variable)) return false;
    if(slope.l>Q(1,4096).h) return cut(variable,(target-constant)/slope);
    if(slope.h<-Q(1,4096).h) return cut(variable,(constant-target)/(-slope));
    return true;
}

Iv inverse_weight(Iv w) {return (Q(17)+sqrt_iv(Q(289)-times(w,40)))/Q(20);}

bool mass_newton(Iv &r, Iv a3, Iv s3, Iv a4, Iv s4, Iv rest) {
    // The full interval derivative is checked positive before interval Newton.
    // This includes extensions beyond the physically admissible mass root.
    Iv w=ward_weight(r);
    Iv d3=Q(289)-times(a3+s3*w,40);
    Iv d4=Q(289)-times(a4+s4*w,40);
    if(d3.l<=Q(1,100).h || d4.l<=Q(1,100).h) return true;
    Iv derivative=ONE+s3*(times(r,20)-Q(17))/sqrt_iv(d3)
                     +s4*(times(r,20)-Q(17))/sqrt_iv(d4);
    if(derivative.l<=Q(1,8).h) return true;
    I64 midpoint=(r.l+r.h)/2;
    Iv center{midpoint,midpoint}, wc=ward_weight(center);
    Iv dc3=Q(289)-times(a3+s3*wc,40);
    Iv dc4=Q(289)-times(a4+s4*wc,40);
    if(dc3.l<=Q(1,100).h || dc4.l<=Q(1,100).h) return true;
    Iv f=center+(Q(34)+sqrt_iv(dc3)+sqrt_iv(dc4))/Q(20)-rest;
    return cut(r,center-f/derivative);
}

Outcome radial_profile(Box &x,int side,bool newton) {
    auto direction=[](Iv c) {return C{c,sqrt_iv(ONE-sq(c))};};
    C g1=direction(x[3]), g2=direction(x[4]), g3=direction(x[5]);
    C h12=g1*g2, h23=g2*g3, last=h12*g3;
    if(!cut(last.x,RING_COS)||!cut(last.y,-RING_SIN) || !cut(h12.x,TWOSTEP_COS) || !cut(h23.x,TWOSTEP_COS)) return Geometry;
    Iv t=x[2], t2=sq(t), div=ONE+t2;
    std::array<C,5> u{};
    u[0]={ONE,ZERO};
    u[2]={times((ONE-t2)/div,side),times(t,2)/div};
    u[1]=u[2]*C{g1.x,-g1.y};
    u[3]=u[2]*g2;
    u[4]=u[2]*h23;
    for(int i=1;i<5;++i) if(!clipunit(u[i])) return Geometry;
    Iv b0=(Q(17)-times(x[0],10))*x[0], b1=ward_weight(x[1]);
    Iv a3=(-b0*u[4].y-b1*last.y)/g3.y, a4=(b0*u[3].y+b1*h12.y)/g3.y;
    Iv s3=-h23.y/g3.y, s4=g2.y/g3.y;
    Iv b2=ward_weight(Iv{OUTER.l,x[1].h});
    Iv bband=ward_weight(OUTER);
    if(!relation(b2,bband,a3,s3)||!relation(b2,bband,a4,s4)) return Geometry;
    Iv min3=ONE-s3, min4=ONE-s4;
    if((min3*b2-a3).h<0 || (min4*b2-a4).h<0) return Geometry;
    if(min3.l>Q(1,4096).h && !cut(b2,Iv{(a3/min3).l,LIMIT})) return Geometry;
    if(min4.l>Q(1,4096).h && !cut(b2,Iv{(a4/min4).l,LIMIT})) return Geometry;
    if(min4.h<-Q(1,4096).h && !cut(b2,Iv{-LIMIT,(-a4/-min4).h})) return Geometry;
    std::array<Iv,5> r{x[0],x[1],OUTER,OUTER,OUTER};
    Iv rest=Q(9,2)-r[0]-r[1];
    for(int iter=0;iter<3;++iter) {
        Iv b3=a3+s3*b2,b4=a4+s4*b2;
        if(!cut(b3,bband)||!cut(b4,bband)) return Geometry;
        if(!cut(r[2],inverse_weight(b2))||!cut(r[3],inverse_weight(b3))||!cut(r[4],inverse_weight(b4))) return Geometry;
        if(!cut(r[2],Iv{OUTER.l,std::min({r[1].h,r[3].h,r[4].h})}) ||
           !cut(r[3],Iv{r[2].l,OUTER.h})||!cut(r[4],Iv{r[2].l,OUTER.h})) return Geometry;
        if(!cut(r[2],rest-r[3]-r[4]) || !cut(r[3],rest-r[2]-r[4]) || !cut(r[4],rest-r[2]-r[3])) return Geometry;
        if(!mass_newton(r[2],a3,s3,a4,s4,rest)) return Geometry;
        if(!cut(b2,ward_weight(r[2]))||!relation(b2,ward_weight(r[3]),a3,s3)||!relation(b2,ward_weight(r[4]),a4,s4)) return Geometry;
    }
    Angles angles{};
    for(int i=1;i<5;++i) angles[0][i]=u[i];
    angles[1][2]=g1; angles[2][3]=g2; angles[3][4]=g3;
    angles[1][3]=h12; angles[2][4]=h23; angles[1][4]=last;
    Outcome old=physical_test(r,angles);
    if(old!=Split) return old;
    if(newton && !correlated_refine(x,r[2],side)) return Newton;
    return Split;
}

struct Node {Box box; int depth;};
int main(int argc,char**argv) {
    int cap=argc>1?std::atoi(argv[1]):25;
    int requested=argc>2?std::atoi(argv[2]):0;
    int begin=argc>3?std::atoi(argv[3]):18;
    if(cap<0 || cap>35 || std::abs(requested)>1) return 2;
    Box init{Iv{Q(43,400).l,Q(53,400).h},OUTER,{0,S},RING_COS,RING_COS,RING_COS};
    for(int side : {-1,1}) {
        if(requested && side!=requested) continue;
        std::vector<Node> stack{{init,0}};
        std::array<uint64_t,6> counts{};
        uint64_t visited=0;
        while(!stack.empty()) {
            Node b=stack.back(); stack.pop_back();
            ++visited;
            Outcome status=radial_profile(b.box,side,b.depth>=begin);
            if(status==Split && b.depth>=cap) status=Undecided;
            ++counts[status];
            if(status!=Split) continue;
            const std::array<int,6> weights{5,3,2,2,2,2};
            int axis=0;
            for(int j=1;j<6;++j) if((b.box[j].h-b.box[j].l)*weights[j]>(b.box[axis].h-b.box[axis].l)*weights[axis]) axis=j;
            I64 mid=(b.box[axis].l+b.box[axis].h)/2;
            if(mid<=b.box[axis].l || mid>=b.box[axis].h) fail(__LINE__);
            Node l{b.box,b.depth+1},h=l;
            l.box[axis].h=mid; h.box[axis].l=mid;
            stack.push_back(h);stack.push_back(l);
        }
        uint64_t total=0; for(uint64_t n:counts) total+=n;
        if(total!=visited || visited!=2*counts[Split]+1) fail(__LINE__);
        std::cout<<side<<": "<<visited;
        for(uint64_t n:counts) std::cout<<" "<<n;
        std::cout<<std::endl;
    }
}
