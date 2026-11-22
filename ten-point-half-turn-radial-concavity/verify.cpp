#include <algorithm>
#include <array>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <stdexcept>

using I = std::int64_t;
using Z = __int128;
constexpr int N=2048, D=65536, R=10000;
constexpr int K=D/4, MASS=43000;
struct Range {int l,u;};
struct Box {std::array<Range,4> angle; std::array<Range,5> radius;};
using Mat=std::array<std::array<I,5>,5>;
std::array<Range,N+1> cs;
std::int64_t nodes, geometry, accepted, unresolved;
int deepest;

I ceil_ratio(Z n,Z d) {
    if (d<=0) throw std::logic_error("nonpositive denominator");
    return I(n>=0 ? (n+d-1)/d : n/d);
}
Range cosine(Range a) {  // All used intervals are in [0,N].
    int l=std::min(cs[a.l].l,cs[a.u].l), u=std::max(cs[a.l].u,cs[a.u].u);
    if (a.l<=N/2 && N/2<=a.u) l=-D;
    return {l,u};
}
int b_at(int x,int y,int c) {
    Z xy=Z(x)*y, ss=Z(x)*x+Z(y)*y, den=ss*D-2*xy*c;
    return int(ceil_ratio(4*xy*(2*xy*D-ss*c)*D*D,den*den));
}
int stiffness(Range x,Range y,int c) {
    int k;
    if (x.u<=y.l) k=b_at(x.u,y.l,c);
    else if (y.u<=x.l) k=b_at(x.l,y.u,c);
    else k=b_at(1,1,c);
    return std::max({k,b_at(x.l,y.u,c),b_at(x.u,y.l,c)});
}
// Lower enclosure of the Schur complement by subtracting an upper rank-one
// product using directed rounding; diagonal losses dominate off-diagonal errors.
bool definite(Mat a) {
    for (int i=0;i<5;i++) {
        // Necessary two-by-two tests keep every accepted pivot and Schur
        // update inside the original diagonal scale.
        if (a[i][i]<=0) return false;
        for (int j=i+1;j<5;j++) {
            if (a[j][j]<=0 || Z(a[i][j])*a[i][j]>=Z(a[i][i])*a[j][j]) return false;
        }
        for (int j=i+1;j<5;j++)
            for (int k=j;k<5;k++) {
                I q=ceil_ratio(Z(a[i][j])*a[i][k],a[i][i]);
                a[j][k]-=q;
            }
        for (int j=i+1;j<5;j++) a[j][j]-=(4-i);
    }
    return true;
}
bool admissible(const Box& b,std::array<std::array<int,5>,5>& c) {
    const auto& t=b.angle;
    // adjacent doubled gap cosine is in (-.3038,.2168)
    constexpr int cl=(-1519*D)/5000-1, cu=(271*D+1249)/1250;
    for (int i=0;i<4;i++) {
        c[0][i+1]=cosine(t[i]).l;
        for (int j=0;j<i;j++) {
            Range gap{t[i].l-t[j].u,t[i].u-t[j].l};
            if (gap.l<0 || gap.u>N) throw std::logic_error("angular order");
            auto x=cosine(gap);
            if (i-j!=2) {
                if (x.u<cl || x.l>cu) return false;
                c[j+1][i+1]=std::min(cu,x.u);
            } else {
                constexpr int opposite=-(81*D)/100; // upward rounding
                if (x.l>opposite) return false;
                c[j+1][i+1]=std::min(opposite,x.u);
            }
        }
    }
    return true;
}
bool vertex(const Mat& m,const std::array<int,4>& r) {
    Mat a=m;
    for (int i=0;i<4;i++) a[i+1][i+1]+= I(10)*D*r[i]/R;
    return definite(a);
}
bool all_radii(const Box& b,const std::array<std::array<int,5>,5>& c) {
    Mat m{};
    for (int i=0;i<5;i++) {
        m[i][i]-=K;
        for (int j=0;j<i;j++) {
            int k=stiffness(b.radius[i],b.radius[j],c[j][i]);
            m[i][i]-=k; m[j][j]-=k; m[j][i]=k;
        }
    }
    m[0][0]+=I(10)*D*b.radius[0].l/R;
    std::array<int,4> r;
    int minsum=0;
    for (int i=0;i<4;i++) { r[i]=b.radius[i+1].l; minsum+=r[i]; }
    if (vertex(m,r)) return true;
    if (minsum>=MASS) return false;
    for (int j=0;j<4;j++) {
        for (int mask=0;mask<8;mask++) {
            int l=0,sum=0;
            for (int i=0;i<4;i++) if (i!=j) {
                auto h=b.radius[i+1]; r[i]=(mask & (1<<l++))?h.u:h.l; sum+=r[i];
            }
            r[j]=MASS-sum;
            if (r[j]>=b.radius[j+1].l && r[j]<=b.radius[j+1].u && !vertex(m,r)) return false;
        }
    }
    return true;
}
void visit(Box b,int depth) {
    nodes++; deepest=std::max(deepest,depth);
    int sum=0; for (int i=1;i<5;i++) sum+=b.radius[i].u;
    std::array<std::array<int,5>,5> c{};
    if (sum<MASS || !admissible(b,c)) {geometry++;return;}
    for (int i=1;i<5;i++) b.radius[i].l=std::max(b.radius[i].l,MASS-sum+b.radius[i].u);
    if (all_radii(b,c)) {accepted++;return;}
    std::array<Range*,9> p;
    for (int i=0;i<4;i++) p[i]=&b.angle[i];
    for (int i=0;i<5;i++) p[i+4]=&b.radius[i];
    int best=-1,score=0;
    for (int i=0;i<9;i++) {
        int width=p[i]->u-p[i]->l;
        int v=width*(i<4?60:(i==4?3:1));
        if (width>=2 && v>score) { best=i;score=v; }
    }
    if (best<0) {unresolved++;return;}
    Range h=*p[best]; int mid=(h.l+h.u)/2;
    p[best]->u=mid;visit(b,depth+1);
    p[best]->l=mid;p[best]->u=h.u;visit(b,depth+1);
}
int main(int argc,char** argv) {
    if (argc!=2) return 2;
    std::ifstream f(argv[1]);
    for (auto& a:cs) if (!(f>>a.l>>a.u)) return 2;
    Box b; b.radius[0]={1075,1325};
    for (int i=1;i<5;i++) b.radius[i]={8900,12710};

    constexpr int s=64, bins=N/s;
    // Intervals have ordered, separated phases (the certified minimum gap
    // exceeds two bins). A point on the angular cut is assigned angle zero.
    for (int i=0;i<bins;i++) for(int j=i+1;j<bins;j++)
    for(int k=j+1;k<bins;k++) for(int l=k+1;l<bins;l++) {
        std::array<int,4> id{i,j,k,l};
        for (int a=0;a<4;a++) b.angle[a]={s*id[a],s*(id[a]+1)};
        visit(b,0);
    }
    std::cout<<nodes<<" "<<geometry<<" "<<accepted<<" "<<unresolved<<" "<<deepest<<"\n";
    return unresolved!=0;
}
