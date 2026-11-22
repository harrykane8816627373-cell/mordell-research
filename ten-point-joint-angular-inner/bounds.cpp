#include <algorithm>
#include <array>
#include <cstdint>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

using I = std::int64_t;
constexpr I S = 1LL << 20;
struct Interval { I lo,hi; };
I floor_div(I a,I b) { if(b<=0)throw std::runtime_error("denominator"); return a>=0?a/b:-((-a+b-1)/b); }
I ceil_div(I a,I b) { return -floor_div(-a,b); }
Interval rat(I a,I b=1) { return {floor_div(a*S,b),ceil_div(a*S,b)}; }
Interval operator+(Interval a,Interval b){return {a.lo+b.lo,a.hi+b.hi};}
Interval operator-(Interval a){return {-a.hi,-a.lo};}
Interval operator-(Interval a,Interval b){return a+(-b);}
I product(I a,I b) {
    if (std::abs(a)>(1LL<<30)||std::abs(b)>(1LL<<30)) throw std::runtime_error("product range");
    return a*b;
}
Interval operator*(Interval a,Interval b) {
    const I terms[]={product(a.lo,b.lo),product(a.hi,b.lo),product(a.lo,b.hi),product(a.hi,b.hi)};
    return {floor_div(*std::min_element(terms,terms+4),S),ceil_div(*std::max_element(terms,terms+4),S)};
}
Interval square(Interval a) {
    auto v=a*a;
    return {std::max<I>(0,v.lo),v.hi};
}
Interval inv(Interval a) { if(a.lo<=0)throw std::runtime_error("positive inversion");
    return {floor_div(S*S,a.hi),ceil_div(S*S,a.lo)}; }
struct Z {Interval re,im;};
Z operator+(Z a,Z b){return {a.re+b.re,a.im+b.im};}
Z operator-(Z a,Z b){return {a.re-b.re,a.im-b.im};}
Z operator*(Z a,Z b){return {a.re*b.re-a.im*b.im,a.im*b.re+a.re*b.im};}
Z operator*(Z a,Interval x){return {a.re*x,a.im*x};}
Z scalar(I n){return {rat(n),rat(0)};}
Z sq(Z a){return {square(a.re)-square(a.im),a.re*a.im*rat(2)};}
Z conj(Z a){return {a.re,-a.im};}
Interval cell(int l,int h,int denom,int count,int i){
    return {rat(l*count+(h-l)*i,denom*count).lo,rat(l*count+(h-l)*(i+1),denom*count).hi};
}
int main() {
 try {
    constexpr int num_x=12, num_y=16, num_r=32, num_w=8;
    int n,lower,upper;
    if(!(std::cin>>n>>lower>>upper))throw std::runtime_error("phase cells");
    struct Phase {Z s,s2,s4;};
    std::vector<Phase> phases(n);
    for(auto &p:phases) for(Z *v:{&p.s,&p.s2,&p.s4}){
      if(!(std::cin>>v->re.lo>>v->re.hi>>v->im.lo>>v->im.hi))throw std::runtime_error("phase table");
    }
    // Bin contributions use one midpoint cell throughout each cyclic walk.
    I total=-10000*S;
    long long eval=0;
    for (int ix=0;ix<num_x;ix++)for(int iy=0;iy<num_y;iy++) {
        Z m={cell(0,57,1000,num_x,ix),cell(0,88,1000,num_y,iy)};
        if((square(m.re)+square(m.im)).lo>square(rat(88,1000)).hi)continue;
        std::vector<I> caps(n,-10000*S);
        for(int ia=0;ia<n;ia++)for(int ir=0;ir<num_r;ir++) {
            Interval r=cell(940,1130,1000,num_r,ir),u=square(r),ri=inv(r),ui=square(ri);
            const auto ph=phases[ia];
            Z s=ph.s*ri,s2=ph.s2*ui,s4=ph.s4*square(ui);
            Z A=m*s,T=scalar(1)-A,T2=sq(T);
            for(int iw=0;iw<num_w;iw++) {
                const Interval w=cell(107,133,1000,num_w,iw);
                // E = (1-m/z)^2-w/z^2. These three entries are the second
                // derivatives of the sum of the two logarithmic inner-outer edges.
                Z B=s2*w,E=T2-B;
                Interval U=inv(square(E.re)+square(E.im));
                Z J=sq(conj(E)*U);
                Interval H=( ((scalar(2)-A)*B+A*T2)*J ).re*rat(4);
                Interval cross=(s2*T*J).im*rat(-4);
                Interval ww=(s4*J).re*rat(-2);
                Interval gap=rat(547,100)-rat(1,2)-H;
                const Interval cost=ww+square(cross)*inv(gap)+(u-rat(109,100))*rat(8);
                caps[ia]=std::max(caps[ia],cost.hi);
                eval++;
            }
        }
        // Every ordered octagon starts in a bin below upper after a cyclic
        // relabeling. Enforce the eighth edge back to this initial bin.
        const I none=-1000000*S;
        for(int start=0;start<upper;start++) {
            std::vector<I> prev(start+n+1,none);
            prev[start]=caps[start];
            for(int edge=1;edge<8;edge++) {
                std::vector<I> next(prev.size(),none);
                for(int k=start+edge*lower;k<=std::min<int>(start+n-lower,start+edge*upper);k++) {
                    I best=none;
                    for(int t=std::max(start,k-upper);t<=k-lower;t++)best=std::max(best,prev[t]);
                    if(best>none)next[k]=best+caps[k%n];
                }
                next.swap(prev);
            }
            for(int k=start+n-upper;k<=start+n-lower;k++)total=std::max(total,prev[k]);
        }
    }
    std::cout<<eval<<" "<<total<<" "<<S<<"\n";
 }catch(const std::exception &s){std::cerr<<s.what()<<"\n";return 1;}
}
