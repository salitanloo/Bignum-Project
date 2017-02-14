// Salitanloo
class Bignum{
public:
    static const int BMAX = 1e4+10;
    static const int bmod = 10;
    static const int outbase = 10;
    int d[BMAX],size=0;
    inline void operator +=(const Bignum& b){*this=*this+b;}
    inline void operator *=(const Bignum& b){*this=*this*b;}
    inline void operator -=( Bignum& b){*this=*this-b;}
    inline bool operator <(const Bignum& b)const{return !(*this>b);}
    inline bool operator >=(const Bignum& b)const{return (*this>b)||(*this==b);}
    inline bool operator <=(const Bignum& b)const{return (*this<b)||(*this==b);}
    inline Bignum operator /(const Bignum& b)const{return ((*this).Divide(b)).F;}
    inline Bignum operator %(const Bignum& b)const{return ((*this).Divide(b)).S;}
    inline void operator /=(const Bignum& b){(*this)=(*this)/b;}
    inline void operator %=(const Bignum& b){(*this)=(*this)/b;}
    inline Bignum pow(const Bignum& a,const int& p)const{
        Bignum ans=1;
        Bignum t=a;
        int n=p;
        while(n>0){
            if(n%2==1)
                ans*=t;
            t*=t;
            n>>=1;
        }
        return ans;
    }
    void relax(){
        for(int i=0;i<=size;i++)
            d[i+1]+=(d[i]/bmod),d[i]%=bmod,size+=(d[size+1]>0);
        while(size&&!d[size]&&size--);
    }
    Bignum(int x=0){
        memset(d,0,sizeof d);
        d[0]=x;
        size=0;
        relax();
    }
    void show(){
        for(int i=size;i>=0;i--)
            cout<<d[i];
        relax();
    }
    inline Bignum operator +(const Bignum& b)const{
        Bignum ans(0);
        int n=max(b.size,size);
        ans.size=n;
        for(int i=0;i<=n;i++)
            ans.d[i]=d[i]+b.d[i];
        ans.relax();
        return ans;
    }
    inline Bignum operator *(const Bignum& b)const{
        Bignum ans(0);
        Bignum temp=0;
        if(b==temp||*this==temp)
            return temp;
        ans.size=size*b.size;
        for(int i=0;i<=size;i++)
            for(int j=0;j<=b.size;j++)
                ans.d[i+j]+=d[i]*(b.d[j]);
        ans.relax();
        return ans;
    }
    inline Bignum operator -(const Bignum& t)const{
        Bignum ans=0;
        Bignum a=*this;
        Bignum b=t;
        ans.size=a.size;
        for(int i=0;i<=size;i++){
            if(a.d[i]-b.d[i]<0)
                a.d[i]+=10,b.d[i+1]++;
            ans.d[i]=a.d[i]-b.d[i];
        }
        ans.relax();
        return ans;
    }
    inline bool operator >(const Bignum& b)const{
        if(size<b.size)
            return false;
        if(size>b.size)
            return true;
        for(int i=max(size,b.size);i>=0;i--){
            if(d[i]>b.d[i])
                return true;
            else if(d[i]<b.d[i])
                return false;
        }
        return false;
    }
    inline bool operator ==(const Bignum& b)const{
        for(int i=max(size,b.size);i>=0;i--)
            if(d[i]!=b.d[i])
                return false;
        return true;
    }
    inline Bignum operator *(const int& b)const{
        Bignum ans=0;
        ans.size=size;
        for(int i=0;i<=size;i++)
            ans.d[i]=d[i]*(b);
        ans.relax();
        return ans;
    }
    inline pair<Bignum,Bignum> Divide(const Bignum& b)const{
        Bignum ans=0;
        Bignum rem=0;
        for(int i=size;i>=0;i--){
            rem=rem*10+d[i];
            for(int j=9;j>=0;j--){
                Bignum t=b*j;
                if(t<=rem){
                    rem-=t;
                    ans=ans*10+j;
                    break;
                }
            }
        }
        return {ans,rem};
    }
};
