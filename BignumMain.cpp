class Bignum{
public:
    static const int BMAX = 1e4+10;
    static const int bmod = 10;
    static const int outbase = 10;
    int d[BMAX],size=0;
    inline void operator +=(const Bignum& b){*this=*this+b;}
    inline void operator *=(const Bignum& b){*this=*this*b;}
    inline void operator -=( Bignum& b){*this=*this-b;}
    inline bool operator <(const Bignum& b)const{return !((*this>b)||(*this==b));}
    inline bool operator >=(const Bignum& b)const{return (*this>b)||(*this==b);}
    inline bool operator <=(const Bignum& b)const{return (*this<b)||(*this==b);}
    inline Bignum operator /(const Bignum& b)const{return ((*this).Divide(b)).F;}
    inline Bignum operator %(const Bignum& b)const{return ((*this).Divide(b)).S;}
    inline void operator /=(const Bignum& b){(*this)=(*this)/b;}
    inline void operator %=(const Bignum& b){(*this)=(*this)/b;}
    inline Bignum pow(const int& n)const{
        if(n==0)
            return 1;
        else if(n==1)
            return *this;
        else if(n%2==0)
            return ((*this)*(*this)).pow(n/2);
        else
            return (*this)*((*this)*(*this)).pow(n/2);
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
    inline Bignum sqrt()const{
        Bignum i(0);
        Bignum j=*this;
        while(i+1<j){
            Bignum mid=(i+j)/2;
            Bignum t=mid*mid;
            if(t<=*this)
                i=mid;
            else
                j=mid-1;
        }
        if(j*j<=*this)
            return j;
        return i;
    }
    inline Bignum ssqrt()const{
        int ans[BMAX],n[BMAX],a[BMAX],temp[BMAX],cnt,t;
        Bignum anss=0;
        memset(n,0,sizeof n);
        memset(ans,0,sizeof ans);
        memset(a,0,sizeof a);
        memset(temp,0,sizeof temp);
        for(int i=0;i<=size;i++)
            n[i+1]=d[size-i];
        cnt=size+1;
        if(cnt&1)
            t=0,a[1]=n[1];
        else
            t=1,a[1]=n[1]*10+n[2];
        for(int i=1;i<=(cnt-1)>>1;i++)
            a[i+1]=n[i*2+t]*10+n[i*2+t+1];
        cnt=(cnt+1)>>1;
        reverse(a+1,a+cnt+1);
        for(int i=cnt;i>=1;i--){
            //cout<<42<<endl;
            for(int j=i+1;j<=cnt;j++){
                int c=j-(i+1);
                if(c&1)
                    temp[i+(c/2)]+=ans[j]*200;
                else
                    temp[i+(c/2)]+=ans[j]*20;
            }
            //cout<<50<<endl;
            for(int j=i+1;j<=cnt+1;j++)
                temp[j]+=temp[j-1]/100,temp[j-1]%=100;
            temp[i]+=10;
            for(int j=9;j>=0;j--){
                fill(n,n+BMAX,0);
                temp[i]--;
                for(int k=i;k<=cnt;k++){
                    n[k]+=temp[k]*j;
                    n[k+1]=n[k]/100;
                    n[k]%=100;
                }
                bool flag=1;
                for(int k=cnt+1;k>=i;k--){
                    if(n[k]>a[k]){
                        flag=0;
                        break;
                    }
                    else if(n[k]<a[k])
                        break;
                }
                if(flag){
                    ans[i]=j;
                    anss=anss*10+j;
                    for(int k=i;k<=cnt;k++){
                        a[k]-=n[k];
                        while(a[k]<0)
                            a[k+1]--,a[k]+=100;
                    }
                    break;
                }
            }
            memset(temp,0,sizeof temp);
        }
        return anss;
    }
};
