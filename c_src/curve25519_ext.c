#define FOR(i,n) for (i = 0;i < n;++i)
#define sv static void

#include <stdio.h> // DEBUG

typedef unsigned char u8;
typedef long long i64;
typedef i64 gf[16];

static const gf _121665 = {0xDB41,1};

sv car25519(gf o)
{
  int i;
  i64 c;
  FOR(i,16) {
    o[i]+=(1LL<<16);
    c=o[i]>>16;
    o[(i+1)*(i<15)]+=c-1+37*(c-1)*(i==15);
    o[i]-=c<<16;
  }
}

sv sel25519(gf p,gf q,int b)
{
  i64 t,i,c=~(b-1);
  FOR(i,16) {
    t= c&(p[i]^q[i]);
    p[i]^=t;
    q[i]^=t;
  }
}

sv pack25519(u8 *o,const gf n)
{
  int i,j,b;
  gf m,t;

  FOR(i,16) t[i]=n[i];
  car25519(t);
  car25519(t);
  car25519(t);
  FOR(j,2) {
    m[0]=t[0]-0xffed;
    for(i=1;i<15;i++) {
      m[i]=t[i]-0xffff-((m[i-1]>>16)&1);
      m[i-1]&=0xffff;
    }
    m[15]=t[15]-0x7fff-((m[14]>>16)&1);
    b=(m[15]>>16)&1;
    m[14]&=0xffff;
    sel25519(t,m,1-b);
  }
  FOR(i,16) {
    o[2*i]=t[i]&0xff;
    o[2*i+1]=t[i]>>8;
  }
}

sv unpack25519(gf o, const u8 *n)
{
  int i;
  FOR(i,16) o[i]=n[2*i]+((i64)n[2*i+1]<<8);
  o[15]&=0x7fff;
}

sv A(gf o,const gf a,const gf b)
{
  int i;
  FOR(i,16) o[i]=a[i]+b[i];
}

sv Z(gf o,const gf a,const gf b)
{
  int i;
  FOR(i,16) o[i]=a[i]-b[i];
}

sv M(gf o,const gf a,const gf b)
{
  i64 i,j,t[31];
  FOR(i,31) t[i]=0;
  FOR(i,16) FOR(j,16) t[i+j]+=a[i]*b[j];
  FOR(i,15) t[i]+=38*t[i+16];
  FOR(i,16) o[i]=t[i];
  car25519(o);
  car25519(o);
}

sv S(gf o,const gf a)
{
  M(o,a,a);
}

sv inv25519(gf o,const gf i)
{
  gf c;
  int a;

  FOR(a,16) c[a]=i[a];
  for(a=253;a>=0;a--) {
    S(c,c);
    if(a!=2&&a!=4) M(c,c,i);
  }
  FOR(a,16) o[a]=c[a];
}


/*
sv dump(const char* msg, i64* x, i64* z) {
    gf tmp;
    u8 r[32];
    int i;

    printf("%s ", msg);
    inv25519(tmp, z);
    M(tmp, tmp, x);
    pack25519(r, tmp);
    for (i=0; i<32; i++) printf(",%d", r[i]);
    printf("\r\n");
}

sv dump(const char* msg, gf x) {
    u8 r[32];
    int i;

    printf("%s ", msg);
    pack25519(r, x);
    for (i=0; i<32; i++) printf(",%d", r[i]);
    printf("\r\n");
}
*/

/*
sv normalize(i64* x, i64* z) {
    gf tmp;
    u8 r[32];
    int i;

    // Infinity check:
    {
        i64 chk=0;
        FOR (i,16) chk |= z[i];
        if (chk==0) return;
    }
    inv25519(tmp, z);
    M(tmp, tmp, x);
    pack25519(r, tmp);
    unpack25519(x, r);
    FOR (i,16) z[i]=0;
    z[0]=1;
}
*/

int curve25519_ext_scalarmult_unrestrained(u8 *q,const u8 *n,const u8 *p)
{
  u8 z[32];
  i64 x[80],r,i;
  gf a,b,c,d,e,f;
  FOR(i,32) z[i]=n[i];

  // -=-Restraining removed-=-
  unpack25519(x,p);
  FOR(i,16) {
    b[i]=x[i];
    d[i]=a[i]=c[i]=0;
  }
  a[0]=d[0]=1;

  for(i=254;i>=0;--i) {
    r=(z[i>>3]>>(i&7))&1;
    sel25519(a,b,r);
    sel25519(c,d,r);
    A(e,a,c);
    Z(a,a,c);
    A(c,b,d);
    Z(b,b,d);
    S(d,e);
    S(f,a);
    M(a,c,a);
    M(c,b,e);
    A(e,a,c);
    Z(a,a,c);
    S(b,a);
    Z(c,d,f);

    M(a,c,_121665);
    A(a,a,d);
    M(c,c,a);
    M(a,d,f);
    M(d,b,x);
    S(b,e);
    sel25519(a,b,r);
    sel25519(c,d,r);
  }
  FOR(i,16) {
    x[i+16]=a[i];
    x[i+32]=c[i];
    x[i+48]=b[i];
    x[i+64]=d[i];
  }
  inv25519(x+32,x+32);
  M(x+16,x+16,x+32);
  pack25519(q,x+16);
  return 0;
}


#define SQRT_DIVISOR_COUNT 3
struct sqrt_divisor_table_entry {
    u8 b_squared[32];
    gf b_inverse;
};
static const struct sqrt_divisor_table_entry sqrt_divisor_table[SQRT_DIVISOR_COUNT] = {
    {{0}, {0}},
    {{1}, {1}},
    {{0xec, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0xff, 0x7f}, // Prime-1
    {0x5f3d, 0xb5f1, 0xe4d8, 0x3b11, 0x1b87, 0x52d0, 0xe7f9, 0xd0bc, 0x2858, 0xc204, 0xff66, 0xd4b2, 0x20f4, 0xb03e, 0xdb7f, 0x547c
    }} // 38214883241950591754978413199355411911188925816896391856984770930832735035197
};

/* Compute x^((Q+1)/2) and x^Q (mod p)
 * where Q = (p-1)/4
 * = ((1**255-19)-1)/4
 * = (1**255-20)/4
 * = 1**253-5
 * ((Q+1)/2) = (1**253-4)/2 = (1**252-2)
 */
static void sqrt_exps(gf r, gf t, const gf x) {
    // Exponent bit masks:
    // Q:         All in [0;252] except #2.
    // ((Q+1)/2): All in [0;251] except #0.
    // Strategy: First compute x^((Q-1)/2).
    // ((Q-1)/2): All in [0;251] except #1.

    gf c;
    int a;
    FOR(a,16) c[a]=x[a];
    for(a=250;a>=0;a--) {
        S(c,c);
        if(a!=1) M(c,c,x);
    }
    {
        gf b;
        M(b, c, x); // Compute x^((Q+1)/2).
        FOR(a,16) r[a]=b[a];
        M(b, b, c); // Compute x^Q.
        FOR(a,16) t[a]=b[a];
    }
}

static int eq(const u8* a, const u8* b)
{
    u8 neq = 0;
    int i;
    FOR(i,32) {
        neq |= (a[i]-b[i]);
    }
    neq |= (neq>>4);
    neq |= (neq>>2);
    neq |= (neq>>1);
    return (~neq)&1;
}

int mod25519_sqrt(u8 *sqrt_x, const u8 *x) {
    //     SqrtExp = (?CURVE_Q+1) div 2,
    // SqrtCand = powmod(X, SqrtExp, ?PRIME25519),
    // Check = powmod(X, ?CURVE_Q, ?PRIME25519),
    // InvDivisor = sqrt_divisor(Check),
    // (InvDivisor * SqrtCand) rem ?PRIME25519.
    int i;
    gf r,t;
    gf a;
    unpack25519(a,x);
    sqrt_exps(r, t, a);

    int good = 0;
    {
        gf b_inv = {0};
        u8 t_packed[32];
        pack25519(t_packed, t);

        FOR(i,SQRT_DIVISOR_COUNT) {
            int match = eq(t_packed, sqrt_divisor_table[i].b_squared);
            gf cand_inv;
            int j;
            FOR(j,16) cand_inv[j] = sqrt_divisor_table[i].b_inverse[j];
            sel25519(b_inv, cand_inv, match);
            good |= match;
        }
        M(r, r, b_inv);
    }
    pack25519(sqrt_x, r);
    return !good;
}

/*
int curve25519_recover_y(const u8 *y, const u8 *x) {
}
*/
