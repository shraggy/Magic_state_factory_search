"""Check that Campbell-Howard Lemma 1 (PRA 95, 022316, Eq. 23: |K^T x + S^T y| ~c F(x) mod 8) agrees with
the factory condition Eq. (4) of the paper at l=3 with unit coefficients, on random circuits (seed 1).
Run: nice -n 19 python3 -I searches/ch_lemma1_check.py  -> outputs/checks/ch_lemma1_check.json"""
import itertools, random
random.seed(1)
def mobius(f, N):  # integer multilinear coefficients of f:{0,1}^N->Z (mod 8)
    a = dict(f)
    for i in range(N):
        for m in range(1<<N):
            if m>>i&1: a[m]=(a[m]-a[m^(1<<i)])%8
    return a
def ch(cols,k,N):
    # P(x,y)-P(x,0) must be ~c 0: linear coeff even, quadratic = 0 mod 4, cubic+ = 0 mod 8
    P={z:sum(bin(z&c).count('1')&1 for c in cols)%8 for z in range(1<<N)}
    D={z:(P[z]-P[z&((1<<k)-1)])%8 for z in range(1<<N)}
    a=mobius(D,N)
    for m,v in a.items():
        w=bin(m).count('1')
        if w==1 and v%2: return False
        if w==2 and v%4: return False
        if w>=3 and v%8: return False
    return True
def eq4(cols,k,N):
    for t in (1,2,3):
        for Q in itertools.combinations(range(N),t):
            m=sum(1<<q for q in Q)
            if m < (1<<k) : continue   # Q inside outputs (outputs are qubits 0..k-1)
            if sum(1 for c in cols if c&m==m)%2: return False
    return True
agree=tot=pos=0
for trial in range(20000):
    k=random.randint(1,3); s=random.randint(1,3); N=k+s
    pool=[c for c in range(1,1<<N) if c>>k]  # every rotation touches a check
    n=random.randint(1,len(pool)); cols=random.sample(pool,n)
    if trial%2==0:  # bias toward solutions: add symmetric-difference fixups not needed; also test known 832
        pass
    a,b=ch(cols,k,N),eq4(cols,k,N); tot+=1; agree+=(a==b); pos+=a
# [[8,3,2]]: check q3 (index 3), outputs 0..2, all 8 columns containing check
c832=[ (1<<3)|m for m in range(8)]
r832=(ch(c832,3,4), eq4(c832,3,4)); print("832", *r832)
print("agree",agree,"of",tot,"; CH-true",pos)
import json, os
out = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'outputs', 'checks', 'ch_lemma1_check.json')
json.dump({'trials': tot, 'agree': agree, 'ch_lemma1_true': pos, 'seed': 1,
           'factory_832_ch': r832[0], 'factory_832_eq4': r832[1],
           'note': 'k=1..3 outputs, s=1..3 checks, random sets of rotations that each touch a check'},
          open(out, 'w'), indent=1)
