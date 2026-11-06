# CS12CS13 as parity rotations (units pi/4) on 3 qubits, add one T on parity {1,2}; list odd-coefficient parities and their F2 rank.
import itertools
def table(f): return [f(x)%8 for x in itertools.product((0,1),repeat=3)]
cs=table(lambda x: 2*(x[0]*x[1]+x[0]*x[2]))
def par(c): # parity-rotation coefficients c_L with phase sum_L c_L [L.x] reproducing table mod 8 up to diagonal Clifford? exact: solve via Walsh
    import fractions
    out={}
    for L in itertools.product((0,1),repeat=3):
        if not any(L): continue
        s=sum(((-1)**(sum(a*b for a,b in zip(L,x))))*v for x,v in zip(itertools.product((0,1),repeat=3),cs_cur))
        out[L]=fractions.Fraction(-s,4)  # phase = c0 + sum c_L*[L.x]; c_L = -(1/4)*sum (-1)^{L.x} f(x) for 3 qubits
    return out
for add in [None]+[L for L in itertools.product((0,1),repeat=3) if any(L)]:
    cs_cur=[(v+(1 if add and sum(a*b for a,b in zip(add,x))%2 else 0))%8 for x,v in zip(itertools.product((0,1),repeat=3),cs)]
    c=par(None)
    odd=[L for L,v in c.items() if v.denominator==1 and v.numerator%2]
    nonint=[L for L,v in c.items() if v.denominator!=1]
    # rank
    rows=[int(''.join(map(str,L)),2) for L in odd]; basis=[]
    for v in rows:
        for b in basis: v=min(v,v^b)
        if v: basis.append(v)
    print('add T on',add,'odd parities',odd,'count',len(odd),'rank',len(basis), 'nonint',nonint)
