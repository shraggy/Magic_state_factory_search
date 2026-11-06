import sys,os,json,glob; sys.path[:0]=['searches','classification']
from classify import phase_polynomial
from k78_orbits import fingerprint, genuine
import unified_search as us
out={}
for k in (7,8,9):
    rows=[]
    for f in glob.glob(f"outputs/staged_B/l3_k{k}/c*.json"):
        for o,b in json.load(open(f))["best"].items():
            g=[(frozenset(S),c) for S,c in b["gates"]]
            P=phase_polynomial([S for S,_ in g],[c for _,c in g],set(range(k)),3)
            pure=all((len(T)==1 and v%2==0) or (len(T)==2 and v%4==0) or (len(T)==3 and v%4==0) or len(T)>3 or len(T)==0 for T,v in P.items())
            if not pure: continue
            tri=sorted(sorted(T) for T,v in P.items() if len(T)==3 and v%8==4)
            ok=bool(tri) and genuine(tri,k)
            rows.append(dict(label=o,file=os.path.basename(f),n=b["n"],N=b["N"],G=b.get("G"),H=b.get("H"),ccz=tri,genuine=ok,fp=fingerprint(tri,k) if ok and k<=8 else None))
    out[k]=rows
    print(k,"B outputs that are CCZ layers:",len(rows),"genuine:",sum(r["genuine"] for r in rows),"n:",sorted((r["n"],r["N"]) for r in rows if r["genuine"]))
json.dump(out,open("outputs/k78_pipeline/staged_B_ccz.json","w"))
