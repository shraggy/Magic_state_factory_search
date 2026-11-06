"""chain_iii3.py plus stage kind x2:W:C = free XOR (x:W:C) run a second time later in the chain under its own name, so
pass files and per-stage gains stay distinct (merge_A_v3.py names pass files by stage). Nothing else changes.
    python searches/chain_iii4.py k I M SEEDDIR OUT CHAIN [--classes o1,o2]
    e.g. CHAIN xs:None:2e6,xs:6:2e6,x:5:2e6,x:6:2e6,x:None:2e6,xsb:8:None:2e6,x2:5:2e6,x2:6:2e6,x2:None:2e6
"""
import os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import chain_iii3 as C3
import chain_iii2 as C2
_s3 = C3.stage


def stage(name, k, spec, cap):
    if name.startswith("x2:"):
        return _s3("x:" + name[3:], k, spec, cap)
    return _s3(name, k, spec, cap)


C2.stage = stage

if __name__ == "__main__":
    a = sys.argv[1:]; only = None
    if "--classes" in a:
        i = a.index("--classes"); only = a[i + 1].split(","); a = a[:i] + a[i + 2:]
    C2.run(int(a[0]), int(a[1]), int(a[2]), a[3], a[4], a[5].split(","), only)
