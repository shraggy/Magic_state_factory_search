"""Replay sweep A's restricted-SAT pass (staged._r_one: check-pattern SAT then synthesis SAT, in one process) on one
class, with the seed sweep A gave it, to see whether the pass itself returns the factory that synthesis SAT run alone
finds. Run from the repository root.
    python outputs/staged_A_rerun/replay_r.py ORBIT SPEC K CAP
"""
import sys, json, time


def main():
    sys.path[:0] = ["searches", "classification"]
    import staged as STG
    o, spec, k, cap = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
    t = time.time(); _, out = STG._r_one((k, o, spec, cap, STG.DEFAULT["iii_r"]))
    print(json.dumps(dict(orbit=o, cap=cap, out=[(kind, b["n"], b["N"]) for kind, b in out],
                          seconds=round(time.time() - t, 1))), flush=True)


if __name__ == "__main__":                     # required: synthesis SAT spawns a subprocess per query (spawn)
    main()
