#!/usr/bin/env python3
"""Independent exact-support oracle on 250 randomized small DAGs (R381)."""
import json, random, struct, subprocess, tempfile, time
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
R357_SRC = ROOT / "reproduction" / "R399" / "src" / "aug_tropical_multiwitness.cpp"
R354_SRC = ROOT / "reproduction" / "R399" / "src" / "exact_best_scenario_support.cpp"
OUT = Path(__file__).resolve().parent / "results"


def allpaths(n, h):
    out = []
    def rec(s, p):
        if s == n:
            out.append(tuple(p)); return
        for t in range(s + h, n + 1):
            rec(t, p + [(s, t)])
    rec(0, [])
    return out


def shortest(n, h, c, beta):
    dp = [None] * (n + 1); pp = [None] * (n + 1)
    dp[0] = 0; pp[0] = []
    for t in range(1, n + 1):
        best = None
        for s in range(t - h + 1):
            if dp[s] is None: continue
            z = dp[s] + c[s, t] + beta
            if best is None or z < best:
                best = z; pp[t] = pp[s] + [(s, t)]
        dp[t] = best
    return tuple(pp[n])


def exact_support(P, paths, L, U, b0, b1):
    E = set(P); kp = len(P); A = sum(L[e] for e in P); Ms = {}
    for Q in paths:
        k = len(Q)
        cost = sum((L[e] if e in E else U[e]) for e in Q)
        Ms[k] = min(Ms.get(k, 10**100), cost)
    if Ms.get(kp, 10**100) != A: return None
    lo, hi = Fraction(b0), Fraction(b1)
    for k, m in Ms.items():
        if k > kp: lo = max(lo, Fraction(A - m, k - kp))
        elif k < kp: hi = min(hi, Fraction(m - A, kp - k))
    return (lo, hi) if lo <= hi else None


def pathstr(P):
    return " ".join(f"{s}:{t}" for s, t in P)


def parse_fraction(s):
    a = s.split("/")
    return Fraction(int(a[0]), int(a[1]) if len(a) > 1 else 1)


def compile_solver(src, out):
    subprocess.run(["g++", "-O3", "-std=c++17", str(src), "-o", str(out)], check=True)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time(); results = []
    with tempfile.TemporaryDirectory(prefix="sac_r381_") as td:
        td = Path(td); r357 = td / "r357"; r354 = td / "r354"
        compile_solver(R357_SRC, r357); compile_solver(R354_SRC, r354)
        random.seed(381)
        count_paths = count_fractional = count_singleton = 0
        for case in range(250):
            n = random.choice([7, 8, 9, 10]); h = random.choice([1, 2, 3])
            ps = allpaths(n, h)
            if len(ps) > 10000: continue
            b0, b1 = 0, random.randint(2, 7)
            L = {}; U = {}
            for s in range(n + 1):
                for t in range(s + h, n + 1):
                    L[s, t] = random.randint(0, 15)
                    U[s, t] = L[s, t] + random.randint(0, 5)
            true = {P: exact_support(P, ps, L, U, b0, b1) for P in ps}
            true = {P: x for P, x in true.items() if x is not None}
            count_paths += len(true)
            count_fractional += sum(1 for lo, hi in true.values() if lo.denominator > 1 or hi.denominator > 1)
            count_singleton += sum(1 for lo, hi in true.values() if lo == hi)
            N = n + 1; arrL = [0] * (N * N); arrU = [0] * (N * N)
            for (s, t), v in L.items(): arrL[s * N + t] = v; arrU[s * N + t] = U[s, t]
            lp = td / "L.bin"; up = td / "U.bin"; wp = td / "w.txt"; cand = td / "cand.txt"; out = td / "support.txt"
            lp.write_bytes(struct.pack("<" + "q" * len(arrL), *arrL))
            up.write_bytes(struct.pack("<" + "q" * len(arrU), *arrU))
            ws = set()
            for scenario in ["L", "M", "U"]:
                c = {e: (L[e] if scenario == "L" else U[e] if scenario == "U" else (L[e] + U[e]) // 2) for e in L}
                for beta in [b0, b1, (b0 + b1) // 2]: ws.add(shortest(n, h, c, beta))
            wp.write_text("\n".join(" ".join(f"{s}:{t}" for s, t in P) for P in sorted(ws)) + "\n")
            run = subprocess.run([str(r357), str(n), str(h), str(b0), str(b1), str(lp), str(up), "1000000", "10", str(wp), "4"], capture_output=True, text=True)
            if "status COMPLETE" not in run.stdout:
                raise RuntimeError(f"R357 did not complete case={case}: {run.stdout[-1000:]} {run.stderr}")
            candpaths = [tuple(tuple(map(int, x.split(":"))) for x in line.split()[1:])
                         for line in run.stdout.splitlines() if line.startswith("PATH ")]
            if not set(true).issubset(set(candpaths)):
                raise RuntimeError(f"R357 missed true candidates case={case}")
            cand.write_text("\n".join("PATH " + pathstr(P) for P in candpaths) + "\n")
            filt = subprocess.run([str(r354), str(n), str(h), str(b0), str(b1), str(lp), str(up), str(cand), str(out)], capture_output=True, text=True)
            if filt.returncode:
                raise RuntimeError(f"R354 failed case={case}: {filt.stdout} {filt.stderr}")
            got = {}
            for line in out.read_text().splitlines():
                if not line.startswith("INDEX "): continue
                z = line.split(); d = dict(zip(z[2::2], z[3::2]))
                if d.get("POSSIBLE") == "1":
                    P = candpaths[int(z[1])]
                    got[P] = (parse_fraction(d["LO"]), parse_fraction(d["HI"]))
            if got != true:
                raise RuntimeError(f"R354 identity/support mismatch case={case}, true={len(true)}, got={len(got)}")
            results.append({"case": case, "n": n, "h": h, "all_paths": len(ps), "possible_paths": len(true),
                            "fractional_endpoint_paths": sum(1 for lo, hi in true.values() if lo.denominator > 1 or hi.denominator > 1),
                            "singleton_paths": sum(1 for lo, hi in true.values() if lo == hi), "status": "PASS"})
            if (case + 1) % 50 == 0: print(f"PASS {case + 1}/250", flush=True)
        summary = {"cases": len(results), "all_pass": len(results) == 250,
                   "possible_path_instances": count_paths, "fractional_endpoint_instances": count_fractional,
                   "singleton_support_instances": count_singleton, "elapsed_seconds": time.time() - t0}
        (OUT / "summary_250.json").write_text(json.dumps({"summary": summary, "cases": results}, indent=2))
        print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
