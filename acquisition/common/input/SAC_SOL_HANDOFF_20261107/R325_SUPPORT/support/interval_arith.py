"""Outward decimal interval arithmetic used for the SAC qualification run.

Contract: every stored IEEE-754 feature, label-independent transform, lambda,
beta and Newton anchor is interpreted as its exact binary value (Decimal.from_float).
The Decimal basic operations are directed at 70 significant decimal digits.
Decimal exp/ln are specified as correctly rounded to nearest; their result is
widened by one representable Decimal on both sides. This is an auditable
software enclosure, not a proof of the Python/libmpdec implementation itself.
"""
from __future__ import annotations

from decimal import Context, Decimal, ROUND_FLOOR, ROUND_CEILING, ROUND_HALF_EVEN
from dataclasses import dataclass
from typing import Iterable, Sequence

PREC = 70
_KW = dict(prec=PREC, Emin=-999999999, Emax=999999999)
DOWN = Context(rounding=ROUND_FLOOR, **_KW)
UP = Context(rounding=ROUND_CEILING, **_KW)
NEAR = Context(rounding=ROUND_HALF_EVEN, **_KW)
ZERO = Decimal(0)
ONE = Decimal(1)
TWO = Decimal(2)
Interval = tuple[Decimal, Decimal]


def dfloat(x: float) -> Decimal:
    """Exact decimal representation of the stored binary64 input."""
    return Decimal.from_float(float(x))


def add_dn(a: Decimal, b: Decimal) -> Decimal:
    return DOWN.add(a, b)


def add_up(a: Decimal, b: Decimal) -> Decimal:
    return UP.add(a, b)


def sub_dn(a: Decimal, b: Decimal) -> Decimal:
    return DOWN.subtract(a, b)


def sub_up(a: Decimal, b: Decimal) -> Decimal:
    return UP.subtract(a, b)


def mul_dn(a: Decimal, b: Decimal) -> Decimal:
    return DOWN.multiply(a, b)


def mul_up(a: Decimal, b: Decimal) -> Decimal:
    return UP.multiply(a, b)


def div_dn(a: Decimal, b: Decimal) -> Decimal:
    if b == 0:
        raise ZeroDivisionError("interval denominator contains zero")
    return DOWN.divide(a, b)


def div_up(a: Decimal, b: Decimal) -> Decimal:
    if b == 0:
        raise ZeroDivisionError("interval denominator contains zero")
    return UP.divide(a, b)


def add_i(a: Interval, b: Interval) -> Interval:
    return add_dn(a[0], b[0]), add_up(a[1], b[1])


def sub_i(a: Interval, b: Interval) -> Interval:
    return sub_dn(a[0], b[1]), sub_up(a[1], b[0])


def scale_i(x: Decimal, a: Interval) -> Interval:
    if x >= 0:
        return mul_dn(x, a[0]), mul_up(x, a[1])
    return mul_dn(x, a[1]), mul_up(x, a[0])


def exp_point(x: Decimal) -> Interval:
    y = NEAR.exp(x)
    if not y.is_finite():
        raise ArithmeticError(f"non-finite Decimal exp at {x}")
    return DOWN.next_minus(y), UP.next_plus(y)


def exp_i(x: Interval) -> Interval:
    return exp_point(x[0])[0], exp_point(x[1])[1]


def ln_i_positive(x: Interval) -> Interval:
    if x[0] <= 0:
        raise ArithmeticError("log interval is not strictly positive")
    lo = NEAR.ln(x[0])
    hi = NEAR.ln(x[1])
    return DOWN.next_minus(lo), UP.next_plus(hi)


def softplus_point(x: Decimal) -> Interval:
    """Enclose log(1+exp(x)) with a stable monotone formulation."""
    if x >= 0:
        # Decimal.__neg__ applies the ambient context (often 28 digits).
        # copy_negate is exact and context-free, so the 70-digit exp enclosure
        # receives the exact sign-flipped endpoint.
        e = exp_point(x.copy_negate())
        arg = add_i((ONE, ONE), e)
        l = ln_i_positive(arg)
        return add_dn(x, l[0]), add_up(x, l[1])
    e = exp_point(x)
    arg = add_i((ONE, ONE), e)
    return ln_i_positive(arg)


def softplus_i(x: Interval) -> Interval:
    return softplus_point(x[0])[0], softplus_point(x[1])[1]


def sigmoid_point(x: Decimal) -> Interval:
    # Keep sign change exact before entering the explicit directed contexts.
    e = exp_point(x.copy_negate())
    denom = add_i((ONE, ONE), e)
    return div_dn(ONE, denom[1]), div_up(ONE, denom[0])


def sigmoid_i(x: Interval) -> Interval:
    return sigmoid_point(x[0])[0], sigmoid_point(x[1])[1]


def dot_interval(x: Sequence[Decimal], theta: Sequence[Decimal]) -> Interval:
    lo = ZERO
    hi = ZERO
    for a, b in zip(x, theta):
        lo = add_dn(lo, mul_dn(a, b))
        hi = add_up(hi, mul_up(a, b))
    return lo, hi


def _ridge_interval(m: int, lam: Decimal, theta: Sequence[Decimal]) -> Interval:
    norm_lo = ZERO
    norm_hi = ZERO
    for t in theta:
        sqlo = mul_dn(t, t)
        sqhi = mul_up(t, t)
        norm_lo = add_dn(norm_lo, sqlo)
        norm_hi = add_up(norm_hi, sqhi)
    coeff_lo = div_dn(mul_dn(Decimal(m), lam), TWO)
    coeff_hi = div_up(mul_up(Decimal(m), lam), TWO)
    return mul_dn(coeff_lo, norm_lo), mul_up(coeff_hi, norm_hi)


@dataclass
class PrefixBank:
    loss_lo: list[Decimal]
    loss_hi: list[Decimal]
    grad_lo: list[list[Decimal]]
    grad_hi: list[list[Decimal]]
    n: int
    p: int
    row_evaluations: int
    exp_evaluations: int
    ln_evaluations: int

    @classmethod
    def build(cls, X: Sequence[Sequence[Decimal]], y: Sequence[int],
              theta_float: Sequence[float]) -> "PrefixBank":
        n = len(X)
        p = len(X[0])
        theta = [dfloat(v) for v in theta_float]
        loss_lo = [ZERO]
        loss_hi = [ZERO]
        grad_lo = [[ZERO] for _ in range(p)]
        grad_hi = [[ZERO] for _ in range(p)]
        exp_calls = 0
        ln_calls = 0
        for i in range(n):
            zlo, zhi = dot_interval(X[i], theta)
            # softplus uses one exp and one log per endpoint.
            sp_lo = softplus_point(zlo)
            sp_hi = softplus_point(zhi)
            exp_calls += 2
            ln_calls += 2
            sp = (sp_lo[0], sp_hi[1])
            if int(y[i]) == 1:
                loss = sub_i(sp, (zlo, zhi))
                ydec = ONE
            else:
                loss = sp
                ydec = ZERO
            loss_lo.append(add_dn(loss_lo[-1], loss[0]))
            loss_hi.append(add_up(loss_hi[-1], loss[1]))

            siglo = sigmoid_point(zlo)
            sighi = sigmoid_point(zhi)
            exp_calls += 2
            sig = (siglo[0], sighi[1])
            residual = sub_i(sig, (ydec, ydec))
            for j in range(p):
                rowg = scale_i(X[i][j], residual)
                grad_lo[j].append(add_dn(grad_lo[j][-1], rowg[0]))
                grad_hi[j].append(add_up(grad_hi[j][-1], rowg[1]))
        # The counters are deterministic operation counts, not FLOPs.
        return cls(loss_lo, loss_hi, grad_lo, grad_hi, n, p, n,
                   exp_calls, ln_calls)

    def query(self, s: int, t: int, lam_float: float,
              theta_float: Sequence[float]) -> tuple[Decimal, Decimal]:
        m = t - s
        lam = dfloat(lam_float)
        theta = [dfloat(v) for v in theta_float]
        loss = (sub_dn(self.loss_lo[t], self.loss_hi[s]),
                sub_up(self.loss_hi[t], self.loss_lo[s]))
        ridge = _ridge_interval(m, lam, theta)
        f = add_i(loss, ridge)
        g: list[Interval] = []
        for j, tj in enumerate(theta):
            part = (sub_dn(self.grad_lo[j][t], self.grad_hi[j][s]),
                    sub_up(self.grad_hi[j][t], self.grad_lo[j][s]))
            rg_lo = mul_dn(mul_dn(Decimal(m), lam), tj)
            rg_hi = mul_up(mul_up(Decimal(m), lam), tj)
            g.append(add_i(part, (rg_lo, rg_hi)))
        norm_sq_hi = ZERO
        for lo, hi in g:
            # Decimal.__abs__ applies the ambient thread context (normally 28
            # digits), which can round a directed gradient endpoint inward.
            # copy_abs is exact and context-free; the square is then rounded
            # upward by the explicit UP context below.
            a = max(lo.copy_abs(), hi.copy_abs())
            norm_sq_hi = add_up(norm_sq_hi, mul_up(a, a))
        mu_lo = mul_dn(Decimal(m), lam)
        denom_lo = mul_dn(TWO, mu_lo)
        residual_penalty_hi = div_up(norm_sq_hi, denom_lo)
        lower = max(ZERO, sub_dn(f[0], residual_penalty_hi))
        upper = f[1]
        if lower > upper:
            raise ArithmeticError(f"invalid segment enclosure [{lower}, {upper}] at {(s,t)}")
        return lower, upper


def direct_segment_bound(X: Sequence[Sequence[Decimal]], y: Sequence[int],
                         s: int, t: int, lam_float: float,
                         theta_float: Sequence[float]) -> tuple[Decimal, Decimal]:
    """One-segment certified bound; no full-window prefix is built."""
    subx = X[s:t]
    suby = y[s:t]
    bank = PrefixBank.build(subx, suby, theta_float)
    return bank.query(0, t - s, lam_float, theta_float)


def dp_add_lower(a: Decimal, b: Decimal, c: Decimal) -> Decimal:
    return add_dn(add_dn(a, b), c)


def dp_add_upper(a: Decimal, b: Decimal, c: Decimal) -> Decimal:
    return add_up(add_up(a, b), c)
