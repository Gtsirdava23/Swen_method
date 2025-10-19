from math import log2, sqrt
from typing import Callable, List, Optional, Tuple


def f_raw(x: float) -> float:
    arg = x**4 + 2.0 * (x**2) + 2.0
    if arg <= 0 or abs(arg - 1.0) < 1e-15:
        return float("inf")
    return 2.0 - 1.0 / log2(arg)


class FCounter:
 
    def __init__(self, f: Callable[[float], float]):
        self.f = f
        self.count = 0

    def __call__(self, x: float) -> float:
        self.count += 1
        return self.f(x)


def svenn_search(
    f: Callable[[float], float],
    x0: float,
    h0: float,
    N: int = 200,
) -> Tuple[List[float], List[float], float, float, float, float, float, float, float]:
   
    X: List[float] = []
    F: List[float] = []

    def logrow(k: int, h: float, x: float, fx: float) -> None:
        print(f"{k:3d} | h={h:10.6g}  x={x:12.6g}  f={fx: .12g}")

    print("Svenn method — bracketing an uncertainty interval")
    print("  k |        h            x              f(x)")
    print("------------------------------------------------")

    c = float(x0)
    h = float(h0)
    fc = f(c)
    X.append(c)
    F.append(fc)
    logrow(1, h, c, fc)

    b = c + h
    fb = f(b)
    X.append(b)
    F.append(fb)
    logrow(2, h, b, fb)
    k = 2

    if fb < fc:
        a, fa = c, fc
        c, fc = b, fb
        while True:
            h *= 2.0
            b = c + h
            fb = f(b)
            k += 1
            X.append(b)
            F.append(fb)
            logrow(k, h, b, fb)
            if k >= N or fb >= fc:
                break
            a, fa = c, fc
            c, fc = b, fb
    else:
        a = c - h
        fa = f(a)
        k += 1
        X.append(a)
        F.append(fa)
        logrow(k, h, a, fa)
        if fa < fc:
            b, fb = c, fc
            c, fc = a, fa
            while True:
                h *= 2.0
                a = c - h
                fa = f(a)
                k += 1
                X.append(a)
                F.append(fa)
                logrow(k, h, a, fa)
                if k >= N or fa >= fc:
                    break
                b, fb = c, fc
                c, fc = a, fa

    if a > b:
        a, b, fa, fb = b, a, fb, fa

    print(f"\nSvenn result: [a,b]=[{a}, {b}], c={c}")
    return X, F, a, b, fa, fb, c, fc, h


def header_table(title: str) -> None:
    print(f"\n{title}")
    print("  i |        L_i           x_i             f(x_i)")
    print("-------------------------------------------------------")


def final_block(f_calls: int, a: float, b: float, xbest: float, fbest: float) -> None:
    print("\nSUMMARY:")
    print("Function evaluations:", f_calls)
    print("Final interval length:", abs(b - a))
    print("Best point x*:", xbest)
    print("Minimum f(x*):", fbest)


def dichotomy(
    f: Callable[[float], float],
    a: float,
    b: float,
    eps: float,
    delta: Optional[float] = None,
) -> Tuple[float, float, float, float]:
  
    if delta is None:
        delta = eps / 2
    header_table("Dichotomy method")
    i = 0
    while (b - a) / 2 > eps:
        i += 1
        m = (a + b) / 2
        x1 = m - delta
        x2 = m + delta
        f1, f2 = f(x1), f(x2)
        if f1 <= f2:
            b = x2
            x_i, f_i = x1, f1
        else:
            a = x1
            x_i, f_i = x2, f2
        L_i = b - a
        print(f"{i:3d} | {L_i:12.6g}  {x_i:12.6g}  {f_i: .12g}")
    xbest = (a + b) / 2
    return a, b, xbest, f(xbest)


def halve_interval(
    f: Callable[[float], float],
    a: float,
    b: float,
    eps: float,
    delta: Optional[float] = None,
) -> Tuple[float, float, float, float]:
    """Interval halving method."""
    if delta is None:
        delta = eps / 2
    header_table("Interval halving method")
    i = 0
    while (b - a) / 2 > eps:
        i += 1
        m = (a + b) / 2
        fL = f(m - delta)
        fR = f(m + delta)
        if fL < fR:
            b = m + delta
            x_i, f_i = m - delta, fL
        else:
            a = m - delta
            x_i, f_i = m + delta, fR
        L_i = b - a
        print(f"{i:3d} | {L_i:12.6g}  {x_i:12.6g}  {f_i: .12g}")
    xbest = (a + b) / 2
    return a, b, xbest, f(xbest)


def golden_section(
    f: Callable[[float], float], a: float, b: float, eps: float
) -> Tuple[float, float, float, float]:
    """Golden-section search."""
    header_table("Golden-section method")
    phi = (1 + sqrt(5)) / 2
    resphi = 2 - phi  # 1/phi^2
    x1 = b - resphi * (b - a)
    x2 = a + resphi * (b - a)
    f1, f2 = f(x1), f(x2)
    i = 0
    while (b - a) / 2 > eps:
        i += 1
        if f1 <= f2:
            b, x2, f2 = x2, x1, f1
            x1 = b - resphi * (b - a)
            f1 = f(x1)
            x_i, f_i = x1, f1
        else:
            a, x1, f1 = x1, x2, f2
            x2 = a + resphi * (b - a)
            f2 = f(x2)
            x_i, f_i = x2, f2
        L_i = b - a
        print(f"{i:3d} | {L_i:12.6g}  {x_i:12.6g}  {f_i: .12g}")
    xbest = (a + b) / 2
    return a, b, xbest, f(xbest)


def step_adaptation(
    f: Callable[[float], float], a: float, b: float, eps: float
) -> Tuple[float, float, float, float]:

    header_table("Adaptive-step method")
    x = (a + b) / 2
    s = (b - a) / 4
    fx = f(x)
    i = 0
    while s > eps:
        i += 1
        xL, xR = x - s, x + s
        fL, fR = f(xL), f(xR)
        if fL < fx and fL <= fR:
            x, fx = xL, fL
        elif fR < fx and fR < fL:
            x, fx = xR, fR
        else:
            s /= 2
        a_i, b_i = x - s, x + s
        L_i = b_i - a_i
        print(f"{i:3d} | {L_i:12.6g}  {x:12.6g}  {fx: .12g}")
    return x - s, x + s, x, fx


if __name__ == "__main__":
    # Initial guess for Svenn
    x0, h0 = -3.0, 0.4
    eps = 1e-3

    f = FCounter(f_raw)

    Xs, Fs, a, b, fa, fb, c, fc, hf = svenn_search(f, x0, h0, N=200)

    a1, b1, x1, f1 = dichotomy(f, a, b, eps, delta=eps / 2)
    final_block(f.count, a1, b1, x1, f1)

    a2, b2, x2, f2 = halve_interval(f, a, b, eps, delta=eps / 2)
    final_block(f.count, a2, b2, x2, f2)

    a3, b3, x3, f3 = golden_section(f, a, b, eps)
    final_block(f.count, a3, b3, x3, f3)

    a4, b4, x4, f4 = step_adaptation(f, a, b, eps)
    final_block(f.count, a4, b4, x4, f4)
