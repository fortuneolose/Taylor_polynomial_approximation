"""
simulate_taylor.py
Bit-accurate Python model of taylor_poly.v + the taylor_poly_tb.v test vectors.

Mirrors every fixed-point operation in the RTL:
  - Signed Q(IW.FRAC_BITS) arithmetic  (default Q15.16, scale=65536)
  - Full 2*DATA_WIDTH multiply followed by truncation of FRAC_BITS LSBs
  - ORDER+1 pipeline stages (Horner's method, one register per stage)
  - active-low async reset, valid handshake

Run with:
    python simulate_taylor.py
"""

import sys

# ---------------------------------------------------------------------------
# Fixed-point helpers
# ---------------------------------------------------------------------------

def to_signed(v, bits):
    """Reinterpret an unsigned integer as a signed two's-complement integer."""
    if v >= (1 << (bits - 1)):
        v -= (1 << bits)
    return v


def to_unsigned(v, bits):
    """Mask a (possibly negative) integer to its unsigned bits-wide value."""
    return v & ((1 << bits) - 1)


def fp_mul_trunc(a, b, data_width, frac_bits):
    """
    Signed multiply of two Q(IW.frac_bits) values and truncate back.

    Matches the RTL:
        prod_full  = a * b;                              // 2*data_width bits
        prod_trunc = prod_full[data_width+frac_bits-1 : frac_bits]
    """
    product = a * b                                      # exact signed product
    # Arithmetic right-shift by frac_bits then mask to data_width bits
    shifted = product >> frac_bits
    return to_signed(to_unsigned(shifted, data_width), data_width)


def real_to_fp(r, frac_bits, data_width):
    """Convert a Python float to signed Q fixed-point integer."""
    raw = int(round(r * (1 << frac_bits)))
    return to_signed(to_unsigned(raw, data_width), data_width)


def fp_to_real(v, frac_bits):
    """Convert a signed fixed-point integer to float."""
    return v / (1 << frac_bits)


# ---------------------------------------------------------------------------
# Taylor polynomial evaluator – combinational Horner step
# (called once per clock cycle per stage)
# ---------------------------------------------------------------------------

class TaylorPolyPipeline:
    """
    RTL-accurate pipeline model of taylor_poly.v.

    Parameters
    ----------
    coeffs   : list of signed fixed-point integers, coeffs[k] = c_k
               (length must be >= order+1)
    order    : polynomial degree
    x0       : expansion point in fixed-point
    data_width, frac_bits : fixed-point format
    """

    def __init__(self, coeffs, order, x0=0, data_width=32, frac_bits=16):
        self.coeffs     = coeffs
        self.order      = order
        self.x0         = x0
        self.DW         = data_width
        self.FB         = frac_bits
        self.latency    = order + 1

        # Pipeline registers: dx_r[i], acc_r[i], vld_r[i]  for i = 0..order
        self.dx_r  = [0] * (order + 1)
        self.acc_r = [0] * (order + 1)
        self.vld_r = [False] * (order + 1)

    def reset(self):
        self.dx_r  = [0] * (self.order + 1)
        self.acc_r = [0] * (self.order + 1)
        self.vld_r = [False] * (self.order + 1)

    def clock(self, valid_in, x_in):
        """
        Simulate one rising clock edge.
        Returns (valid_out, y_out).
        """
        DW = self.DW
        FB = self.FB
        N  = self.order

        # ── Stage 0 (combinational inputs → stage-0 registers) ────────────────
        # Must compute new stage-0 values BEFORE overwriting, since later
        # stages read previous-cycle values of lower-numbered stages.
        # We therefore compute all next-state values first, then commit.

        # Compute next state for stage 0
        dx_next_0  = to_signed(to_unsigned(x_in - self.x0, DW), DW)
        acc_next_0 = self.coeffs[N]   # seed with c_ORDER
        vld_next_0 = valid_in

        # Compute next state for stages 1..ORDER
        dx_next  = [None] * (N + 1)
        acc_next = [None] * (N + 1)
        vld_next = [None] * (N + 1)

        dx_next[0]  = dx_next_0
        acc_next[0] = acc_next_0
        vld_next[0] = vld_next_0

        for i in range(1, N + 1):
            prod_trunc   = fp_mul_trunc(self.dx_r[i-1], self.acc_r[i-1], DW, FB)
            coeff_k      = self.coeffs[N - i]
            acc_new      = to_signed(
                to_unsigned(coeff_k + prod_trunc, DW), DW)
            dx_next[i]   = self.dx_r[i-1]
            acc_next[i]  = acc_new
            vld_next[i]  = self.vld_r[i-1]

        # Commit
        self.dx_r  = dx_next
        self.acc_r = acc_next
        self.vld_r = vld_next

        return self.vld_r[N], self.acc_r[N]


# ---------------------------------------------------------------------------
# Helper: coloured terminal output (works on most terminals)
# ---------------------------------------------------------------------------

GREEN  = "\033[92m"
RED    = "\033[91m"
CYAN   = "\033[96m"
YELLOW = "\033[93m"
RESET  = "\033[0m"


def _pass(label, got, expected):
    print(f"  {GREEN}[PASS]{RESET}  {label:20s}  "
          f"y = {got:12.7f}   (expected {expected:12.7f})")


def _fail(label, got, expected, err):
    print(f"  {RED}[FAIL]{RESET}  {label:20s}  "
          f"y = {got:12.7f}   (expected {expected:12.7f})   err = {err:.3e}")


# ---------------------------------------------------------------------------
# Test runner – drives a TaylorPolyPipeline exactly as taylor_poly_tb.v does
# ---------------------------------------------------------------------------

def run_test(dut, x_values, expected_reals, label, tol_lsb):
    """
    Apply x_values one per clock, then drain the pipeline.
    Check each output against expected_reals with tolerance tol_lsb LSBs.
    """
    FB     = dut.FB
    scale  = 1 << FB
    tol    = tol_lsb / scale

    # Convert x values to fixed-point
    xfp = [real_to_fp(x, FB, dut.DW) for x in x_values]

    # Expected value shift-register (mirrors testbench)
    sr = [0.0] * dut.latency

    pass_cnt = fail_cnt = 0
    total_cycles = len(x_values) + dut.latency

    for cycle in range(total_cycles):
        # Drive input
        if cycle < len(x_values):
            vin  = True
            xin  = xfp[cycle]
            exp  = expected_reals[cycle]
        else:
            vin  = False
            xin  = 0
            exp  = 0.0

        # Shift expected SR
        sr = [exp] + sr[:-1]

        # Clock the DUT
        vout, yout = dut.clock(vin, xin)

        # Check output when valid
        if vout:
            got      = fp_to_real(yout, FB)
            expected = sr[-1]
            err      = abs(got - expected)
            lbl      = f"{label}  x={x_values[cycle - dut.latency]:.3f}"
            if err <= tol:
                _pass(lbl, got, expected)
                pass_cnt += 1
            else:
                _fail(lbl, got, expected, err)
                fail_cnt += 1

    return pass_cnt, fail_cnt


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    DW    = 32
    FB    = 16
    SCALE = 1 << FB

    # Q16.16 coefficient literals (matching taylor_poly_tb.v)
    FP_1    =  65536   # 1.0
    FP_HALF =  32768   # 0.5
    FP_1_6  =  10923   # ≈1/6  (0x00002AAB)
    FP_1_24 =   2731   # ≈1/24 (0x00000AAB)
    FP_0    =      0

    print("=" * 65)
    print("  Taylor Polynomial Evaluator – Bit-Accurate Python Simulation")
    print(f"  Fixed-point: Q{DW-FB-1}.{FB}   scale = {SCALE}")
    print("=" * 65)

    total_pass = 0
    total_fail = 0

    # ── Test 1: f(x) = x²   ORDER=2, X0=0 ───────────────────────────────────
    print(f"\n{CYAN}--- DUT 1: f(x) = x^2   (ORDER=2, latency=3 cycles) ---{RESET}")

    # coeffs[k] = c_k  →  c0=0, c1=0, c2=1.0
    quad_coeffs = [FP_0, FP_0, FP_1, FP_0, FP_0, FP_0, FP_0, FP_0]

    dut_quad = TaylorPolyPipeline(
        coeffs=quad_coeffs, order=2, x0=0, data_width=DW, frac_bits=FB)

    x_quad    = [-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0]
    exp_quad  = [x**2 for x in x_quad]

    p, f = run_test(dut_quad, x_quad, exp_quad, "x^2", tol_lsb=2)
    total_pass += p;  total_fail += f

    # ── Test 2: e^x  5-term  ORDER=4, X0=0 ───────────────────────────────────
    print(f"\n{CYAN}--- DUT 2: e^x  5-term approx (ORDER=4, latency=5 cycles) ---{RESET}")

    # coeffs[k] = c_k  →  c0=1, c1=1, c2=1/2, c3=1/6, c4=1/24
    exp_coeffs = [FP_1, FP_1, FP_HALF, FP_1_6, FP_1_24,
                  FP_0, FP_0, FP_0]

    dut_exp = TaylorPolyPipeline(
        coeffs=exp_coeffs, order=4, x0=0, data_width=DW, frac_bits=FB)

    x_exp = [-1.0, 0.0, 0.5, 1.0, 2.0, 0.25]
    # True polynomial value (5-term), NOT exp() — this is what the hardware computes
    exp_exp = [
        sum((xv**k) / _factorial(k) for k in range(5))
        for xv in x_exp
    ]

    p, f = run_test(dut_exp, x_exp, exp_exp, "e^x", tol_lsb=16)
    total_pass += p;  total_fail += f

    # ── Summary ───────────────────────────────────────────────────────────────
    print("\n" + "=" * 65)
    print(f"  Results:  {total_pass} PASSED  /  {total_fail} FAILED")
    if total_fail == 0:
        print(f"  {GREEN}ALL TESTS PASSED{RESET}")
    else:
        print(f"  {RED}*** FAILURES DETECTED ***{RESET}")
    print("=" * 65)

    # ── Extra: show approximation vs exact e^x ────────────────────────────────
    print(f"\n{YELLOW}--- Approximation quality: 5-term e^x vs math.exp() ---{RESET}")
    import math
    print(f"  {'x':>6}  {'5-term':>12}  {'math.exp':>12}  {'abs err':>12}")
    print(f"  {'-'*6}  {'-'*12}  {'-'*12}  {'-'*12}")
    for xv in [-1.0, -0.5, 0.0, 0.5, 1.0, 1.5, 2.0]:
        poly5 = sum((xv**k) / _factorial(k) for k in range(5))
        exact = math.exp(xv)
        print(f"  {xv:6.2f}  {poly5:12.7f}  {exact:12.7f}  {abs(poly5-exact):12.2e}")

    return 0 if total_fail == 0 else 1


def _factorial(n):
    r = 1
    for i in range(2, n + 1):
        r *= i
    return r


if __name__ == "__main__":
    sys.exit(main())
