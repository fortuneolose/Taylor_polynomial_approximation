// =============================================================================
// taylor_poly_tb.v
// Testbench for taylor_poly.v   (Verilog-2001 compatible)
//
// Two DUT instances are exercised:
//
//  dut_quad  -  Exact quadratic  f(x) = x^2  (ORDER=2, around x0=0)
//               c0=0, c1=0, c2=1.0 -> P(x)=x^2, no approximation error.
//
//  dut_exp   -  e^x, 5-term series (ORDER=4, around x0=0)
//               c_k = 1/k!  for k=0..4
//
// Both DUTs share clock and reset.
// Results are verified using a real-valued expected-value shift register that
// accounts for the pipeline latency (ORDER + 1 cycles per DUT).
//
// Fixed-point format: signed Q(15.16)  ->  scale = 65536.
//
// --- Corrections from review --------------------------------------------------
//   FIX A  real_to_fp used $rtoi() which truncates toward zero, introducing
//          a systematic +-1 LSB error in every stimulus value.  Replaced with
//          round-half-away-from-zero: add +0.5 for positive, -0.5 for negative
//          before calling $rtoi(), consistent with the DUT's round-half-up
//          policy.
//
//   FIX B  The e^x expected values (ee[]) were computed from exact mathematical
//          coefficients (1/6, 1/24, ...) rather than the actual fixed-point
//          coefficient values loaded into the DUT (FP_1_6 ~= 10923/65536,
//          FP_1_24 ~= 2731/65536).  For x=2 this created a systematic ~24 LSB
//          gap between ee[] and the DUT output -- guaranteed false failure.
//          Fixed by deriving ee[] from fp_to_real(EXP_COEFF_k) at runtime,
//          so expected values match what the DUT was actually programmed with.
//
//   FIX C  The check_exp tolerance comment was misleading: after FIX B the
//          only residual error is round-half-up rounding across 4 multiply
//          stages.  Worst case (dx=2, ORDER=4) is 0.5*(1+2+4+8)=7.5 LSBs, so
//          16 LSBs remains the correct tolerance -- but the comment now
//          documents this derivation explicitly rather than saying "a few LSBs".
// =============================================================================

`timescale 1ns/1ps

module taylor_poly_tb;

    // --- Parameters & constants -----------------------------------------------
    localparam CLK_PERIOD = 10;    // ns
    localparam DW         = 32;    // DATA_WIDTH
    localparam FB         = 16;    // FRAC_BITS
    localparam SCALE      = 65536; // 1 << FB

    // Q16.16 coefficient literals
    localparam signed [DW-1:0] FP_1    = 32'h00010000; // 1.0
    localparam signed [DW-1:0] FP_HALF = 32'h00008000; // 0.5
    localparam signed [DW-1:0] FP_1_6  = 32'h00002AAB; // ~1/6   10923/65536
    localparam signed [DW-1:0] FP_1_24 = 32'h00000AAB; // ~1/24   2731/65536
    localparam signed [DW-1:0] FP_0    = 32'h00000000;

    // --- Clock / reset --------------------------------------------------------
    reg clk;
    reg rst_n;
    initial clk = 0;
    always #(CLK_PERIOD/2) clk = ~clk;

    // =========================================================================
    // DUT 1: f(x) = x^2   ORDER=2, MAX_ORDER=7
    // =========================================================================
    localparam QUAD_ORD = 2;
    localparam QUAD_MAX = 7;
    localparam QUAD_LAT = QUAD_ORD + 1;

    localparam [(QUAD_MAX+1)*DW-1:0] QUAD_COEFFS =
        { FP_0,   // c7
          FP_0,   // c6
          FP_0,   // c5
          FP_0,   // c4
          FP_0,   // c3
          FP_1,   // c2 = 1.0
          FP_0,   // c1 = 0
          FP_0 }; // c0 = 0

    reg                      quad_vin;
    reg  signed [DW-1:0]     quad_xin;
    wire                     quad_vout;
    wire signed [DW-1:0]     quad_yout;

    taylor_poly #(
        .DATA_WIDTH ( DW          ),
        .FRAC_BITS  ( FB          ),
        .ORDER      ( QUAD_ORD    ),
        .MAX_ORDER  ( QUAD_MAX    ),
        .X0         ( FP_0        ),
        .COEFFS     ( QUAD_COEFFS )
    ) dut_quad (
        .clk       ( clk       ),
        .rst_n     ( rst_n     ),
        .valid_in  ( quad_vin  ),
        .x_in      ( quad_xin  ),
        .valid_out ( quad_vout ),
        .y_out     ( quad_yout )
    );

    // =========================================================================
    // DUT 2: e^x  5-term  ORDER=4, MAX_ORDER=7
    // =========================================================================
    localparam EXP_ORD = 4;
    localparam EXP_MAX = 7;
    localparam EXP_LAT = EXP_ORD + 1;

    localparam [(EXP_MAX+1)*DW-1:0] EXP_COEFFS =
        { FP_0,      // c7
          FP_0,      // c6
          FP_0,      // c5
          FP_1_24,   // c4
          FP_1_6,    // c3
          FP_HALF,   // c2
          FP_1,      // c1
          FP_1 };    // c0

    reg                      exp_vin;
    reg  signed [DW-1:0]     exp_xin;
    wire                     exp_vout;
    wire signed [DW-1:0]     exp_yout;

    taylor_poly #(
        .DATA_WIDTH ( DW         ),
        .FRAC_BITS  ( FB         ),
        .ORDER      ( EXP_ORD    ),
        .MAX_ORDER  ( EXP_MAX    ),
        .X0         ( FP_0       ),
        .COEFFS     ( EXP_COEFFS )
    ) dut_exp (
        .clk       ( clk       ),
        .rst_n     ( rst_n     ),
        .valid_in  ( exp_vin   ),
        .x_in      ( exp_xin   ),
        .valid_out ( exp_vout  ),
        .y_out     ( exp_yout  )
    );

    // --- Fixed-point conversion helpers ---------------------------------------

    // Convert signed Q16.16 to real.
    function real fp_to_real;
        input signed [DW-1:0] fp;
        begin
            fp_to_real = $itor(fp) / $itor(SCALE);
        end
    endfunction

    // FIX A: Convert real to nearest signed Q16.16 using round-half-away-from-
    // zero.  The original $rtoi() truncated toward zero, introducing a
    // systematic +-1 LSB stimulus error.  Adding +-0.5 before truncation gives
    // correct nearest-integer rounding, matching the DUT's round-half-up
    // policy for positive products.
    function signed [DW-1:0] real_to_fp;
        input real r;
        begin
            if (r >= 0.0)
                real_to_fp = $rtoi(r * $itor(SCALE) + 0.5);
            else
                real_to_fp = $rtoi(r * $itor(SCALE) - 0.5);
        end
    endfunction

    // --- FIX B: real-valued DUT coefficient values ----------------------------
    // Derived at runtime from the same localparam literals used to configure
    // the DUT, so expected values track the actual fixed-point quantisation
    // rather than the exact mathematical values.
    real ec0, ec1, ec2, ec3, ec4; // e^x DUT coefficients as reals

    // --- Stimulus and expected-value tables -----------------------------------
    localparam QUAD_N = 7;
    reg  signed [DW-1:0] qx [0:QUAD_N-1];
    real                 qe [0:QUAD_N-1];  // expected: x^2

    localparam EXP_N = 6;
    reg  signed [DW-1:0] ex [0:EXP_N-1];
    real                 ee [0:EXP_N-1];  // expected: 5-term poly (FP coeffs)

    // --- Expected-value shift registers ---------------------------------------
    real quad_sr [0:QUAD_LAT-1];
    real exp_sr  [0:EXP_LAT-1];

    // --- Pass / fail counters -------------------------------------------------
    integer pass_cnt, fail_cnt;

    // --- Loop variables -------------------------------------------------------
    integer qi, ei, j;

    // --- Checker tasks --------------------------------------------------------
    task check_quad;
        input real expected;
        real got, err;
        begin
            got = fp_to_real(quad_yout);
            err = (got > expected) ? (got - expected) : (expected - got);
            // Tolerance: 2 LSBs covers the round-half-up error on both
            // pipeline stages for ORDER=2.
            if (err <= (2.0 / $itor(SCALE))) begin
                $display("  [QUAD PASS]  x^2 :  y = %10.6f  (expected %10.6f)",
                         got, expected);
                pass_cnt = pass_cnt + 1;
            end else begin
                $display("  [QUAD FAIL]  x^2 :  y = %10.6f  (expected %10.6f)  err=%e",
                         got, expected, err);
                fail_cnt = fail_cnt + 1;
            end
        end
    endtask

    task check_exp;
        input real expected;
        real got, err;
        begin
            got = fp_to_real(exp_yout);
            err = (got > expected) ? (got - expected) : (expected - got);
            // FIX B+C: now that expected is computed from the actual FP
            // coefficients, the only residual error is round-half-up rounding
            // across ORDER=4 multiply stages.  Worst case (dx=2):
            //   0.5 * (1 + dx + dx^2 + dx^3) = 0.5 * (1+2+4+8) = 7.5 LSBs.
            // 16 LSBs gives a comfortable margin for all test inputs.
            if (err <= (16.0 / $itor(SCALE))) begin
                $display("  [EXP  PASS]  e^x :  y = %10.6f  (expected %10.6f)",
                         got, expected);
                pass_cnt = pass_cnt + 1;
            end else begin
                $display("  [EXP  FAIL]  e^x :  y = %10.6f  (expected %10.6f)  err=%e",
                         got, expected, err);
                fail_cnt = fail_cnt + 1;
            end
        end
    endtask

    // --- Stimulus -------------------------------------------------------------
    initial begin
        // -- FIX B: resolve DUT e^x coefficients to real once at startup -------
        // Evaluating fp_to_real() on the same localparam literals that were
        // packed into EXP_COEFFS means ee[] will match exactly what the DUT
        // was programmed with, irrespective of quantisation.
        ec0 = fp_to_real(FP_1);
        ec1 = fp_to_real(FP_1);
        ec2 = fp_to_real(FP_HALF);
        ec3 = fp_to_real(FP_1_6);
        ec4 = fp_to_real(FP_1_24);

        // -- Populate stimulus tables ------------------------------------------

        // Quadratic: expected = x^2  (exact, no coefficient approximation)
        qx[0] = real_to_fp(-3.0); qe[0] =  9.0;
        qx[1] = real_to_fp(-2.0); qe[1] =  4.0;
        qx[2] = real_to_fp(-1.0); qe[2] =  1.0;
        qx[3] = real_to_fp( 0.0); qe[3] =  0.0;
        qx[4] = real_to_fp( 1.0); qe[4] =  1.0;
        qx[5] = real_to_fp( 2.0); qe[5] =  4.0;
        qx[6] = real_to_fp( 3.0); qe[6] =  9.0;

        // FIX B: e^x expected = P(x) using the actual FP coefficient values,
        // not the mathematically exact 1/k! values.  For x=2, the original
        // code's use of exact 1/6 and 1/24 created a ~24 LSB gap (c3 error
        // ~13.6 LSBs, c4 error ~10.7 LSBs) that always tripped check_exp.
        begin : exp_stim
            real xv;

            xv = fp_to_real(real_to_fp(-1.0));
            ex[0] = real_to_fp(xv);
            ee[0] = ec0 + ec1*xv + ec2*xv*xv + ec3*xv*xv*xv + ec4*xv*xv*xv*xv;

            xv = 0.0;
            ex[1] = real_to_fp(xv);
            ee[1] = ec0;

            xv = fp_to_real(real_to_fp(0.5));
            ex[2] = real_to_fp(xv);
            ee[2] = ec0 + ec1*xv + ec2*xv*xv + ec3*xv*xv*xv + ec4*xv*xv*xv*xv;

            xv = fp_to_real(real_to_fp(1.0));
            ex[3] = real_to_fp(xv);
            ee[3] = ec0 + ec1*xv + ec2*xv*xv + ec3*xv*xv*xv + ec4*xv*xv*xv*xv;

            xv = fp_to_real(real_to_fp(2.0));
            ex[4] = real_to_fp(xv);
            ee[4] = ec0 + ec1*xv + ec2*xv*xv + ec3*xv*xv*xv + ec4*xv*xv*xv*xv;

            xv = fp_to_real(real_to_fp(0.25));
            ex[5] = real_to_fp(xv);
            ee[5] = ec0 + ec1*xv + ec2*xv*xv + ec3*xv*xv*xv + ec4*xv*xv*xv*xv;
        end

        // -- VCD waveform dump -------------------------------------------------
        $dumpfile("taylor_poly.vcd");
        $dumpvars(0, taylor_poly_tb);

        // -- Initialise shift registers and counters ---------------------------
        for (j = 0; j < QUAD_LAT; j = j + 1) quad_sr[j] = 0.0;
        for (j = 0; j < EXP_LAT;  j = j + 1) exp_sr[j]  = 0.0;
        pass_cnt = 0;
        fail_cnt = 0;

        // -- Reset -------------------------------------------------------------
        rst_n    = 0;
        quad_vin = 0;  quad_xin = 0;
        exp_vin  = 0;  exp_xin  = 0;
        repeat (4) @(posedge clk);
        @(negedge clk); rst_n = 1;

        // ----------------------------------------------------------------------
        $display("=============================================================");
        $display("  Taylor Polynomial Evaluator - Simulation");
        $display("  Format: signed Q%0d.%0d   scale = %0d", DW-FB-1, FB, SCALE);
        $display("  e^x FP coefficients:");
        $display("    c0=%.8f  c1=%.8f  c2=%.8f", ec0, ec1, ec2);
        $display("    c3=%.8f  c4=%.8f", ec3, ec4);
        $display("=============================================================");

        // -- Run DUT 1 (quadratic) ---------------------------------------------
        $display("\n--- DUT 1: f(x) = x^2  (ORDER=%0d, latency=%0d cycles) ---",
                 QUAD_ORD, QUAD_LAT);

        exp_vin = 0; exp_xin = 0;

        for (qi = 0; qi < QUAD_N + QUAD_LAT; qi = qi + 1) begin
            @(negedge clk);

            if (qi < QUAD_N) begin
                quad_xin = qx[qi];
                quad_vin = 1;
            end else begin
                quad_vin = 0;
            end

            for (j = QUAD_LAT-1; j > 0; j = j - 1)
                quad_sr[j] = quad_sr[j-1];
            quad_sr[0] = (qi < QUAD_N) ? qe[qi] : 0.0;

            @(posedge clk); #1;
            if (quad_vout)
                check_quad(quad_sr[QUAD_LAT-1]);
        end

        quad_vin = 0;

        // -- Run DUT 2 (e^x) ---------------------------------------------------
        $display("\n--- DUT 2: e^x (5-term, ORDER=%0d, latency=%0d cycles) ---",
                 EXP_ORD, EXP_LAT);

        for (ei = 0; ei < EXP_N + EXP_LAT; ei = ei + 1) begin
            @(negedge clk);

            if (ei < EXP_N) begin
                exp_xin = ex[ei];
                exp_vin = 1;
            end else begin
                exp_vin = 0;
            end

            for (j = EXP_LAT-1; j > 0; j = j - 1)
                exp_sr[j] = exp_sr[j-1];
            exp_sr[0] = (ei < EXP_N) ? ee[ei] : 0.0;

            @(posedge clk); #1;
            if (exp_vout)
                check_exp(exp_sr[EXP_LAT-1]);
        end

        // -- Summary -----------------------------------------------------------
        $display("\n=============================================================");
        $display("  Results:  %0d PASSED  /  %0d FAILED", pass_cnt, fail_cnt);
        if (fail_cnt == 0)
            $display("  ALL TESTS PASSED");
        else
            $display("  *** FAILURES DETECTED ***");
        $display("=============================================================");

        $finish;
    end

    // --- Timeout watchdog -----------------------------------------------------
    initial begin
        #100000;
        $display("TIMEOUT");
        $finish;
    end

endmodule
