// =============================================================================
// taylor_poly_tb.v
// Testbench for taylor_poly.v   (Verilog-2001 compatible)
//
// Two DUT instances are exercised:
//
//  dut_quad  –  Exact quadratic  f(x) = x²  (ORDER=2, around x0=0)
//               c0=0, c1=0, c2=1.0 → P(x)=x², no approximation error.
//
//  dut_exp   –  e^x, 5-term series (ORDER=4, around x0=0)
//               c_k = 1/k!  for k=0..4
//
// Both DUTs share clock and reset.
// Results are verified using a real-valued expected-value shift register that
// accounts for the pipeline latency (ORDER + 1 cycles per DUT).
//
// Fixed-point format: signed Q(15.16)  →  scale = 65536.
// =============================================================================

`timescale 1ns/1ps

module taylor_poly_tb;

    // ─── Parameters & constants ───────────────────────────────────────────────
    localparam CLK_PERIOD = 10;   // ns
    localparam DW         = 32;   // DATA_WIDTH
    localparam FB         = 16;   // FRAC_BITS
    localparam SCALE      = 65536; // 1 << FB

    // Q16.16 coefficient literals
    localparam signed [DW-1:0] FP_1    = 32'h00010000; // 1.0
    localparam signed [DW-1:0] FP_HALF = 32'h00008000; // 0.5
    localparam signed [DW-1:0] FP_1_6  = 32'h00002AAB; // ≈1/6  (error<2^-16)
    localparam signed [DW-1:0] FP_1_24 = 32'h00000AAB; // ≈1/24 (error<2^-16)
    localparam signed [DW-1:0] FP_0    = 32'h00000000;

    // ─── Clock / reset ────────────────────────────────────────────────────────
    reg clk;
    reg rst_n;
    initial clk = 0;
    always #(CLK_PERIOD/2) clk = ~clk;

    // =========================================================================
    // DUT 1: f(x) = x²   ORDER=2, MAX_ORDER=7
    //   COEFFS: c0=0, c1=0, c2=1.0, c3..c7=0
    //   Packing: {c7, c6, c5, c4, c3, c2, c1, c0} → c0 at bits [31:0]
    // =========================================================================
    localparam QUAD_ORD = 2;
    localparam QUAD_MAX = 7;
    localparam QUAD_LAT = QUAD_ORD + 1; // pipeline latency in cycles

    localparam [(QUAD_MAX+1)*DW-1:0] QUAD_COEFFS =
        { FP_0,   // c7
          FP_0,   // c6
          FP_0,   // c5
          FP_0,   // c4
          FP_0,   // c3
          FP_1,   // c2 = 1.0
          FP_0,   // c1 = 0
          FP_0 }; // c0 = 0  (LSBs)

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
    //   c0=1, c1=1, c2=1/2, c3=1/6, c4=1/24
    // =========================================================================
    localparam EXP_ORD = 4;
    localparam EXP_MAX = 7;
    localparam EXP_LAT = EXP_ORD + 1;

    localparam [(EXP_MAX+1)*DW-1:0] EXP_COEFFS =
        { FP_0,      // c7  bits[255:224]
          FP_0,      // c6  bits[223:192]
          FP_0,      // c5  bits[191:160]
          FP_1_24,   // c4  bits[159:128]
          FP_1_6,    // c3  bits[127:96]
          FP_HALF,   // c2  bits[95:64]
          FP_1,      // c1  bits[63:32]
          FP_1 };    // c0  bits[31:0]

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

    // ─── Fixed-point conversion helpers ──────────────────────────────────────
    // Convert signed Q16.16 to real
    function real fp_to_real;
        input signed [DW-1:0] fp;
        begin
            fp_to_real = $itor(fp) / $itor(SCALE);
        end
    endfunction

    // Convert real to nearest signed Q16.16
    function signed [DW-1:0] real_to_fp;
        input real r;
        begin
            real_to_fp = $rtoi(r * $itor(SCALE));
        end
    endfunction

    // ─── Module-level stimulus and expected-value tables ─────────────────────
    // ── Quadratic test vectors (7 samples) ───────────────────────────────────
    localparam QUAD_N = 7;
    reg signed [DW-1:0] qx [0:QUAD_N-1];
    real                qe [0:QUAD_N-1];  // expected: x²

    // ── e^x test vectors (6 samples) ─────────────────────────────────────────
    localparam EXP_N = 6;
    reg signed [DW-1:0] ex [0:EXP_N-1];
    real                ee [0:EXP_N-1];  // expected: 5-term polynomial value

    // ─── Expected-value shift registers ──────────────────────────────────────
    real quad_sr [0:QUAD_LAT-1]; // depth = latency
    real exp_sr  [0:EXP_LAT-1];

    // ─── Pass / fail counters ─────────────────────────────────────────────────
    integer pass_cnt, fail_cnt;

    // ─── Loop variables (module scope for Verilog-2001 compatibility) ─────────
    integer qi, ei, j;

    // ─── Checker tasks ────────────────────────────────────────────────────────
    task check_quad;
        input real expected;
        real got, err;
        begin
            got = fp_to_real(quad_yout);
            err = (got > expected) ? (got - expected) : (expected - got);
            // Tolerate 1 LSB of rounding (2^-16)
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
            // Accumulated rounding over 5 stages: allow a few LSBs
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

    // ─── Stimulus ─────────────────────────────────────────────────────────────
    initial begin
        // ── Populate stimulus tables ──────────────────────────────────────────
        // Quadratic: expected = x²
        qx[0] = real_to_fp(-3.0); qe[0] =  9.0;
        qx[1] = real_to_fp(-2.0); qe[1] =  4.0;
        qx[2] = real_to_fp(-1.0); qe[2] =  1.0;
        qx[3] = real_to_fp( 0.0); qe[3] =  0.0;
        qx[4] = real_to_fp( 1.0); qe[4] =  1.0;
        qx[5] = real_to_fp( 2.0); qe[5] =  4.0;
        qx[6] = real_to_fp( 3.0); qe[6] =  9.0;

        // e^x 5-term: expected = 1 + x + x²/2 + x³/6 + x⁴/24
        ex[0] = real_to_fp(-1.0);
        ee[0] = 1.0 + (-1.0) + (1.0/2.0) + (-1.0/6.0) + (1.0/24.0); // ≈0.375

        ex[1] = real_to_fp(0.0);
        ee[1] = 1.0;

        ex[2] = real_to_fp(0.5);
        ee[2] = 1.0 + 0.5 + 0.125 + (0.125/6.0) + (0.0625/24.0);    // ≈1.6484

        ex[3] = real_to_fp(1.0);
        ee[3] = 1.0 + 1.0 + 0.5 + (1.0/6.0) + (1.0/24.0);            // ≈2.7083

        ex[4] = real_to_fp(2.0);
        ee[4] = 1.0 + 2.0 + 2.0 + (8.0/6.0) + (16.0/24.0);           // ≈7.0

        ex[5] = real_to_fp(0.25);
        ee[5] = 1.0 + 0.25 + 0.03125 + (0.015625/6.0) + (0.00390625/24.0);

        // ── VCD waveform dump ─────────────────────────────────────────────────
        $dumpfile("taylor_poly.vcd");
        $dumpvars(0, taylor_poly_tb);

        // ── Initialise shift registers and counters ───────────────────────────
        for (j = 0; j < QUAD_LAT; j = j + 1) quad_sr[j] = 0.0;
        for (j = 0; j < EXP_LAT;  j = j + 1) exp_sr[j]  = 0.0;
        pass_cnt = 0;
        fail_cnt = 0;

        // ── Reset ─────────────────────────────────────────────────────────────
        rst_n    = 0;
        quad_vin = 0;  quad_xin = 0;
        exp_vin  = 0;  exp_xin  = 0;
        repeat (4) @(posedge clk);
        @(negedge clk); rst_n = 1;

        // ─────────────────────────────────────────────────────────────────────
        $display("=============================================================");
        $display("  Taylor Polynomial Evaluator – Simulation");
        $display("  Format: signed Q%0d.%0d   scale = %0d", DW-FB-1, FB, SCALE);
        $display("=============================================================");

        // ── Run DUT 1 (quadratic) ─────────────────────────────────────────────
        $display("\n--- DUT 1: f(x) = x^2  (ORDER=%0d, latency=%0d cycles) ---",
                 QUAD_ORD, QUAD_LAT);

        exp_vin = 0; exp_xin = 0;

        for (qi = 0; qi < QUAD_N + QUAD_LAT; qi = qi + 1) begin
            @(negedge clk);

            // Drive input while samples remain
            if (qi < QUAD_N) begin
                quad_xin = qx[qi];
                quad_vin = 1;
            end else begin
                quad_vin = 0;
            end

            // Shift expected-value pipeline (oldest at index QUAD_LAT-1)
            for (j = QUAD_LAT-1; j > 0; j = j - 1)
                quad_sr[j] = quad_sr[j-1];
            quad_sr[0] = (qi < QUAD_N) ? qe[qi] : 0.0;

            // Check output one delta after posedge (register values settled)
            @(posedge clk); #1;
            if (quad_vout)
                check_quad(quad_sr[QUAD_LAT-1]);
        end

        quad_vin = 0;

        // ── Run DUT 2 (e^x) ──────────────────────────────────────────────────
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

        // ── Summary ───────────────────────────────────────────────────────────
        $display("\n=============================================================");
        $display("  Results:  %0d PASSED  /  %0d FAILED", pass_cnt, fail_cnt);
        if (fail_cnt == 0)
            $display("  ALL TESTS PASSED");
        else
            $display("  *** FAILURES DETECTED ***");
        $display("=============================================================");

        $finish;
    end

    // ─── Timeout watchdog ─────────────────────────────────────────────────────
    initial begin
        #100000;
        $display("TIMEOUT");
        $finish;
    end

endmodule
