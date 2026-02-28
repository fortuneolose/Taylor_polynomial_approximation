// =============================================================================
// taylor_poly.v
// Pipelined Taylor-Series Polynomial Evaluator (Horner's Method)
//
// Evaluates the degree-ORDER polynomial
//
//   P(x) = c0 + c1*(x-x0) + c2*(x-x0)^2 + ... + cN*(x-x0)^N
//
// where every c_k = f^(k)(x0) / k! is a pre-computed, fixed-point constant.
//
// Implementation uses Horner's method to avoid explicit power computation:
//
//   P(x) = c0 + dx*(c1 + dx*(c2 + ... + dx*(c_{N-1} + dx*cN)...))
//   dx    = x - x0
//
// Each Horner iteration occupies one pipeline stage, so the overall latency
// is LATENCY = ORDER + 1 clock cycles (one extra cycle to register dx and
// seed the accumulator with the highest-degree coefficient).
//
// ─── Fixed-point format ───────────────────────────────────────────────────────
//   All ports and parameters use signed Q(IW.FRAC_BITS) representation where
//   IW = DATA_WIDTH - FRAC_BITS - 1  (1 sign bit, IW magnitude bits, FRAC_BITS
//   fractional bits).
//
//   For the default DATA_WIDTH=32, FRAC_BITS=16 the representable range is
//   [-32768 , 32767.9999847...] with resolution 2^-16 ≈ 1.526e-5.
//
// ─── Coefficient packing ──────────────────────────────────────────────────────
//   COEFFS is a flat (MAX_ORDER+1)*DATA_WIDTH-bit vector.
//   Coefficient c_k occupies bits  [k*DATA_WIDTH +: DATA_WIDTH].
//   Pack the vector as  {c_{MAX_ORDER}, ..., c_1, c_0}  (c_0 at the LSBs).
//
// ─── Parameters ───────────────────────────────────────────────────────────────
//   DATA_WIDTH  Total bit-width of every data path (default 32).
//   FRAC_BITS   Fractional bits in the Q format     (default 16).
//   ORDER       Polynomial degree, 0 ≤ ORDER ≤ MAX_ORDER.
//   MAX_ORDER   Upper bound on ORDER; sizes the COEFFS array (default 7).
//   X0          Expansion point in fixed-point format.
//   COEFFS      Packed Taylor coefficients (see above).
//
// ─── Corrections from review ──────────────────────────────────────────────────
//   1. Reset changed to synchronous active-low (removed negedge rst_n from
//      sensitivity lists).
//   2. dx_r reduced to [0:ORDER-1]; the dead dx_r[ORDER] register is removed.
//      Stage ORDER reads dx_r[ORDER-1] directly (ORDER=0 special-cased).
//   3. Accumulator addition widened to DATA_WIDTH+1 bits with saturation back
//      to DATA_WIDTH to prevent silent two's-complement overflow.
//   4. Truncation replaced with round-half-up: the full product is biased by
//      0.5 LSB (bit FRAC_BITS-1) before the fractional bits are discarded,
//      removing the systematic negative bias introduced by plain truncation.
//   5. LATENCY localparameter added for self-documenting pipeline depth.
// =============================================================================

`timescale 1ns/1ps

module taylor_poly #(
    parameter                                        DATA_WIDTH = 32,
    parameter                                        FRAC_BITS  = 16,
    parameter                                        ORDER      = 3,
    parameter                                        MAX_ORDER  = 7,
    parameter signed [DATA_WIDTH-1:0]                X0         = {DATA_WIDTH{1'b0}},
    parameter        [(MAX_ORDER+1)*DATA_WIDTH-1:0]  COEFFS     = {(MAX_ORDER+1)*DATA_WIDTH{1'b0}}
)(
    input  wire                         clk,
    input  wire                         rst_n,    // active-low synchronous reset
    input  wire                         valid_in,
    input  wire signed [DATA_WIDTH-1:0] x_in,
    output wire                         valid_out,
    output wire signed [DATA_WIDTH-1:0] y_out
);

    // -------------------------------------------------------------------------
    // Derived constants
    // -------------------------------------------------------------------------

    // Pipeline depth exposed for downstream consumers.
    localparam LATENCY = ORDER + 1;

    // Saturation limits for signed DATA_WIDTH arithmetic.
    localparam signed [DATA_WIDTH-1:0] SAT_MAX =  {1'b0, {(DATA_WIDTH-1){1'b1}}};  //  2^(DATA_WIDTH-1) - 1
    localparam signed [DATA_WIDTH-1:0] SAT_MIN =  {1'b1, {(DATA_WIDTH-1){1'b0}}};  // -2^(DATA_WIDTH-1)

    // Rounding constant: 0.5 LSB in the full-precision Q(2*IW . 2*FRAC_BITS)
    // domain equals 1 << (FRAC_BITS - 1), applied before the fractional
    // truncation step.  Declared as a wire so it can be used in generate scope.
    localparam [2*DATA_WIDTH-1:0] ROUND_HALF =
        (FRAC_BITS > 0) ? ({{(2*DATA_WIDTH-1){1'b0}}, 1'b1} << (FRAC_BITS - 1))
                        : {2*DATA_WIDTH{1'b0}};

    // -------------------------------------------------------------------------
    // Pipeline register arrays
    //   dx_r  : dx = x - x0, forwarded through stages 0 ... ORDER-1.
    //           Only ORDER stages need it; stage ORDER consumes dx_r[ORDER-1]
    //           so the dead dx_r[ORDER] register from the original is removed.
    //           The array is guarded for the ORDER=0 degenerate case (no Horner
    //           iterations, dx is never read) by keeping a minimum size of 1.
    //   acc_r : Horner accumulator, stages 0 ... ORDER.
    //   vld_r : valid token, stages 0 ... ORDER.
    // -------------------------------------------------------------------------
    localparam DX_STAGES = (ORDER > 0) ? ORDER : 1;  // avoid zero-size array

    reg signed [DATA_WIDTH-1:0] dx_r  [0:DX_STAGES-1];
    reg signed [DATA_WIDTH-1:0] acc_r [0:ORDER];
    reg                         vld_r [0:ORDER];

    // ─── Stage 0: register dx and seed accumulator with c_ORDER ──────────────
    // FIX 1: synchronous reset — rst_n removed from sensitivity list.
    always @(posedge clk) begin
        if (!rst_n) begin
            dx_r[0]  <= {DATA_WIDTH{1'b0}};
            acc_r[0] <= {DATA_WIDTH{1'b0}};
            vld_r[0] <= 1'b0;
        end else begin
            dx_r[0]  <= x_in - X0;
            acc_r[0] <= COEFFS[ORDER*DATA_WIDTH +: DATA_WIDTH];
            vld_r[0] <= valid_in;
        end
    end

    // ─── Stages 1 ... ORDER: Horner iterations ────────────────────────────────
    // At stage i:
    //   acc = saturate( c_{ORDER-i} + round( dx * acc_prev ) )
    //
    // FIX 4 (rounding): a 0.5-LSB bias is added to the full 2*DATA_WIDTH-bit
    //   product before discarding the lower FRAC_BITS bits, replacing the
    //   original truncation that had a systematic negative bias.
    //
    // FIX 3 (saturation): the addition is performed at DATA_WIDTH+1 bits and
    //   clamped to [SAT_MIN, SAT_MAX] before registering, preventing the silent
    //   two's-complement wrap that the original code was susceptible to.
    //
    // FIX 2 (dead register): dx_r is forwarded only up to index ORDER-1.
    //   The last iteration (i = ORDER) reads dx_r[ORDER-1] without writing a
    //   new dx_r entry.
    genvar i;
    generate
        for (i = 1; i <= ORDER; i = i + 1) begin : g_horner
            // ── Combinational datapath ────────────────────────────────────────
            // Full-precision signed product: Q(2*IW . 2*FRAC_BITS)
            wire signed [2*DATA_WIDTH-1:0]  prod_full;
            // Biased product for round-half-up
            wire        [2*DATA_WIDTH-1:0]  prod_biased;
            // Product rounded and truncated back to Q(IW.FRAC_BITS)
            wire signed [DATA_WIDTH-1:0]    prod_round;
            // Taylor coefficient for this stage: c_{ORDER - i}
            wire signed [DATA_WIDTH-1:0]    coeff_k;
            // Widened sum before saturation (DATA_WIDTH+1 bits)
            wire signed [DATA_WIDTH:0]      sum_wide;

            // Coefficient and multiplication
            assign coeff_k     = COEFFS[(ORDER-i)*DATA_WIDTH +: DATA_WIDTH];
            assign prod_full   = $signed(dx_r[i-1]) * $signed(acc_r[i-1]);

            // FIX 4: round-half-up — bias by 2^(FRAC_BITS-1) then truncate.
            assign prod_biased = prod_full + ROUND_HALF;
            assign prod_round  = prod_biased[DATA_WIDTH+FRAC_BITS-1 : FRAC_BITS];

            // FIX 3: widen addition to detect overflow, then saturate.
            assign sum_wide    = {coeff_k[DATA_WIDTH-1], coeff_k}
                               + {prod_round[DATA_WIDTH-1], prod_round};

            // ── Pipeline registers ────────────────────────────────────────────
            // FIX 1: synchronous reset.
            // FIX 2: only forward dx when there is a next stage to consume it.
            always @(posedge clk) begin
                if (!rst_n) begin
                    acc_r[i] <= {DATA_WIDTH{1'b0}};
                    vld_r[i] <= 1'b0;
                    if (i < ORDER)
                        dx_r[i] <= {DATA_WIDTH{1'b0}};
                end else begin
                    // Saturate accumulator
                    if (sum_wide[DATA_WIDTH] != sum_wide[DATA_WIDTH-1])
                        // Sign bits differ → overflow occurred
                        acc_r[i] <= sum_wide[DATA_WIDTH] ? SAT_MIN : SAT_MAX;
                    else
                        acc_r[i] <= sum_wide[DATA_WIDTH-1:0];

                    vld_r[i] <= vld_r[i-1];

                    // FIX 2: suppress dx forwarding at the final stage
                    if (i < ORDER)
                        dx_r[i] <= dx_r[i-1];
                end
            end
        end
    endgenerate

    // ─── Output ───────────────────────────────────────────────────────────────
    assign y_out     = acc_r[ORDER];
    assign valid_out = vld_r[ORDER];

endmodule
