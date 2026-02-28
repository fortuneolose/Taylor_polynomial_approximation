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
// is ORDER + 1 clock cycles (one extra cycle to register dx and seed the
// accumulator with the highest-degree coefficient).
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
    // Pipeline register arrays  (index 0 = stage 0, index ORDER = output stage)
    // -------------------------------------------------------------------------
    reg signed [DATA_WIDTH-1:0] dx_r  [0:ORDER];  // dx = x - x0 forwarded each stage
    reg signed [DATA_WIDTH-1:0] acc_r [0:ORDER];  // Horner accumulator
    reg                         vld_r [0:ORDER];  // valid token

    // ─── Stage 0: register dx and seed accumulator with c_ORDER ──────────────
    always @(posedge clk or negedge rst_n) begin
        if (!rst_n) begin
            dx_r[0]  <= {DATA_WIDTH{1'b0}};
            acc_r[0] <= {DATA_WIDTH{1'b0}};
            vld_r[0] <= 1'b0;
        end else begin
            dx_r[0]  <= x_in - X0;
            // Seed with the highest-degree Taylor coefficient c_ORDER.
            // The part-select is a constant expression at elaboration time.
            acc_r[0] <= COEFFS[ORDER*DATA_WIDTH +: DATA_WIDTH];
            vld_r[0] <= valid_in;
        end
    end

    // ─── Stages 1 … ORDER: Horner iterations ─────────────────────────────────
    // At stage i:
    //   acc = c_{ORDER-i}  +  dx * acc_prev
    //
    // Multiplication produces a 2*DATA_WIDTH-bit full-precision product in
    // Q(2*IW . 2*FRAC_BITS) format.  Truncating back to Q(IW.FRAC_BITS) is
    // achieved by discarding the lower FRAC_BITS bits:
    //   product_trunc = product[DATA_WIDTH+FRAC_BITS-1 : FRAC_BITS]
    // This retains the sign bit and DATA_WIDTH-1 lower magnitude/fraction bits.
    genvar i;
    generate
        for (i = 1; i <= ORDER; i = i + 1) begin : g_horner
            // Full-width signed product (no overflow loss at this width)
            wire signed [2*DATA_WIDTH-1:0] prod_full;
            // Product rounded back to the original Q format
            wire signed [DATA_WIDTH-1:0]   prod_trunc;
            // Taylor coefficient for this stage: c_{ORDER - i}
            wire signed [DATA_WIDTH-1:0]   coeff_k;

            assign coeff_k    = COEFFS[(ORDER-i)*DATA_WIDTH +: DATA_WIDTH];
            assign prod_full  = dx_r[i-1] * acc_r[i-1];
            assign prod_trunc = prod_full[DATA_WIDTH+FRAC_BITS-1 : FRAC_BITS];

            always @(posedge clk or negedge rst_n) begin
                if (!rst_n) begin
                    dx_r[i]  <= {DATA_WIDTH{1'b0}};
                    acc_r[i] <= {DATA_WIDTH{1'b0}};
                    vld_r[i] <= 1'b0;
                end else begin
                    dx_r[i]  <= dx_r[i-1];
                    acc_r[i] <= coeff_k + prod_trunc;
                    vld_r[i] <= vld_r[i-1];
                end
            end
        end
    endgenerate

    // ─── Output ───────────────────────────────────────────────────────────────
    assign y_out     = acc_r[ORDER];
    assign valid_out = vld_r[ORDER];

endmodule
