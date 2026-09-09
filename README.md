# Optimization Techniques, Assignment 1

Two solvers in plain Python 3. No external libraries. Every quantity is a
`fractions.Fraction`, so the results are exact rather than rounded.

    python3 assignment/big_m_simplex.py
    python3 assignment/transportation_vam_modi.py
    python3 assignment/transportation_vam_modi.py --vam-only

Saved runs of all three commands are in `outputs/`.

## 1. Big-M simplex

The case study is the diet, or feed-blending, problem. A dairy blends two
supplements into 1 kg batches of cattle feed. Supplement A costs Rs 4/kg and
B costs Rs 1/kg. The blend has to hit a protein specification exactly, clear
a vitamin floor, and fit inside the mixing drum.

Decision variables:

    x1 = kg of supplement A per batch
    x2 = kg of supplement B per batch

The LPP:

    minimise    Z = 4*x1 + x2

    subject to  3*x1 +  x2  = 3      protein, exact specification
                4*x1 + 3*x2 >= 6     vitamin premix floor
                 x1 + 2*x2  <= 4     mixing drum capacity
                x1, x2 >= 0

The program converts this to standard form itself. An equality gets one
artificial variable, a `>=` row gets a surplus and an artificial, and a `<=`
row gets a slack:

| constraint | variables added |
|---|---|
| `3*x1 + x2 = 3` | artificial `a1` |
| `4*x1 + 3*x2 >= 6` | surplus `e2`, artificial `a2` |
| `x1 + 2*x2 <= 4` | slack `s3` |

Each row then owns exactly one basic variable, which gives the simplex a
starting basis. The origin is not feasible here, which is why the artificials
are there.

The objective penalises each artificial by `+M`. Rather than substituting
some large number for M, the program keeps the objective row as two parallel
Fraction rows, one holding the M coefficients and one holding the constants.
Comparing `(M coefficient, constant)` as a pair reproduces the rule that any
M term outweighs any constant. That avoids floating-point error, and it avoids
having to guess a numeric M large enough for the data.

Result, after 3 simplex iterations:

    x1 = 2/5 = 0.4 kg
    x2 = 9/5 = 1.8 kg
    Z  = 17/5 = Rs 3.40 per batch

## 2. Transportation problem (VAM and MODI)

The second case study is cement distribution. Three plants ship to four regional
warehouses. Output and requirement are fixed for the month, and freight per
tonne depends on which plant supplies which warehouse.

    capacity (tonnes)  S1 = 20   S2 = 30   S3 = 50            total 100
    demand   (tonnes)  D1 = 10   D2 = 40   D3 = 20   D4 = 30  total 100

    freight c_ij, Rs per tonne

            D1    D2    D3    D4
      S1     8     6    10     9
      S2     9    12    13     7
      S3    14     9    16     5

With `x_ij` as the tonnes shipped from plant i to warehouse j, the LPP is:

    minimise    Z = sum_i sum_j  c_ij * x_ij

    subject to  sum_j x_ij = a_i    each plant ships its full output
                sum_i x_ij = b_j    each warehouse gets its requirement
                x_ij >= 0

Supply and demand balance in this case study. When they do not, the program
inserts a zero-cost dummy row or column before solving.

### Phase 1, Vogel's Approximation Method

VAM builds the initial basic feasible solution. It finds the row or column
with the largest penalty, meaning the gap between that line's two cheapest
cells, then fills the cheapest cell in that line as far as the remaining
supply and demand allow. It repeats until everything is allocated. Pass
`--vam-only` to stop here.

Initial BFS cost: Rs 800.

### Phase 2, MODI

MODI tests that solution for optimality. It solves `u_i + v_j = c_ij` across
the basic cells, anchored at `u_0 = 0`, to get shadow prices. Any non-basic
cell whose opportunity cost `c_ij - u_i - v_j` is negative would lower the
total, so that cell enters the basis and the allocation shifts around a
closed loop until no negative opportunity cost is left.

Optimal shipment plan:

            D1    D2    D3    D4
      S1     -    20     -     -
      S2    10     -    20     -
      S3     -    20     -    30

    minimum total cost = Rs 800

VAM happened to reach the optimum directly here, so MODI spends one
degenerate pivot confirming it rather than reducing the cost. The program
asserts the row sums, the column sums, and the m+n-1 basis count at the end
of every run.
