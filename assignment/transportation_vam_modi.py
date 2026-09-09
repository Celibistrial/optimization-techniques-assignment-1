"""transportation problem. VAM gets the first BFS then MODI improves it"""

import sys
from fractions import Fraction


class Transportation:
    def __init__(self, supply, demand, cost):
        self.supply = [Fraction(v) for v in supply]
        self.demand = [Fraction(v) for v in demand]
        self.cost = [[Fraction(v) for v in row] for row in cost]
        self._balance()
        self.m = len(self.cost)
        self.n = len(self.cost[0])

    def _balance(self):                         # zero cost dummy row/col
        s, d = sum(self.supply), sum(self.demand)
        if s > d:
            self.demand.append(s - d)
            for row in self.cost:
                row.append(Fraction(0))
        elif d > s:
            self.supply.append(d - s)
            self.cost.append([Fraction(0)] * len(self.cost[0]))

    def vam(self):
        m, n = self.m, self.n
        supply = list(self.supply)
        demand = list(self.demand)
        alloc = [[Fraction(0)] * n for _ in range(m)]
        basis = set()

        active_rows = set(range(m))
        active_cols = set(range(n))

        # only one cell left so no second cheapest. just use the cost
        def row_penalty(r):
            costs = sorted(self.cost[r][c] for c in active_cols)
            return costs[1] - costs[0] if len(costs) > 1 else costs[0]

        def col_penalty(c):
            costs = sorted(self.cost[r][c] for r in active_rows)
            return costs[1] - costs[0] if len(costs) > 1 else costs[0]

        while active_rows and active_cols:
            best = None
            for r in active_rows:
                p = row_penalty(r)
                if best is None or p > best[0]:
                    best = (p, "row", r)
            for c in active_cols:
                p = col_penalty(c)
                if best is None or p > best[0]:
                    best = (p, "col", c)

            _, kind, idx = best

            if kind == "row":
                r = idx
                c = min(active_cols, key=lambda cc: self.cost[r][cc])
            else:
                c = idx
                r = min(active_rows, key=lambda rr: self.cost[rr][c])

            x = min(supply[r], demand[c])
            alloc[r][c] = x
            basis.add((r, c))
            supply[r] -= x
            demand[c] -= x
            if supply[r] == 0:
                active_rows.discard(r)
            if demand[c] == 0:
                active_cols.discard(c)

        self._repair_degenerate(alloc, basis)
        return alloc, basis

    def _repair_degenerate(self, alloc, basis):
        # cell cant close a loop or the basis isnt a tree and modi stalls
        m, n = self.m, self.n
        while len(basis) < m + n - 1:
            target = None
            for r in range(m):
                for c in range(n):
                    if (r, c) in basis or alloc[r][c] != 0:
                        continue
                    try:
                        self._find_cycle((r, c), basis)
                    except RuntimeError:
                        target = (r, c)          # no loop so its independant
                        break
                if target:
                    break
            if target is None:
                break
            basis.add(target)

    def modi(self, alloc, basis, verbose=True):
        m, n = self.m, self.n
        steps = 0

        while True:
            # u_i + v_j = c_ij over the basis. u_0 = 0. zeros below are a guard
            u = [None] * m
            v = [None] * n
            u[0] = Fraction(0)
            changed = True
            while changed:
                changed = False
                for (i, j) in basis:
                    if u[i] is not None and v[j] is None:
                        v[j] = self.cost[i][j] - u[i]; changed = True
                    elif v[j] is not None and u[i] is None:
                        u[i] = self.cost[i][j] - v[j]; changed = True
            for i in range(m):
                if u[i] is None:
                    u[i] = Fraction(0)
            for j in range(n):
                if v[j] is None:
                    v[j] = Fraction(0)

            entering = None
            worst = Fraction(0)
            for r in range(m):
                for c in range(n):
                    if (r, c) in basis:
                        continue
                    d = self.cost[r][c] - u[r] - v[c]
                    if d < worst:
                        worst, entering = d, (r, c)

            if entering is None:                 # all d >= 0 so were done
                break

            # push theta round the loop. even gain odd give up
            loop = self._find_cycle(entering, basis)
            theta = None
            for k in range(1, len(loop), 2):
                q = alloc[loop[k][0]][loop[k][1]]
                if theta is None or q < theta:
                    theta = q
            for k, (r, c) in enumerate(loop):
                if k % 2 == 0:
                    alloc[r][c] += theta
                else:
                    alloc[r][c] -= theta
            basis.add(entering)
            # only one leaves even if the pivot empties a few
            for k in range(1, len(loop), 2):
                r, c = loop[k]
                if alloc[r][c] == 0:
                    basis.discard((r, c))
                    break
            steps += 1
            if verbose:
                print(f"  MODI iter {steps}: enter ({entering[0]+1},{entering[1]+1})"
                      f"  theta = {theta}")

        return alloc, basis, steps

    def _find_cycle(self, start, basis):
        # row close only. vertical lets in odd loops that break the totals
        path = [start]

        def dfs(cell, came_horizontal):
            r, c = cell
            if came_horizontal:                 # go vertical now
                for i in range(self.m):
                    nxt = (i, c)
                    if nxt in basis and nxt not in path:
                        path.append(nxt)
                        if dfs(nxt, False):
                            return True
                        path.pop()
            else:                               # horizontal now
                for j in range(self.n):
                    nxt = (r, j)
                    if nxt == start and len(path) > 2:
                        return True
                    if nxt in basis and nxt not in path:
                        path.append(nxt)
                        if dfs(nxt, True):
                            return True
                        path.pop()
            return False

        if dfs(start, True):
            return path
        raise RuntimeError("no closed loop found for entering cell")

    def total_cost(self, alloc):
        total = Fraction(0)
        for r in range(self.m):
            for c in range(self.n):
                total += self.cost[r][c] * alloc[r][c]
        return total

    def show(self, alloc, basis=None):
        print("  allocation matrix (rows = sources, cols = destinations):")
        for r in range(self.m):
            cells = []
            for c in range(self.n):
                tag = "*" if basis and (r, c) in basis else " "
                cells.append(f"{alloc[r][c]}{tag}")
            print("    " + "  ".join(f"{x:>6}" for x in cells))


if __name__ == "__main__":
    # case study cement distribution. 3 plants to 4 warehouses
    demand = [10, 40, 20, 30]
    supply = [20, 30, 50]
    cost = [[8, 6, 10, 9],      # plant S1
            [9, 12, 13, 7],     # plant S2
            [14, 9, 16, 5]]     # plant S3

    tp = Transportation(supply, demand, cost)
    only_vam = "--vam-only" in sys.argv

    print("Transportation problem: 3 sources x 4 destinations")
    print("  supply  :", supply)
    print("  demand  :", demand)
    print("  costs   :")
    for row in cost:
        print("           ", row)

    print("\n[Phase 1] Vogel's Approximation Method")
    alloc, basis = tp.vam()
    tp.show(alloc, basis)
    vam_cost = tp.total_cost(alloc)
    print(f"  VAM initial BFS cost = {vam_cost}  ( {float(vam_cost):.2f} )")

    if only_vam:
        sys.exit(0)

    print("\n[Phase 2] MODI optimality test / improvement")
    alloc, basis, iters = tp.modi(alloc, basis)
    tp.show(alloc, basis)
    opt = tp.total_cost(alloc)
    print(f"  MODI iterations      = {iters}")
    print(f"  optimal total cost   = {opt}  ( {float(opt):.2f} )")

    row_sum = [sum(alloc[r][c] for c in range(tp.n)) for r in range(tp.m)]
    col_sum = [sum(alloc[r][c] for r in range(tp.m)) for c in range(tp.n)]
    lst = lambda xs: "[" + ", ".join(str(v) for v in xs) + "]"
    print("\n  check: row sums  =", lst(row_sum), " (supply:", supply, ")")
    print("  check: col sums  =", lst(col_sum), " (demand:", demand, ")")
    print("  basis cells      =", len(basis), " (m+n-1 =", tp.m + tp.n - 1, ")")

    assert row_sum == tp.supply and col_sum == tp.demand, "allocation infeasible"
    assert all(alloc[r][c] >= 0 for r in range(tp.m) for c in range(tp.n))
    assert opt == 800, f"known optimum for this case study is 800, got {opt}"
    print("  self-check       = OK")
