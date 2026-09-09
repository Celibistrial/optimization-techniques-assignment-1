"""big-M simplex. M kept as two fraction rows zM and zc"""

from fractions import Fraction


class BigMSimplex:
    def __init__(self, c, A, relations, b, sense="min"):
        self.n_orig = len(c)
        self.A = [[Fraction(v) for v in row] for row in A]
        self.relations = list(relations)
        self.b = [Fraction(v) for v in b]
        self.sense = sense
        self.c = [Fraction(v) for v in c]
        if sense == "max":  # max c*x is just min -c*x
            self.c = [-v for v in self.c]

    def _build(self):
        m = len(self.A)
        slack_col = [None] * m
        surp_col = [None] * m
        art_col = [None] * m

        n = self.n_orig
        for i, rel in enumerate(self.relations):
            if rel == "<=":
                slack_col[i] = n
                n += 1
            elif rel == ">=":
                surp_col[i] = n
                art_col[i] = n + 1
                n += 2
            else:  # "="
                art_col[i] = n
                n += 1

        names = [f"x{j + 1}" for j in range(self.n_orig)]
        basic = [None] * m
        for i, rel in enumerate(self.relations):
            if rel == "<=":
                names.append(f"s{i + 1}")
                basic[i] = slack_col[i]
            elif rel == ">=":
                names.append(f"e{i + 1}")
                names.append(f"a{i + 1}")
                basic[i] = art_col[i]
            else:
                names.append(f"a{i + 1}")
                basic[i] = art_col[i]

        rows = [[Fraction(0)] * (n + 1) for _ in range(m)]
        for i in range(m):
            for j in range(self.n_orig):
                rows[i][j] = self.A[i][j]
            if slack_col[i] is not None:
                rows[i][slack_col[i]] = Fraction(1)
            if surp_col[i] is not None:
                rows[i][surp_col[i]] = Fraction(-1)
            if art_col[i] is not None:
                rows[i][art_col[i]] = Fraction(1)
            rows[i][-1] = self.b[i]

        # min c*x + M*(artificals)
        zc = [Fraction(0)] * (n + 1)
        zM = [Fraction(0)] * (n + 1)
        for j in range(self.n_orig):
            zc[j] = self.c[j]
        for i in range(m):
            if art_col[i] is not None:
                zM[art_col[i]] = Fraction(1)

        # artificals start basic so clear em from the z row
        for i in range(m):
            if art_col[i] is not None:
                self._eliminate(rows[i], zM, zc, art_col[i])

        return rows, zM, zc, basic, names

    def solve(self, verbose=False):
        rows, zM, zc, basic, names = self._build()
        m = len(self.A)
        n = len(names)

        if verbose:
            self._show(rows, zM, zc, names, basic)

        iterations = 0
        while True:
            # most negative enters. M beats any constant
            entering = None
            best = (Fraction(0), Fraction(0))
            for j in range(n):
                if (zM[j], zc[j]) < best:
                    best, entering = (zM[j], zc[j]), j
            if entering is None:
                break

            leaving = None
            ratio = None
            for i in range(m):
                a = rows[i][entering]
                if a > 0:
                    r = rows[i][-1] / a
                    # on a tie kick the artifical out
                    if (
                        leaving is None
                        or r < ratio
                        or (r == ratio and names[basic[i]].startswith("a"))
                    ):
                        leaving, ratio = i, r
            if leaving is None:
                raise ArithmeticError(
                    "LP unbounded: no positive pivot in entering column."
                )

            if verbose:
                print(
                    f"[iter {iterations}] enter {names[entering]}  "
                    f"leave {names[basic[leaving]]}"
                )
            self._pivot(rows, zM, zc, leaving, entering)
            basic[leaving] = entering
            iterations += 1
            if verbose:
                self._show(rows, zM, zc, names, basic)

        x = [Fraction(0)] * self.n_orig
        for i in range(m):
            if basic[i] < self.n_orig:
                x[basic[i]] = rows[i][-1]

        # artifical still basic and nonzero = infeasible
        feasible = all(
            not (names[basic[i]].startswith("a") and rows[i][-1] != 0) for i in range(m)
        )

        # z row holds -Z. max got negated in __init__ so flip back
        z = zc[-1] if self.sense == "max" else -zc[-1]

        return {"x": x, "z": z, "feasible": feasible, "iterations": iterations}

    @staticmethod
    def _eliminate(pivot_row, zM, zc, col):  # pivot_row has to be M free
        for zrow in (zM, zc):
            f = zrow[col]
            if f:
                for j in range(len(zrow)):
                    zrow[j] -= f * pivot_row[j]

    @classmethod
    def _pivot(cls, rows, zM, zc, leaving, entering):
        p = rows[leaving][entering]
        rows[leaving] = [v / p for v in rows[leaving]]
        for i, row in enumerate(rows):
            if i != leaving and row[entering]:
                f = row[entering]
                rows[i] = [v - f * w for v, w in zip(row, rows[leaving])]
        cls._eliminate(rows[leaving], zM, zc, entering)

    @staticmethod
    def _show(rows, zM, zc, names, basic):
        def cell(mm, cc):
            if not mm:
                return str(cc)
            return f"{mm}M{cc:+}" if cc else f"{mm}M"

        print("  ".join(f"{nm:>7}" for nm in names) + "  RHS")
        print("  z| " + "  ".join(f"{cell(a, b):>7}" for a, b in zip(zM, zc)))
        for i, row in enumerate(rows):
            print(f"{names[basic[i]]:>3}| " + "  ".join(f"{str(v):>7}" for v in row))
        print()


if __name__ == "__main__":
    # case study diet blend. min 4x1+x2 st 3x1+x2=3 and 4x1+3x2>=6 and x1+2x2<=4
    c = [4, 1]
    A = [[3, 1], [4, 3], [1, 2]]
    relations = ["=", ">=", "<="]
    b = [3, 6, 4]

    solver = BigMSimplex(c, A, relations, b, sense="min")
    res = solver.solve(verbose=True)

    print("=" * 46)
    if not res["feasible"]:
        print("Infeasible: an artificial variable could not be driven to 0.")
    else:
        for k in range(len(c)):
            print(f"  x{k + 1} = {res['x'][k]}   ( {float(res['x'][k]):.6f} )")
        print(f"  Z   = {res['z']}   ( {float(res['z']):.6f} )")
        print(f"  simplex iterations = {res['iterations']}")
        assert (res["x"][0], res["x"][1], res["z"]) == (
            Fraction(2, 5),
            Fraction(9, 5),
            Fraction(17, 5),
        ), "known optimum changed"
        print("  self-check         = OK")
