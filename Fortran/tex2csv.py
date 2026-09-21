"""Convert the roots-fortran LaTeX results table into CSV.

`fpm test` regenerates `table_<kind>.tex` -- a wide, colour-coded longtable of
function evaluation counts for every method on every test problem. This script
turns that table into a plain CSV for analysis, dropping the f(x) column (the
LaTeX maths is not useful as data).

Usage, from the Fortran directory:

    python tex2csv.py                      # results/root-fortran output/table_real64.tex
    python tex2csv.py in.tex [out.csv]

Output columns:

    no               problem number
    n                problem parameter, blank when the problem has none
    interval         bracket, de-LaTeX'd, e.g. "[1 + 1e-9, 4 - 1e-9]"
    root             the root, as printed by the test program
    bisect_fallback  ";"-separated methods that failed on this row and fell
                     back to bisection (the "*" footnote in the table)
    <method>...      evaluation count, one column per method, in table order

Three things worth knowing about the source table:

  * `no` is not a key. Parameterised problems repeat the same number with
    different `n`, so 157 problem numbers span 223 rows. Key on
    (no, n, interval).
  * The "*" failure marker cannot be recovered from the count. Most starred
    cells are >= 1000, but row 49 `rbp` = 1802 is *not* starred, so the flag
    is kept in its own column.
  * An f(x) cell may hold a \\multirow with an `array` environment, whose own
    & and \\\\ are not row structure -- hence the brace-aware splitter below.
"""

import csv
import re
import sys
from pathlib import Path

BS = chr(92)
CELLCOLOR = re.compile(re.escape(BS) + r"cellcolor\{[^}]*\}")
VERB = re.compile(re.escape(BS) + r"verb\|([^|]*)\|")
STAR = "$^*$"

N_LEADING = 5  # no, f(x), n, interval, root -- before the per-method columns

DEFAULT_TEX = Path("results") / "root-fortran output" / "table_real64.tex"


def split_row(line):
    """Split a longtable row on top-level `&`, dropping any trailing rule.

    Tracks brace depth and skips escaped characters so that `&` and `\\\\`
    inside a cell's own maths do not split the row.
    """
    line = line.rstrip()
    for tail in (BS + BS + BS + "hline", BS + "hline", BS + BS):
        if line.endswith(tail):
            line = line[: -len(tail)]
            break

    cells, buf, depth, i = [], [], 0, 0
    while i < len(line):
        ch = line[i]
        if ch == BS and i + 1 < len(line):
            buf.append(line[i:i + 2])
            i += 2
            continue
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        elif ch == "&" and depth == 0:
            cells.append("".join(buf))
            buf = []
            i += 1
            continue
        buf.append(ch)
        i += 1
    cells.append("".join(buf))
    return cells


def clean_interval(cell):
    """Reduce a LaTeX interval to plain text."""
    s = cell.strip().strip("$").strip()
    s = s.replace(BS + "times 10^{", "e{")        # 1 \times 10^{-9} -> 1 e{-9}
    s = re.sub(r"\s*e\{(-?\d+)\}", r"e\1", s)      # 1 e{-9}          -> 1e-9
    s = s.replace(BS + "pi", "pi")
    s = s.replace(BS + ",", " ").replace(BS + " ", " ")
    return re.sub(r"\s+", " ", s).strip()


def find_methods(lines):
    """Read the method names out of the table header."""
    for line in lines:
        if line.lstrip().startswith("No.") and BS + "verb|" in line:
            names = VERB.findall(line)
            if names:
                return names
    raise SystemExit("could not find the table header row")


def convert(tex_path, csv_path):
    lines = Path(tex_path).read_text(encoding="utf-8", errors="replace").splitlines()
    methods = find_methods(lines)

    rows, n_flagged = [], 0
    for line in lines:
        if not re.match(r"\s*\d+\s*&", line):
            continue
        cells = split_row(line)
        want = N_LEADING + len(methods)
        if len(cells) < want:
            raise SystemExit("row %s: %d cells, expected %d"
                             % (cells[0].strip(), len(cells), want))

        no = cells[0].strip()
        evals, failed = [], []
        for name, cell in zip(methods, cells[N_LEADING:want]):
            cell = CELLCOLOR.sub("", cell)
            if STAR in cell:
                failed.append(name)
                cell = cell.replace(STAR, "")
            cell = cell.strip()
            if not re.fullmatch(r"\d+", cell):
                raise SystemExit("row %s: cannot parse cell %r" % (no, cell))
            evals.append(int(cell))

        n_flagged += len(failed)
        rows.append([no,
                     cells[2].strip(),
                     clean_interval(cells[3]),
                     cells[4].strip().strip("$").strip(),
                     ";".join(failed)] + evals)

    with open(csv_path, "w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["no", "n", "interval", "root", "bisect_fallback"] + methods)
        writer.writerows(rows)

    return len(rows), len(methods), n_flagged


def main(argv):
    base = Path(__file__).resolve().parent
    tex = Path(argv[1]) if len(argv) > 1 else base / DEFAULT_TEX
    csv_out = Path(argv[2]) if len(argv) > 2 else tex.with_suffix(".csv")

    if not tex.is_file():
        raise SystemExit("no such file: %s" % tex)

    n_rows, n_methods, n_flagged = convert(tex, csv_out)
    print("%s -> %s" % (tex.name, csv_out.name))
    print("  %d rows, %d methods, %d bisection-fallback cells"
          % (n_rows, n_methods, n_flagged))


if __name__ == "__main__":
    main(sys.argv)
