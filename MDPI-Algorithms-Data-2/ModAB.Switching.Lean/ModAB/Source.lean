import ModAB.Expressions

/- Arithmetic expressions extracted from the accompanying C# source.
Regenerate with tools/translate_source.py; check with --check. -/
namespace ModAB.Source
open IEEE64

def e1Expr : Expr := (.abs (.input 0))
def e2Expr : Expr := (.abs (.input 1))
def rhoExpr : Expr := (.op .div (.min e1Expr e2Expr) (.max e1Expr e2Expr))
def aExpr : Expr := (.op .sub .one rhoExpr)
def bExpr : Expr := (.op .add .one rhoExpr)
def vExpr : Expr := (.op .div (.op .div aExpr .two) bExpr)
def rExpr : Expr := (.op .sub .one vExpr)
def kExpr : Expr := (.op .mul rExpr rExpr)
def scaleExpr : Expr := (.max (.max e1Expr e2Expr) (.abs (.input 2)))
def p1Expr : Expr := (.op .div (.input 0) scaleExpr)
def p2Expr : Expr := (.op .div (.input 1) scaleExpr)
def p3Expr : Expr := (.op .div (.input 2) scaleExpr)
def hExpr : Expr := (.op .div (.op .add p1Expr p2Expr) .two)
def leftExpr : Expr := (.abs (.op .sub p3Expr hExpr))
def term1Expr : Expr := (.op .mul kExpr (.abs hExpr))
def term2Expr : Expr := (.op .mul kExpr (.abs p3Expr))
def rightExpr : Expr := (.op .add term1Expr term2Expr)
def gapExpr : Expr := (.op .sub rightExpr leftExpr)

noncomputable def decision (gap : ℝ) : Nat :=
  if 1/(2:ℝ)^48 < gap then 3 else if gap < -(1/(2:ℝ)^48) then 2 else 1

end ModAB.Source
