Attribute VB_Name = "ModAB"
Option Explicit

' ============================================================
' find_root � Modified Anderson�Bj�rck root finder (UDF)
' Usage: = find_root(x_cell, y_cell, a, b, [target], [xtol], [ytol], [maxiter])
'   x_cell  : cell acting as the variable x
'   y_cell  : cell whose FORMULA computes f(x) using x_cell
'   a, b    : bracket bounds (numbers or cell refs)
'   target  : right-hand side y in f(x) = y (default 0)
' Reference: Ganchovski & Traykov, 2023; Algorithms 2026, 19, 332.
' ============================================================
Public Function find_root(x_cell As Range, y_cell As Range, _
                           ByVal a As Double, ByVal b As Double, _
                  Optional ByVal target As Double = 0#, _
                  Optional ByVal xtol As Double = 0.00000000000001, _
                  Optional ByVal ytol As Double = 0#, _
                  Optional ByVal maxiter As Long = 200) As Variant
    Dim tmpl As String, ws As Worksheet
    Dim x1 As Double, x2 As Double, t As Double
    Dim epsy As Double, epsx As Double
    Dim y1 As Double, y2 As Double, y3 As Double
    Dim x3 As Double, ym As Double, dy As Double
    Dim r As Double, k As Double, m As Double
    Dim threshold As Double
    Dim side As Long, i As Long
    Dim bisection As Boolean

    On Error GoTo failed

    tmpl = BuildTemplate(x_cell.Cells(1, 1), y_cell.Cells(1, 1))
    If Len(tmpl) = 0 Then GoTo failed
    Set ws = y_cell.Worksheet

    x1 = a: x2 = b
    If x2 < x1 Then t = x1: x1 = x2: x2 = t

    epsy = ytol * MaxD(Abs(target), 1#)

    y1 = fEval(tmpl, ws, x1) - target
    If Abs(y1) <= epsy Then find_root = x1: Exit Function

    y2 = fEval(tmpl, ws, x2) - target
    If Abs(y2) <= epsy Then find_root = x2: Exit Function
    
    If (Sgn(y1) = Sgn(y2)) Then find_root = CVErr(xlErrNA): Exit Function
    
    side = 0
    bisection = True
    threshold = x2 - x1

    For i = 1 To maxiter
        If bisection Then
            x3 = (x1 + x2) * 0.5
        Else
            x3 = (x1 * y2 - y1 * x2) / (y2 - y1)
        End If

        epsx = xtol * MaxD(Abs(x3), 1#)
        If (x2 - x1) <= epsx Then find_root = x3: Exit Function

        If bisection Then
            y3 = fEval(tmpl, ws, x3) - target
            ym = (y1 + y2) * 0.5
            dy = y2 - y1
            r = 1# - Abs(ym / dy)
            k = r * r
            If Abs(ym - y3) < k * (Abs(y3) + Abs(ym)) Then
                bisection = False
                threshold = (x2 - x1) * 8#
            End If
        Else
            If x3 <= x1 Then
                x3 = x1: y3 = y1
            ElseIf x3 >= x2 Then
                x3 = x2: y3 = y2
            Else
                y3 = fEval(tmpl, ws, x3) - target
            End If
            threshold = threshold * 0.5
        End If

        If Abs(y3) <= epsy Then find_root = x3: Exit Function

        If (y1 > 0) = (y3 > 0) Then
            If side = 1 Then
                m = 1# - y3 / y1
                If m > 0# Then y2 = y2 * m Else y2 = y2 * 0.5
            ElseIf Not bisection Then
                side = 1
            End If
            x1 = x3: y1 = y3
        Else
            If side = -1 Then
                m = 1# - y3 / y2
                If m > 0# Then y1 = y1 * m Else y1 = y1 * 0.5
            ElseIf Not bisection Then
                side = -1
            End If
            x2 = x3: y2 = y3
        End If

        If (x2 - x1) > threshold Then
            bisection = True
            side = 0
        End If
    Next i

    find_root = CVErr(xlErrNum)         ' did not converge
    Exit Function
failed:
    find_root = CVErr(xlErrValue)       ' bad formula / eval error
End Function

' ---- replace every reference to x_cell in y_cell's formula
'      with the placeholder "{X}" -----------------------------
Private Function BuildTemplate(x_cell As Range, y_cell As Range) As String
    Dim s As String
    s = y_cell.Formula
    If Len(s) = 0 Then Exit Function
    If Left$(s, 1) <> "=" Then Exit Function
    s = Mid$(s, 2)

    ' split x_cell address into column letters + row number
    Dim addr As String, i As Long, ch As String
    Dim col As String, rw As String
    addr = x_cell.Address(False, False)            ' e.g. "AB12"
    For i = 1 To Len(addr)
        ch = Mid$(addr, i, 1)
        If ch >= "0" And ch <= "9" Then
            col = Left$(addr, i - 1)
            rw = Mid$(addr, i)
            Exit For
        End If
    Next i

    Static re As Object
    If re Is Nothing Then
        Set re = CreateObject("VBScript.RegExp")
        re.Global = True
        re.IgnoreCase = True
    End If
    ' (^|non-id char) + optional $ anchors + (lookahead: not id char)
    re.Pattern = "(^|[^A-Za-z0-9_])(\$?" & col & "\$?" & rw & ")(?![0-9A-Za-z_])"

    BuildTemplate = re.Replace(s, "$1{X}")
End Function

Private Function fEval(ByVal tmpl As String, ws As Worksheet, _
                       ByVal x As Double) As Double
    Dim expr As String
    expr = Replace$(tmpl, "{X}", "(" & Trim$(Str(x)) & ")")
    fEval = CDbl(ws.Evaluate(expr))               ' eval in y_cell's sheet
End Function

Private Function MaxD(ByVal a As Double, ByVal b As Double) As Double
    If a > b Then MaxD = a Else MaxD = b
End Function