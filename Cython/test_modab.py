"""
Test script for Cython ModAB implementation
Build first: python setup.py build_ext --inplace
"""
import math

try:
    import cymodab
    print("[OK] Cython module imported successfully")
except ImportError as e:
    print(f"[FAIL] Failed to import Cython module: {e}")
    print("  Build it first: python setup.py build_ext --inplace")
    exit(1)

# Test cases
def test_sqrt2():
    """Find sqrt(2) using x^2 - 2 = 0"""
    def f(x):
        return x * x - 2

    root = cymodab.modAB_root(f, 0.0, 2.0)
    expected = math.sqrt(2)
    error = abs(root - expected)
    print(f"Test sqrt(2): root={root:.15f}, expected={expected:.15f}, error={error:.2e}")
    assert error < 1e-14, f"Error too large: {error}"

def test_sin():
    """Find pi using sin(x) = 0 in [3, 4]"""
    root = cymodab.modAB_root(math.sin, 3.0, 4.0)
    expected = math.pi
    error = abs(root - expected)
    print(f"Test sin(x)=0: root={root:.15f}, expected={expected:.15f}, error={error:.2e}")
    assert error < 1e-14, f"Error too large: {error}"

def test_cubic():
    """Find root of x^3 - x - 2 = 0 (root ≈ 1.521379707...)"""
    def f(x):
        return x**3 - x - 2

    root = cymodab.modAB_root(f, 1.0, 2.0)
    # Verify by checking f(root) ≈ 0
    y = f(root)
    print(f"Test cubic: root={root:.15f}, f(root)={y:.2e}")
    assert abs(y) < 1e-14, f"Function value too large: {y}"

def test_target_value():
    """Find x where x^2 = 5 (i.e., x^2 - 5 = 0, target y=0)"""
    def f(x):
        return x * x

    root = cymodab.modAB_root(f, 0.0, 3.0, y=5.0)  # Find where f(x) = 5
    expected = math.sqrt(5)
    error = abs(root - expected)
    print(f"Test target y=5: root={root:.15f}, expected={expected:.15f}, error={error:.2e}")
    assert error < 1e-14, f"Error too large: {error}"

if __name__ == "__main__":
    print("Running tests...\n")
    test_sqrt2()
    test_sin()
    test_cubic()
    test_target_value()
    print("\n[OK] All tests passed!")
