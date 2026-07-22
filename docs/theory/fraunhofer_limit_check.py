"""Validate the stable small-kxx/kyy Fraunhofer edge-sum formulas against
brute-force quadrature and against the current GOAD clamped formula.

A(p,q) = integral over polygon of exp(i(p x + q y)) dx dy

Forms tested:
  - reference: dense triangle quadrature (float64)
  - current:   GOAD's symmetrized edge sum with KXY_EPSILON clamps (float32)
  - stable:    Green's-theorem single-coordinate form with phi series (float32)
  - limit:     Area + i(p Mx + q My) small-argument expansion
"""
import numpy as np

KXY_EPSILON = 1e-3
DIFF_DMIN = 1e-5


def hexagon(radius=3.0, rot=0.3):
    ang = rot + np.arange(6) * np.pi / 3
    return np.stack([radius * np.cos(ang), radius * np.sin(ang)], axis=1)


def reference(verts, p, q, n=400):
    """Dense quadrature: fan-triangulate from centroid, midpoint rule per subtriangle."""
    c = verts.mean(axis=0)
    total = 0.0 + 0.0j
    nv = len(verts)
    for j in range(nv):
        a, b = verts[j], verts[(j + 1) % nv]
        # subdivide triangle (c, a, b) into n^2 congruent pieces via barycentric grid
        for u in range(n):
            for v in range(n - u):
                # two micro-triangles per (u,v) cell except on the diagonal
                for (du, dv, w) in (((u + 1 / 3, v + 1 / 3), None, 1),
                                    ((u + 2 / 3, v + 2 / 3), None, 1) if u + v < n - 1 else (None, None, 0)):
                    if w == 0:
                        continue
                    uu, vv = du
                    x = c + (a - c) * uu / n + (b - c) * vv / n
                    total += np.exp(1j * (p * x[0] + q * x[1]))
        # scale by micro-triangle area at the end
    tri_area = 0.5 * abs(np.cross(verts[0] - c, verts[1] - c))
    # count of micro triangles per big triangle = n^2, each area tri_area/n^2
    # accumulate properly instead: redo cleanly
    return None  # replaced below


def reference2(verts, p, q, n=300):
    """Dense quadrature via uniform barycentric sampling of fan triangles."""
    c = verts.mean(axis=0)
    nv = len(verts)
    total = 0.0 + 0.0j
    for j in range(nv):
        a, b = verts[j], verts[(j + 1) % nv]
        area = 0.5 * abs(np.cross(a - c, b - c))
        # midpoint-rule grid over the triangle
        u, v = np.meshgrid(np.arange(n), np.arange(n))
        mask_low = (u + v) < n
        # lower micro-triangles centroids
        uu1 = (u + 1 / 3)[mask_low] / n
        vv1 = (v + 1 / 3)[mask_low] / n
        # upper micro-triangles centroids
        mask_up = (u + v) < (n - 1)
        uu2 = (u + 2 / 3)[mask_up] / n
        vv2 = (v + 2 / 3)[mask_up] / n
        uu = np.concatenate([uu1, uu2])
        vv = np.concatenate([vv1, vv2])
        xs = c[0] + (a[0] - c[0]) * uu + (b[0] - c[0]) * vv
        ys = c[1] + (a[1] - c[1]) * uu + (b[1] - c[1]) * vv
        total += np.exp(1j * (p * xs + q * ys)).sum() * (area / (n * n))
    return total


def current_goad(verts, p, q):
    """Replicate diff2.rs edge sum incl. clamps, in float32."""
    f32 = np.float32
    p32, q32 = f32(p), f32(q)
    if abs(p32) < KXY_EPSILON:
        p32 = f32(KXY_EPSILON)
    if abs(q32) < KXY_EPSILON:
        q32 = f32(KXY_EPSILON)
    nv = len(verts)
    total = np.complex64(0)
    for j in range(nv):
        xj, yj = f32(verts[j][0]), f32(verts[j][1])
        dx = f32(verts[(j + 1) % nv][0]) - xj
        dy = f32(verts[(j + 1) % nv][1]) - yj
        if abs(dx) < DIFF_DMIN:
            mj = f32(1e6 if dy * dx >= 0 else -1e6)
        else:
            mj = dy / dx
        nj = f32(1e6 if mj >= 0 else -1e6) if abs(mj) < 1e-6 else f32(1.0) / mj
        if abs(dx) < DIFF_DMIN:
            dx = f32(DIFF_DMIN) * np.sign(dx) if dx != 0 else f32(DIFF_DMIN)
        if abs(dy) < DIFF_DMIN:
            dy = f32(DIFF_DMIN) * np.sign(dy) if dy != 0 else f32(DIFF_DMIN)
        delta = p32 * xj + q32 * yj
        delta1 = q32 * mj + p32
        delta2 = p32 * nj + q32
        omega1 = dx * delta1
        omega2 = dy * delta2
        alpha = f32(1.0) / (f32(2.0) * q32 * delta1)
        beta = f32(1.0) / (f32(2.0) * p32 * delta2)
        if not (np.isfinite(alpha) and np.isfinite(beta)):
            continue
        sumim = alpha * (np.cos(delta) - np.cos(delta + omega1)) - beta * (
            np.cos(delta) - np.cos(delta + omega2))
        sumre = -alpha * (np.sin(delta) - np.sin(delta + omega1)) + beta * (
            np.sin(delta) - np.sin(delta + omega2))
        s = np.complex64(complex(sumre, sumim))
        if np.isnan(s):
            continue
        total += s
    # code multiplies by inv_denom and phase later; here compare raw sum,
    # which equals -i * A(p,q)
    return complex(total) * 1j  # undo the -i to compare with A directly


def phi(omega):
    """(e^{i w} - 1)/(i w), series for small w, float32 arithmetic."""
    f32 = np.float32
    w = f32(omega)
    if abs(w) < f32(1e-2):
        return np.complex64(complex(1.0 - w * w / 6.0, w / 2.0 - w**3 / 24.0))
    ew = np.complex64(complex(np.cos(w), np.sin(w)))
    return (ew - np.complex64(1.0)) / np.complex64(1j * w)


def stable(verts, p, q, tiny=None):
    """Single-coordinate Green form in float32; picks the larger of |p|,|q|.

    tiny: threshold on |p|*L,|q|*L below which the Area limit is used.
    """
    f32 = np.float32
    p32, q32 = f32(p), f32(q)
    L = f32(np.abs(verts).max())
    if tiny is None:
        tiny = f32(1e-4)
    if max(abs(p32), abs(q32)) * L < tiny:
        area = 0.0
        mx = my = 0.0
        nv = len(verts)
        for j in range(nv):
            x0, y0 = verts[j]
            x1, y1 = verts[(j + 1) % nv]
            cr = x0 * y1 - x1 * y0
            area += cr
            mx += (x0 + x1) * cr
            my += (y0 + y1) * cr
        area *= 0.5
        mx /= 6.0
        my /= 6.0
        return complex(area) + 1j * complex(p32 * mx + q32 * my)
    nv = len(verts)
    total = np.complex64(0)
    if abs(p32) >= abs(q32):
        for j in range(nv):
            xj, yj = f32(verts[j][0]), f32(verts[j][1])
            dy = f32(verts[(j + 1) % nv][1]) - yj
            dx = f32(verts[(j + 1) % nv][0]) - xj
            delta = p32 * xj + q32 * yj
            omega = p32 * dx + q32 * dy
            ed = np.complex64(complex(np.cos(delta), np.sin(delta)))
            total += np.complex64(dy) * ed * phi(omega)
        return complex(total / np.complex64(1j * p32))
    else:
        for j in range(nv):
            xj, yj = f32(verts[j][0]), f32(verts[j][1])
            dy = f32(verts[(j + 1) % nv][1]) - yj
            dx = f32(verts[(j + 1) % nv][0]) - xj
            delta = p32 * xj + q32 * yj
            omega = p32 * dx + q32 * dy
            ed = np.complex64(complex(np.cos(delta), np.sin(delta)))
            total += np.complex64(dx) * ed * phi(omega)
        return complex(-total / np.complex64(1j * q32))


def main():
    verts = hexagon()
    area = 0.5 * abs(sum(
        verts[j][0] * verts[(j + 1) % 6][1] - verts[(j + 1) % 6][0] * verts[j][1]
        for j in range(6)))
    print(f"hexagon area = {area:.6f}\n")
    cases = [
        ("generic", 1.3, 0.7),
        ("q small (in clamp zone)", 1.3, 5e-4),
        ("q tiny", 1.3, 1e-6),
        ("q zero", 1.3, 0.0),
        ("both in clamp zone", 5e-4, 5e-4),
        ("both tiny (backscatter@180)", 2e-5, -8e-8),
        ("both zero", 0.0, 0.0),
        ("p tiny sign flip A", -2e-5, 3e-4),
        ("p tiny sign flip B", 2e-5, 3e-4),
    ]
    hdr = f"{'case':<30} {'reference':>24} {'current(clamped)':>24} {'stable':>24} {'cur err%':>9} {'stab err%':>10}"
    print(hdr)
    for name, p, q in cases:
        ref = reference2(verts, p, q)
        cur = current_goad(verts, p, q)
        stb = stable(verts, p, q)
        ec = abs(cur - ref) / abs(ref) * 100
        es = abs(stb - ref) / abs(ref) * 100
        print(f"{name:<30} {ref:>24.6f} {cur:>24.6f} {stb:>24.6f} {ec:>8.3f}% {es:>9.5f}%")


if __name__ == "__main__":
    main()
