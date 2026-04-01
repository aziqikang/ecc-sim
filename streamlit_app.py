import streamlit as st
import matplotlib.pyplot as plt
import numpy as np
from ecc import PrimeField, Curve, group_order, group_structure, point_order

# Page config
st.set_page_config(
    page_title="Elliptic Curve Explorer",
    page_icon="🔐",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("Elliptic Curve Group Structure Explorer")

# Sidebar: Curve Parameters
st.sidebar.markdown("## Curve Parameters")
st.sidebar.markdown("$$y^2 = x^3 + ax + b \\pmod{p}$$")

# Presets
preset = st.sidebar.selectbox(
    "Load preset curve:",
    options=["Custom", "Curve 97 (non-cyclic)", "Curve 71 (cyclic)", "Curve 5 (small)"],
    index=0
)

# Define presets
presets_dict = {
    "Curve 97 (non-cyclic)": {"p": 97, "a": 2, "b": 3},
    "Curve 71 (cyclic)": {"p": 71, "a": 1, "b": 1},
    "Curve 5 (small)": {"p": 5, "a": 1, "b": 0},
}

if preset != "Custom":
    params = presets_dict[preset]
    p_val = params["p"]
    a_val = params["a"]
    b_val = params["b"]
else:
    p_val = st.sidebar.slider("Prime p", min_value=5, max_value=500, value=97, step=1)
    a_val = st.sidebar.number_input("Coefficient a", value=2, step=1)
    b_val = st.sidebar.number_input("Coefficient b", value=3, step=1)

# Validate and compute
try:
    Fp = PrimeField(p_val)
    E = Curve(Fp, a=a_val, b=b_val)

    # Compute group properties
    N = group_order(E)
    (e1, n1), (e2, n2) = group_structure(E)

    is_valid = True
    error_msg = None
except Exception as e:
    is_valid = False
    error_msg = str(e)

# Main display
if not is_valid:
    st.error(f"Invalid curve parameters: {error_msg}")
else:
    # Header with curve info
    col1, col2, col3 = st.columns(3)
    with col1:
        st.metric("Prime Field", f"F_{p_val}")
    with col2:
        st.metric("Group Order", f"N = {N}")
    with col3:
        hasse_lower = p_val + 1 - 2 * int(np.sqrt(p_val))
        hasse_upper = p_val + 1 + 2 * int(np.sqrt(p_val))
        hasse_ok = hasse_lower <= N <= hasse_upper
        st.metric("Hasse Bound", "✓ Valid" if hasse_ok else "✗ Invalid")

    st.divider()

    # Group Structure
    st.markdown("## Group Structure")
    st.markdown(f"**E(F_{p_val}) ≅ ℤ/{n1} × ℤ/{n2}**")

    col1, col2 = st.columns(2)
    with col1:
        st.markdown(f"### Generator e₁")
        st.write(f"Order: {n1}")
        st.code(f"e1 = {e1}", language="text")

    with col2:
        if n2 > 1:
            st.markdown(f"### Generator e₂")
            st.write(f"Order: {n2}")
            st.code(f"e2 = {e2}", language="text")
        else:
            st.markdown("### Cyclic Group")
            st.write("n₂ = 1 → Group is cyclic")

    st.divider()

    # Visualization
    st.markdown("## Point Visualization")

    if p_val <= 200:
        col1, col2 = st.columns(2)

        with col1:
            st.markdown("### All Points on Curve")

            # Collect all affine points
            points = []
            for x in range(p_val):
                rhs = (x**3 + a_val * x + b_val) % p_val

                # Check if rhs is quadratic residue
                if rhs == 0:
                    points.append((x, 0))
                elif pow(rhs, (p_val - 1) // 2, p_val) == 1:
                    # Find y via brute force (small p)
                    for y in range(p_val):
                        if (y * y) % p_val == rhs:
                            points.append((x, y))
                            break

            # Plot
            fig, ax = plt.subplots(figsize=(6, 6))
            if points:
                xs, ys = zip(*points)
                ax.scatter(xs, ys, alpha=0.6, s=50, color='blue', label='Affine points')

                # Highlight generators
                if e1 and not e1.is_infinity:
                    ax.scatter([e1.x.value], [e1.y.value], color='red', s=150,
                              marker='*', label=f'e₁ (order {n1})', zorder=5)

                if e2 and not e2.is_infinity and n2 > 1:
                    ax.scatter([e2.x.value], [e2.y.value], color='green', s=150,
                              marker='*', label=f'e₂ (order {n2})', zorder=5)

            ax.set_xlim(-1, p_val)
            ax.set_ylim(-1, p_val)
            ax.set_xlabel('x')
            ax.set_ylabel('y')
            ax.set_title(f'E(F_{p_val}): y² = x³ + {a_val}x + {b_val}')
            ax.grid(True, alpha=0.3)
            ax.set_aspect('equal', adjustable='box')
            ax.legend()
            st.pyplot(fig)

        with col2:
            st.markdown("### Curve Properties")
            st.write(f"**Number of affine points:** {len(points)}")
            st.write(f"**Total group order (with ∞):** {len(points) + 1}")

            # Factorization of N
            from ecc.ec.finite import _factor
            factors = _factor(N)

            st.write("**Prime factorization of N:**")
            for prime, exp in sorted(factors.items()):
                if exp > 1:
                    st.write(f"  {prime}^{exp} = {prime**exp}")
                else:
                    st.write(f"  {prime}")

            # Group structure details
            st.write("**Group structure decomposition:**")
            if n2 == 1:
                st.write(f"  ℤ/{n1} (cyclic)")
            else:
                st.write(f"  ℤ/{n1} × ℤ/{n2}")
                st.write(f"  n₁ · n₂ = {n1} · {n2} = {n1 * n2} = N ✓")
    else:
        st.info(f"Visualization skipped for p = {p_val} (large). Use p ≤ 200 for curve plot.")

    st.divider()

    # Interactive Tabs
    tab1, tab2, tab3 = st.tabs(["Point Arithmetic", "Scalar Multiplication", "Learn"])

    with tab1:
        st.markdown("### Point Addition")
        st.write("Select two points and compute their sum.")

        col1, col2 = st.columns(2)

        with col1:
            st.write("**Point P**")
            px = st.number_input("P.x", min_value=0, max_value=p_val-1, value=1, step=1, key="px")
            py = st.number_input("P.y", min_value=0, max_value=p_val-1, value=1, step=1, key="py")

        with col2:
            st.write("**Point Q**")
            qx = st.number_input("Q.x", min_value=0, max_value=p_val-1, value=2, step=1, key="qx")
            qy = st.number_input("Q.y", min_value=0, max_value=p_val-1, value=4, step=1, key="qy")

        try:
            P = E.point(px, py)
            Q = E.point(qx, qy)
            R = P + Q

            col1, col2, col3 = st.columns(3)
            with col1:
                st.write("**P**")
                st.code(str(P), language="text")
            with col2:
                st.write("**Q**")
                st.code(str(Q), language="text")
            with col3:
                st.write("**P + Q**")
                st.code(str(R), language="text")
        except Exception as e:
            st.error(f"Invalid points: {e}")

    with tab2:
        st.markdown("### Scalar Multiplication")
        st.write("Compute [k]P for a point P and scalar k.")

        col1, col2 = st.columns(2)

        with col1:
            st.write("**Point P**")
            sx = st.number_input("P.x", min_value=0, max_value=p_val-1, value=1, step=1, key="sx")
            sy = st.number_input("P.y", min_value=0, max_value=p_val-1, value=1, step=1, key="sy")

        with col2:
            st.write("**Scalar k**")
            k = st.number_input("k", min_value=0, value=5, step=1)

        try:
            P = E.point(sx, sy)
            result = k * P

            # Show order of P
            ord_P = point_order(E, P, N)

            st.write(f"**P** = {P}")
            st.write(f"**ord(P)** = {ord_P}")
            st.write(f"**[{k}]P** = {result}")

            if k % ord_P == 0:
                st.success(f"✓ [k]P = ∞ (k is multiple of ord(P))")
        except Exception as e:
            st.error(f"Invalid input: {e}")

    with tab3:
        st.markdown("### About Elliptic Curve Group Structure")

        st.markdown("""
        #### What is Group Structure?

        Every elliptic curve over a finite field F_p has an abelian group structure.
        By Lagrange's theorem, every group of order N can be decomposed as:

        **E(F_p) ≅ ℤ/n₁ × ℤ/n₂**

        where n₂ | n₁ and n₁ · n₂ = N.

        #### Key Concepts

        - **Group Order (N)**: Total number of points on the curve (including point at infinity)
        - **e₁, e₂**: Generators of cyclic factors with orders n₁, n₂
        - **Hasse's Theorem**: |N - (p+1)| ≤ 2√p constrains the group order
        - **Cyclic Groups**: When n₂ = 1, the group is cyclic (most common in cryptography)

        #### Computing Group Structure

        The algorithm:
        1. Compute N = #E(F_p) using baby-step giant-step (BSGS) for large p
        2. Factorize N into prime powers: N = ∏ q_i^{a_i}
        3. For each prime q, find the exponent r such that the cyclic factors have exponents q^r and q^{a-r}
        4. Find generators e₁ (order n₁) and e₂ (order n₂) via random sampling

        #### Educational Value

        - Understand group theory concretely
        - See why cryptographic curves are carefully chosen (large prime order)
        - Foundation for attacks like Pohlig-Hellman (exploit small factors)
        """)

        if preset != "Custom":
            st.markdown("#### Preset Curves")
            for name, params in presets_dict.items():
                st.write(f"- **{name}**: y² = x³ + {params['a']}x + {params['b']} (mod {params['p']})")

if __name__ == "__main__":
    pass
