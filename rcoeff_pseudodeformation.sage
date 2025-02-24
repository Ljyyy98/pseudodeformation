# Compute the pseudodeformation ring using SageMath

# Define a polynomial ring over Z with variables x_i, y_i, z_i, w_i

R.<x1,x2,x3,x4,y1,y2,y3,y4,z1,z2,z3,z4,w1,w2,w3,w4>=ZZ[]
# Variables:
#  - x_i: Entries of the universal matrix representation of generator `a`
#  - y_i: Entries for `b`
#  - z_i: Entries for `c`
#  - w_i: Entries for `d`

J1 = ideal([
    4*x1^3 + x1^4 + 3*x1^2*(2 + x2*x3) + x2*x3*(6 + x2*x3 + 4*x4 + x4^2) 
    + 2*x1*(2 + x2*x3*(4 + x4)),

    x2*(2 + x1 + x4)*(2 + 2*x1 + x1^2 + 2*x2*x3 + 2*x4 + x4^2),
    x3*(2 + x1 + x4)*(2 + 2*x1 + x1^2 + 2*x2*x3 + 2*x4 + x4^2),

    x2^2*x3^2 + x4*(4 + 6*x4 + 4*x4^2 + x4^3) + x2*x3*(6 + x1^2 + 8*x4 + 3*x4^2 + 2*x1*(2 + x4)),

    2*y1 + y1^2 + y2*y3, y2*(2 + y1 + y4), y3*(2 + y1 + y4), 
    y2*y3 + y4*(2 + y4),

    2*w1 + w1^2 + w2*w3, w2*(2 + w1 + w4), w3*(2 + w1 + w4), 
    w2*w3 + w4*(2 + w4)
])

# Relations a^4 = b^2 = d^2 = 1
J2 = ideal([
    (1 + x1)^2 + x2*x3 - (1 + z1)^2 - z2*z3,
    x2*(2 + x1 + x4) - z2*(2 + z1 + z4),
    x3*(2 + x1 + x4) - z3*(2 + z1 + z4),
    x2*x3 + (1 + x4)^2 - z2*z3 - (1 + z4)^2,

    -1 + (1 + x1) * ((1 + y1)*(1 + x1 + y1 + x1*y1 + x3*y2) + (x2 + x2*y1 + y2 + x4*y2)*y3)
    + x3*(x2*(1 + y1)*(1 + y4) + y2*(2 + x1 + x4 + y1 + x1*y1 + x3*y2 + y4 + x4*y4)),

    (x2 + x2*y1 + y2 + x4*y2) * (2 + x1 + x4 + y1 + x1*y1 + x3*y2 + x2*y3 + y4 + x4*y4),
    (x3 + y3 + x1*y3 + x3*y4) * (2 + x1 + x4 + y1 + x1*y1 + x3*y2 + x2*y3 + y4 + x4*y4),

    -1 + (1 + x4) * (y2*(x3 + y3 + x1*y3 + x3*y4) + (1 + y4)*(1 + x4 + x2*y3 + y4 + x4*y4))
    + x2*(x3*(1 + y1)*(1 + y4) + y3*(2 + x1 + x4 + y1 + x1*y1 + x2*y3 + y4 + x4*y4))
])

# Relations c^2 = a^2, baba = 1
J3 = ideal([
    -(x3*z2) + x2*z3, -(x2*z1) + x1*z2 - x4*z2 + x2*z4,
    x3*z1 - x1*z3 + x4*z3 - x3*z4, x3*z2 - x2*z3,

    w3*x2 - w2*x3, w2*x1 - w1*x2 + w4*x2 - w2*x4,
    -(w3*x1) + w1*x3 - w4*x3 + w3*x4, -(w3*x2) + w2*x3,

    -(y3*z2) + y2*z3, -(y2*z1) + y1*z2 - y4*z2 + y2*z4,
    y3*z1 - y1*z3 + y4*z3 - y3*z4, y3*z2 - y2*z3,

    w3*y2 - w2*y3, w2*y1 - w1*y2 + w4*y2 - w2*y4,
    -(w3*y1) + w1*y3 - w4*y3 + w3*y4, -(w3*y2) + w2*y3,

    2*w1 + w1^2 + w2*w3 - 2*x1 - x1^2 - x2*x3 + 2*w1*z1 + w1^2*z1 - 2*x1*z1 - x1^2*z1 - x2*x3*z1 
    + w3*z2 + w1*w3*z2 + w2*z3 + w1*w2*z3 - 2*x2*z3 - x1*x2*z3 - x2*x4*z3 + w2*w3*z4,

    2*w2 + w1*w2 + w2*w4 - 2*x2 - x1*x2 - x2*x4 + w2*z1 + w1*w2*z1 + w1*z2 + w4*z2 + w1*w4*z2 
    - 2*x1*z2 - x1^2*z2 - x2*x3*z2 + w2^2*z3 + w2*z4 + w2*w4*z4 - 2*x2*z4 - x1*x2*z4 - x2*x4*z4
])

# Relations: ac = ca, ad = da, bc = cb, bd = db, dcd = a^2 c
I = J1 + J2 + J3

# Check if tr(a^2) - 2 is in I
assert (x1^2 + x4^2 + 2*x1 + 2*x4 + 2*x2*x3) in I
print("tr(a^2) - 2 is in I: True")

# Double-check: tr(c^2) - 2 is in I
assert (z1^2 + z4^2 + 2*z1 + 2*z4 + 2*z2*z3) in I
print("tr(c^2) - 2 is in I: True")

