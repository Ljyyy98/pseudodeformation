// Compute the pseudodeformation ring over F_2

// Define the extraspecial group of order 32 using its GAP ID
G := SmallGroup(32,49);
num_map := NumberingMap(G); // Numbering map for elements in G
ordG:= #G;

// Define a polynomial ring over GF(2) with 2*ordG variables with R.i as T(g_i) and R.(i+ordG) as D(g_i)
R := PolynomialRing(GF(2), 2*ordG);


// Initialize an ideal I with the relations T(1) = 0 and D(1) = 1
I := ideal<R | R.1, R.(1+ordG)>;

// Initialize the maximal ideal m
m := ideal<R | 0>;

// Generate the maximal ideal in the polynomial ring
for x in [1..2*ordG] do 
    m := m + ideal<R | R.x>; 
end for; 

// Add relations defining the pseudorepresentation structure
for l, n in [1..ordG] do 
    for k in [1..ordG] do 
        if (k @@ num_map) eq (l @@ num_map) * (n @@ num_map) then 
            // Enforce the determinant relation: D(l * n) = D(l) * D(n) + D(l) + D(n)
            I := I + ideal<R | R.(k+ordG) - R.(l+ordG) * R.(n+ordG) - R.(l+ordG) - R.(n+ordG)>;

            // Enforce trace compatibility relations
            for a in [1..ordG] do 
                if (a @@ num_map) eq (n @@ num_map) * (l @@ num_map) then
                    I := I + ideal<R | R.k - R.a>; 
                end if;
            end for;

            // Additional trace relation
            for j in [1..ordG] do 
                if (j @@ num_map) eq (l @@ num_map) / (n @@ num_map) then
                    I := I + ideal<R | R.(n+ordG) * R.j + R.j - R.l * R.n + R.k>;		
                end if;
            end for;
        end if;
    end for;
end for;

// All defining relations for pseudorepresentations are now in I.

// Verify whether the element tr(a^2) (R.2) belongs to certain ideals
print R.2 in m + I;   // True: R.2 is in m + I
print R.2 in m^2 + I; // True: R.2 is in m^2 + I
print R.2 in m^3 + I; // False: R.2 is not in m^3 + I

// Interpretation: R.2 is nonzero in the deformation ring since it is not in m^3 + I.

// Compute dimensions of quotient rings to study the structure of the pseudodeformation ring
print Dimension(quo<R | m + I>);   // Dimension 1
print Dimension(quo<R | m^2 + I>); // Dimension 20
print Dimension(quo<R | m^3 + I>); // Dimension 72
print Dimension(quo<R | m^4 + I>); // Dimension 121
print Dimension(quo<R | m^5 + I>); // Dimension 137

// Check if the fifth power of m is inside I
print m^5 subset I; // True

// Compute the final dimension of the quotient ring
print Dimension(quo<R | I>); // Dimension 137

// The Hilbert Series is already determined by the computations above.
