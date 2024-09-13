-- Code related to Example 6.12
-- We compute the coefficients of the extended ED polynomial of the Veronese conic in PP^2

restart
n = 2 -- projective dimension
KK = QQ;
R = KK[toList(x_0..x_n)|toList(u_0..u_n)|{t}|(flatten apply(n+1, i-> apply(n+1-i, j-> q_(i,i+j))))]
xx = matrix{{x_0..x_n}}
uu = matrix{{u_0..u_n}}
xu = xx-uu;

-- X is the ideal of the Veronese conic
X = ideal(x_0*x_2-x_1^2)

jacX = diff(xx, transpose gens X);
Xsing = trim(X + minors(codim X, jacX));

-- Q is the symmetric matrix of a general inner product in C^n
Q = matrix for i in 0..n list for j in 0..n list if i==j then q_(min(i,j),max(i,j)) else (1/2)*q_(min(i,j),max(i,j))

-- Dxu is the squared distance function (with respect to Q) between x and u
Dxu = (xu*Q*transpose(xu))_(0,0)
distance = ideal(Dxu-t^2)
time J = minors((codim X)+1, diff(xx, Dxu)||jacX);
time sat = saturate(X+J, Xsing);
time el = eliminate((entries xx)#0, sat+distance);
EDpolynomial = el_0;
EDdegree = degree(t, EDpolynomial)//2

for i in 0..EDdegree do coeff_i = sub(contract(t^(2*i), EDpolynomial),t=>0);
for i in 0..EDdegree do fac_i = factor(coeff_i);

#(fac_4) -- one factor only in q_(i,j). This is the polynomial of degree 6 displayed in Example 6.12
#(fac_3) -- one factor in q_(i,j) and u_k
#(fac_2) -- one factor in q_(i,j) and u_k
#(fac_1) -- one factor only in q_(i,j), another one in q_(i,j) and u_k
#(fac_0) -- one factor only in u_k, another one only in q_(i,j), a third one only in u_k

(quot3,rr3) = quotientRemainder(matrix{{coeff_3}},symmetricPower(2,uu));
(quot2,rr2) = quotientRemainder(matrix{{coeff_2}},symmetricPower(4,uu));

W1 = trim(ideal(coeff_4) + ideal(first entries transpose quot3));
codim W1, degree W1

W2 = trim(ideal(coeff_4) + ideal(first entries transpose quot3) + ideal(first entries transpose quot2));
codim W2, degree W2
time decW2 = decompose W2; -- used 2446.45 seconds
#decW1 -- 2

Y1 = decW2#0
Y2 = decW2#1

codim(Y1), degree(Y1)
codim(Y2), degree(Y2)

det(Q)%Y1 -- = 0, hence the elements of Y1 are not inner products
det(Q)%Y2 -- generically nonzero

-- remarks on the code
-- p_4(u,Q) = p_4(Q) is the polynomial of degree 6 displayed in Example 6.12.
-- p_3(u,Q) is irreducible of bi-degree (2,7) with respect to (u,Q)
-- p_2(u,Q) is irreducible of bi-degree (4,8)
-- p_1(u,Q) is the product of det(Q) times a polynomial of bi-degree (6,6) (hence, the total bi-degree is (6,9))
-- p_0(u,Q) has three factors: the equation of the conic (u_0*u_1-u_2^2)^2 squared, det(Q)^2, and another polynomial of bi-degree (4,4) (hence, the total bi-degree is (8,10)). The last factor in p_0(u,Q) defines the variety (X^\vee \cap Q)^\vee, a union of 4 lines.

-- Let us look at the vanishing of the coefficients of the symbolic ED polynomial.

-- The variety {Q | p_4(Q)=0} in PP(Sym²C³)=PP^5 corresponds to conics tangent to X, or such that ED defect is 1.
-- The variety {Q | p_4(Q)=0 and p_3(u,Q)=0 for all u} has codimension 2 and degree 12. This is the locus of conics Q such that EDdegree_Q(X) = 2. It corresponds to the locus of conics bitangent to X. The Bombieri-Weyl inner product belongs to this variety.
-- The variety {Q | p_4(Q)=0, p_3(u,Q)=0, p_2(u,Q)=0 for all u} has two irreducible components Y_1,Y_2. We verified that
-- codim(Y_1) = 2, deg(Y_1) = 6
-- codim(Y_2) = 3, deg(Y_2) = 4
-- Furthermore, det(Q)=0 if Q is taken in Y_1, while det(Q) is not necessarily 0 for Y_2. This implies immediately that the elements of Y_1 are not inner products.
-- Secondly, the elements of Y_2 should correspond to conics Q such that X and Q meet in a triple point and another point.
-- Again, in this way Q cannot be positive definite.
