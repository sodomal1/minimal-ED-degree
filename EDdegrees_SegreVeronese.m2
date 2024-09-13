restart
needsPackage "LinearTruncations"
-- the following function computes the generic ED degree of a Segre-Veronese variety
-- L is a list of nonnegative integers, the projective dimensions of the factors
-- D is a list of nonnegative integers, the degrees of the embedding
-- We use the formula of Theorem 3.5
genEDdegree = (L,D) -> (
    k := #L;
    N := sum(L);
    for s in 0..N do ind_s = diagonalMultidegrees(s,k);
    sum(N+1, s -> (-1)^s*(2^(N+1-s)-1)*(N-s)!*(sum apply(ind_s, g -> product(k, l-> if g#l > L#l then 0 else (binomial(L#l+1,g#l)*(D#l)^(L#l-g#l))/((L#l-g#l)!)))))
    )

-- examples
genEDdegree({1},{3})

-- function genEDdegree simplified for binary format tensors (L = {1,...,1})
genEDdegreeBinary = k -> sum(k+1, i-> (-1)^i*(2^(k+1-i)-1)*(k-i)!*binomial(k,i)*2^i)
-- the following list corresponds to Table 1 in the paper
apply(10, k-> genEDdegreeBinary(k+1))

-- the following matrix corresponds to Table 2 in the paper
matrix for m in 1..10 list for n in 1..10 list genEDdegree({m,n},{1,1})

-- the following function computes the ED degree of a Segre-Veronese variety
-- with respect to a Frobenius inner product, using the Friedland-Ottaviani formula of Theorem 3.2
fEDdegSegre = (L,D) -> (
    k := #L;
    PP := QQ[z_1..z_k];
    f := product(k,i->((sum(k, j -> D#j*z_(j+1))-z_(i+1))^(L#i+1)-z_(i+1)^(L#i+1))//(sum(k, j -> D#j*z_(j+1))-2*z_(i+1)));
    sub(contract(product(k,j->z_(j+1)^(L#j)),f),QQ)
    )

-- examples
fEDdegSegre({1},{3})

