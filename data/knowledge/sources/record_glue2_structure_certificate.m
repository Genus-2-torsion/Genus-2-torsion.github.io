// Proof certificates for the necessary torsion groups above 128 in 2-gluings.
// See notes/record_glue2_classification_2026_09_06.md for the mathematical proof.
// No GRH, BSD, rank-parity assumption, or bounded point enumeration is used
// in the genus-two rank-zero certificates.
SetColumns(0); SetSeed(1); SetMemoryLimit(2*10^9);
SetClassGroupBoundGenerators(map< PowerStructure(RngOrd) -> Integers() | O :-> Ceiling(MinkowskiBound(O)) >);
Q:=Rationals(); R<x>:=PolynomialRing(Q);

// A rational Weierstrass point W embeds C in J by P |-> [P-W].
// If J(Q) has exponent two, every rational P is Weierstrass.
g2polys:=[ x*(x^2-1)*(x^2-9), -x*(x^2-1)*(x^2-9),
    x*(x^2+1)*(x^2-2*x-1), -(x^2-1)*(x^2+1)*(x^2-2*x-1),
    -x*(x^2+1)*(x^2-2*x-1), (x^2-1)*(x^2+1)*(x^2-2*x-1) ];
for f in g2polys do
    C:=HyperellipticCurve(f); J:=Jacobian(C);
    S:=TwoSelmerGroup(J); T:=TorsionSubgroup(J); T2:=TwoTorsionSubgroup(J);
    assert #S eq #T2; // raw 2-Selmer bound proves rank zero
    pts:=Chabauty0(J);
    if f in g2polys[1..4] then
        assert #T eq #T2 and &and[n eq 2 : n in Invariants(T)];
        assert &and[P[2] eq 0 or P[3] eq 0 : P in pts];
    else
        assert Invariants(T) eq [2,10];
        assert &and[P[3] eq 0 or Q!(P[1]/P[3]) in {0,1,-1} : P in pts];
    end if;
    printf "GENUS2_PROOF f=%o selmer_order=%o two_torsion_order=%o torsion=%o ALL_POINTS=%o\n",f,#S,#T2,Invariants(T),pts;
end for;

// Low-genus quotients: their full Mordell-Weil groups are computed.
g1polys:=[ x*(x+3)*(x-1), x*(x-3)*(x+1),
    (x^2-1)*(x^2-9), -(x^2-1)*(x^2-9),
    x^4-6*x^2+1, -(x^4-6*x^2+1),
    -x*(x-1)*(3*x^2-3*x+1), x*(x^2-1) ];
allowed:=[ {Q|0,1,-1,3,-3}, {Q|0,-1,3}, {Q|0,1,-1,3,-3},
    {Q|1,-1,3,-3}, {Q|0}, {Q|1,-1}, {Q|0,1,1/2}, {Q|0,1,-1} ];
for i in [1..#g1polys] do
    C:=HyperellipticCurve(g1polys[i]);
    seed:=Points(C : Bound:=10)[1]; // only used to choose an origin
    E,mp:=EllipticCurve(C,seed); lo,hi:=RankBounds(E); assert hi eq 0;
    T,tm:=TorsionSubgroup(E); back:=Inverse(mp);
    pts:={back(tm(t)) : t in T};
    assert &and[P[3] eq 0 or Q!(P[1]/P[3]) in allowed[i] : P in pts];
    printf "GENUS1_PROOF f=%o rank_bounds=%o torsion=%o ALL_POINTS=%o\n",g1polys[i],[lo,hi],Invariants(T),pts;
end for;

K<t>:=FunctionField(Q);
function SameSC(a,b)
    return IsSquare(a/b);
end function;
function DeltaRoots(rr)
    ans:=[];
    for i in [1..3] do
        row:=[rr[i]-rr[j] : j in [1..3]];
        row[i]:=&*[row[j] : j in [1..3] | j ne i];
        Append(~ans,row);
    end for;
    return ans;
end function;

// Full-2, cyclic-6 family: all three diagonal entries are neither +1
// nor -1 modulo squares, by the genus-one certificates above.
al:=[K|0,(t+3)*(t-1)^3,(t-3)*(t+1)^3]; D:=DeltaRoots(al);
diag:=[(t^2-1)*(t^2-9),t*(t+3)*(t-1),-t*(t-3)*(t+1)];
assert &and[SameSC(D[i][i],diag[i]) : i in [1..3]];
vand:=&*[al[i]-al[j] : i,j in [1..3] | i lt j];
assert SameSC(vand,t*(t^2-1)*(t^2-9));

// Full-2, cyclic-8 family: no entry equals +1, and exactly one equals -1
// modulo rational squares. The other eight are ruled out by the curves above.
be:=[K|0,-16*t^4,-(t^2-1)^4]; px:=-4*(t-1)*t*(t+1)^3;
rawv:=[px-r : r in be]; raww:=DeltaRoots(be)[2];
a:=t^2-2*t-1; b:=t^2+2*t-1; q:=a*b;
v:=[-t*(t^2-1),-t*(t^2+1)*a,(t^2-1)*(t^2+1)*a];
w:=[K|-1,-q,q]; vw:=[t*(t^2-1),t*(t^2+1)*b,(t^2-1)*(t^2+1)*b];
assert &and[SameSC(rawv[i],v[i]) and SameSC(raww[i],w[i]) and
    SameSC(v[i]*w[i],vw[i]) : i in [1..3]];
assert -v[2] eq t*(t^2+1)*(t^2-2*t-1);
assert -vw[2] eq (-t)*((-t)^2+1)*((-t)^2-2*(-t)-1);
assert -v[3] eq -(t^2-1)*(t^2+1)*(t^2-2*t-1);
assert -vw[3] eq -((-t)^2-1)*((-t)^2+1)*((-t)^2-2*(-t)-1);

// Cyclic-12 self-gluing: x(P12) is neither square nor d times a square.
p:=-t*(t-1)*(3*t^2-3*t+1);
d:=(6*t^2-6*t+1)*(2*t^2-2*t+1);
assert SameSC(-4*(t-1)*(3*t^2-3*t+1)*t^5,p);
s:=2*t-1;
assert 64*p*d eq (1-s^4)*(9*s^4-1);
assert -( (3*s^2)^2-1 )*( (3*s^2)^2-9 ) eq 24^2*p*d;
A12:=24*t^8-96*t^7+216*t^6-312*t^5+288*t^4-168*t^3+60*t^2-12*t+1;
B12:=16*(3*t^2-3*t+1)^2*(t-1)^6*t^6;
assert SameSC(A12^2-4*B12,d);
assert Evaluate(Numerator(A12^2-4*B12),Q!1/2) eq 0;
print "SYMBOLIC_KUMMER_IDENTITIES_PASSED";
print "MIXED_192_EXCLUDED: full26 diagonal entries are never +/-1; full28 has one -1 and no +1";

// Universal square-class calculation for the full-2/cyclic-6 block.
// A row permutation different from the root permutation forces either
// a square diagonal entry or a signed square Vandermonde.
V:=VectorSpace(GF(2),7); neg:=V.1;
function FormalDelta(aa,bb,cc,neg)
    return [[neg+aa+cc,aa,neg+cc],[neg+aa,neg+aa+bb,bb],
        [cc,neg+bb,neg+bb+cc]];
end function;
DE:=FormalDelta(V.2,V.3,V.4,neg); DF:=FormalDelta(V.5,V.6,V.7,neg);
for perm in SymmetricGroup(3) do
    if Order(perm) eq 1 then continue; end if;
    rel:=sub<V|[DE[i][j]+DF[i^perm][j] : i,j in [1..3]]>;
    if Order(perm) eq 2 then
        assert &or[DE[i][i] in rel : i in [1..3]];
    else
        assert (V.2+V.3+V.4) in rel or (neg+V.2+V.3+V.4) in rel;
    end if;
end for;
print "ROW_ROOT_SYNCHRONIZATION_PASSED: all five nonidentity permutations excluded";

// Abstract torsion calculation on the half-torsion preimage of the gluing.
// psi acts on E[2]; phi matches bases of the two torsion Kummer images.
function FullMatchGroup(a,b,psi,phi)
    B:=AbelianGroup([2^(a+1),4,2^(b+1),4]);
    e2:=[2^a*B.1,2*B.2]; f2:=[2^b*B.3,2*B.4];
    ker:=[e2[i]+&+[Integers()!psi[i,j]*f2[j] : j in [1..2]] : i in [1..2]];
    gens:=[2*B.i : i in [1..4]] cat
        [B.i+&+[Integers()!phi[i,j]*B.(j+2) : j in [1..2]] : i in [1..2]];
    U:=sub<B|gens>;
    G:=quo<U|[U!z : z in ker]>;
    return Invariants(G);
end function;
GL:=GL(2,2); id:=Identity(GL); e:=Vector(GF(2),[1,0]); f:=Vector(GF(2),[0,1]);
counts:=[0,0,0];
for psi in GL do
    for phi in GL do
        if psi eq phi then
            assert FullMatchGroup(1,1,psi,phi) eq [2,2,2,2]; counts[1]+:=1;
        end if;
        // The unique -1 entry on E_[2,6] is off the diagonal.
        if e*psi^-1 ne f*phi^-1 then
            assert FullMatchGroup(1,3,psi,phi) eq [2,2,16]; counts[2]+:=1;
        end if;
        if e*psi eq e and f*phi eq f then
            assert FullMatchGroup(3,3,psi,phi) eq [2,2,8,8]; counts[3]+:=1;
        end if;
    end for;
end for;
assert counts eq [6,24,4];
printf "GROUP_STRUCTURE_PASSED counts=%o full26=[2,2,6,6] mixed_hypothetical_excluded=[2,2,48] full28=[2,2,8,8]\n",counts;
// Rational changes of the marked parameter allow the two surviving
// full-2 blocks to be written with the identity root matching.
h:=(t+3)/(t-1); scale:=16/(t-1)^4;
ah:=[K|0,(h+3)*(h-1)^3,(h-3)*(h+1)^3];
assert ah eq [scale*(al[3]-al[3]),scale*(al[2]-al[3]),scale*(al[1]-al[3])];
h:=(t+1)/(t-1); scale:=16/(t-1)^8;
bh:=[K|0,-16*h^4,-(h^2-1)^4];
assert bh eq [scale*be[1],scale*be[3],scale*be[2]];
assert IsSquare(scale);
print "PARAMETER_RELABELING_IDENTITIES_PASSED";
// Genus of the normalized fixed-parent multiquadratic covers. Constants
// depending on the fixed parent are squares over Qbar and do not affect this.
function RadicalCoverGenus(fs)
    factors:=Setseq(&join[Seqset([g[1] : g in Factorization(f)]) : f in fs]);
    signs:=Matrix(GF(2),[[Valuation(f,g) mod 2 : g in factors] cat
        [Degree(f) mod 2] : f in fs]);
    assert Rank(signs) eq #fs;
    branch:=&+[Degree(g) : g in factors]+(IsZero(ColumnSubmatrix(signs,#factors+1,1)) select 0 else 1);
    deg:=2^#fs;
    genus:=1+(-2*deg+branch*(deg div 2)) div 2;
    return genus,branch;
end function;
g,b:=RadicalCoverGenus([x,(x-1)*(x+3),(x+1)*(x-3)]);
assert g eq 5 and b eq 6;
printf "NORMALIZED_COVER_GENUS block=full26 degree=8 branches=%o genus=%o\n",b,g;
g,b:=RadicalCoverGenus([x*(x^2-1),x*(x^2+1)*(x^2-2*x-1),x^4-6*x^2+1]);
assert g eq 13 and b eq 10;
printf "NORMALIZED_COVER_GENUS block=full28 degree=8 branches=%o genus=%o\n",b,g;
print "STRUCTURE_CERTIFICATE_DONE";
print "CERTIFIED gluing structure";
quit;
