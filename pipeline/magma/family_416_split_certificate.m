// Split [4,16] infinitude, 2026-10-20.
// Standalone certificate; no external data files or conjectural settings.
// Run from the repository root: magma -b pipeline/magma/family_416_split_certificate.m
SetColumns(0);
SetMemoryLimit(3*10^9);
Q := Rationals();
K<t> := FunctionField(Q);
R<x> := PolynomialRing(K);

function HLP(al, be)
    a1 := (al[3]-al[2])^2/(be[3]-be[2])
        + (al[2]-al[1])^2/(be[2]-be[1])
        + (al[1]-al[3])^2/(be[1]-be[3]);
    b1 := (be[3]-be[2])^2/(al[3]-al[2])
        + (be[2]-be[1])^2/(al[2]-al[1])
        + (be[1]-be[3])^2/(al[1]-al[3]);
    a2 := al[1]*(be[3]-be[2]) + al[2]*(be[1]-be[3])
        + al[3]*(be[2]-be[1]);
    b2 := be[1]*(al[3]-al[2]) + be[2]*(al[1]-al[3])
        + be[3]*(al[2]-al[1]);
    df := ((al[1]-al[2])*(al[1]-al[3])*(al[2]-al[3]))^2;
    dg := ((be[1]-be[2])*(be[1]-be[3])*(be[2]-be[3]))^2;
    A := dg*a1/a2; B := df*b1/b2;
    P<X> := PolynomialRing(Parent(al[1]));
    h := -(A*(al[2]-al[1])*(al[1]-al[3])*X^2
               + B*(be[2]-be[1])*(be[1]-be[3]))
          *(A*(al[3]-al[2])*(al[2]-al[1])*X^2
               + B*(be[3]-be[2])*(be[2]-be[1]))
          *(A*(al[1]-al[3])*(al[3]-al[2])*X^2
               + B*(be[1]-be[3])*(be[3]-be[2]));
    return h, A;
end function;

function Delta(P, roots)
    if IsZero(P) then return [Parent(roots[1])!1 : r in roots]; end if;
    xx := P[1]/P[3];
    return [xx eq roots[i]
        select &*[roots[i]-roots[j] : j in [1..3] | j ne i]
        else xx-roots[i] : i in [1..3]];
end function;

function SameClasses(v, w)
    return &and[IsSquare(v[i]/w[i]) : i in [1..3]];
end function;

m := t^2-1/16;
al := [m,-t/2-1/8,t/2-1/8];
be := [K!(-15/112),K!(-6/49),K!(-10/63)];
assert al[1]*(be[3]-be[2])+al[2]*(be[1]-be[3])
       +al[3]*(be[2]-be[1]) eq -(256*t^2+47*t+16)/7056;
E := EllipticCurve(&*[x-r : r in al]);
F := EllipticCurve(&*[x-r : r in be]);
P := E![0,-m/2,1]; T := E![al[2],0,1];
QF := F![0,5/98,1]; U := F![be[3],0,1];
assert 2*P eq E![al[1],0,1] and 4*P eq E!0;
assert 4*QF eq F![be[2],0,1] and 8*QF eq F!0;
assert QF+U ne F!0 and U ne 4*QF;
assert SameClasses(Delta(P,al),[-m,2*(t+1/4),-2*(t-1/4)]);
assert SameClasses(Delta(T,al),[K!(-1),t,-t]);
assert SameClasses(Delta(QF,be),[K!105,K!6,K!70]);
assert SameClasses(Delta(U,be),[K!(-7),K!(-1),K!7]);
assert SameClasses(Delta(QF+U,be),[K!(-15),K!(-6),K!10]);
print "PASS: elliptic torsion orders, labels and all Kummer triples.";

D := 256*t^2+47*t+16;
f := D*(16*x^2+9*t)*(175*x^2-9*(4*t-1)^2)
     *(9*x^2+(4*t+1)^2);
h,A := HLP(al,be);
assert A eq -25*D/17076987648;
assert Evaluate(h,7*x) eq (A*t*m^2/12096)^2*f;
assert Discriminant(f) eq
    56171432824012800*t*(4*t-1)^2*(4*t+1)^2*D^22;
jt := (256*t^4+224*t^2+1)^3/(t^2*(4*t-1)^4*(4*t+1)^4);
assert jInvariant(E) eq jt;
assert Degree(Numerator(jt)) eq 12;
assert Degree(Denominator(jt)) eq 10;
print "PASS: explicit HLP model, discriminant and nonconstant elliptic j.";

L<s> := FunctionField(Q);
ev := hom<K -> L | -(3*(3*s^2-10*s+5)^2+(3*s^2-5)^2)
                         /(4*(3*s^2-5)^2)>;
aa := 3*(3*s^2-10*s+5)/(2*(3*s^2-5));
bb := -5*(3*s^2-6*s+5)/(2*(3*s^2-5));
ts := ev(t);
assert Evaluate(ts,1) eq -1 and Evaluate(ts,2) eq -19/49;
assert aa^2 eq -3*(ts+1/4);
assert bb^2 eq -5*(ts-1/4);
assert 6*bb^2-10*aa^2 eq 15;
assert ev(m) eq (aa*bb/15)^2*15;
assert SameClasses([ev(c) : c in Delta(P,al)],[L!(-15),L!(-6),L!10]);
assert not IsSquare(L!(-15)) and not IsSquare(L!15);
// First-coordinate classes of the four F torsion cosets.
assert &and[not IsSquare(Q!(-1)/v) : v in [1,-7,105,-15]];
assert &and[not IsSquare(Q!(-15)/v) : v in [1,-7,105]];
print "PASS: conic parametrization, matching and the constant nonsquare obstructions.";

Z := Integers();
relation_matrix := Matrix(Z,3,3,[4,0,0, 0,8,0, -3,3,2]);
smith := SmithForm(relation_matrix);
assert [smith[i,i] : i in [1..3]] eq [1,4,16];
print "PASS: subgroup presentation has invariant factors [4,16].";

QX<X> := PolynomialRing(Q);
FQ := EllipticCurve([Q!(47/42),Q!(5/49),Q!(5/49),Q!0,Q!0]);
TF := TorsionSubgroup(FQ);
assert Invariants(TF) eq [2,8];
for p in [11,13] do
    assert #ChangeRing(FQ,GF(p)) eq 16;
    printf "Fixed elliptic factor: #F(F_%o)=16.\n", p;
end for;
print "PASS: fixed elliptic factor has exact torsion [2,8].";

// Independent exact check of the small seed and its census identification.
seedf := (16*X^2-9)*(7*X^2-9)*(X^2+1);
KX<z> := FunctionField(Q);
assert z^6*Evaluate(X^6-14*X^4-95*X^2+1008,3/z)
       eq 9*Evaluate(seedf,z);
Cseed := HyperellipticCurve(seedf);
Census := HyperellipticCurve(-4*X^4-24*X^2+252,X^3+X);
assert IsIsomorphic(Cseed,Census);
TG := TorsionSubgroup(Jacobian(Cseed));
assert Invariants(TG) eq [4,16];
print "PASS: t=-1 seed is Q-isomorphic to census 4.16.a, exact torsion [4,16].";

// A second parameter is a separate check of the specialization formula.
ts2 := Q!(-19/49);
f2 := QX![Evaluate(c,ts2) : c in Coefficients(f)];
C2 := IntegralModel(HyperellipticCurve(f2));
TG2 := TorsionSubgroup(Jacobian(C2));
assert Invariants(TG2) eq [4,16];
assert not IsIsomorphic(Cseed,C2);
print "PASS: s=2, t=-19/49 has exact [4,16] and is not Q-isomorphic to the seed.";
print "ALL ASSERTIONS PASSED.";
print "CERTIFIED: split [4,16] family identities and exact torsion calculations.";
quit;
