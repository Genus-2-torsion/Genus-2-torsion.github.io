/* attribution_check.m -- which published curve is each census curve of a disputed attribution?
   (errata of 2026-10-14 in data/knowledge/sources/PROVENANCE.md)

   Compares census curves (models as in data/curves/<id>.json) with the curves printed in the cited
   sources, by G2-invariants and by Magma's IsIsomorphic over Q.  Equations of the sources:
     Leprevost 1995 = F. Leprevost, J. Theor. Nombres Bordeaux 7 (1995), 283-306 (scan, numdam):
       C24,1  p. 291: f = (2x^2 - 2x - 1)(2x^4 - 10x^3 + 7x^2 + 4x - 4)
       C24,2  p. 292: f = (x^2 - x + 1)(x^4 - 3x^3 + 8x^2 - 3x + 1)
       C25    p. 293: f = C^2 - 12x^2, C = -6x^3 + 13x^2 - 6x + 3  (= 36x^6 - 156x^5 + ... + 9, p. 302)
       C27    p. 293: f = C^2 - 36x^2, C = 2x^3 - 15x^2 + 6x - 3 = (2x^3 - 15x^2 - 3)(2x^3 - 15x^2 + 12x - 3)
              p. 302: X = (x-1)^2/(2x^3 - 15x^2 - 3), Y = 2(x^3 - 3x^2 + 15x + 3) y/(2x^3 - 15x^2 - 3)^2
                      maps C27 onto E1: Y^2 = 648X^3 + 297X^2 + 18X + 1 (so J(C27) is split over Q)
       C21    p. 300: f = 4x^6 - 12x^5 + 13x^4 - 6x^3 + 3x^2 - 2x + 1 (constructed in section 2.1)
     Howe 2015 = E. W. Howe, Bull. London Math. Soc. 47 (2015), 127-135, arXiv:1407.2654v3, Table 2,
       order 27: y^2 + (6x^3+3x^2+3x-2)y = 6x  [Leprevost 1995, Thm 1.2.1];  y^2 + (x^3-2x+1)y = x^3  New;
                 y^2 + (2x^3+3x^2-3x+2)y = 6x^3+6  New;  y^2 + (6x^3+9x^2+6x-1)y = -3x^2  New
     HLP 2000, Corollary 5 and the Remark after it: the curve is X0(22).
   Run:  magma -b pipeline/magma/attribution_check.m < /dev/null > data/knowledge/sources/attribution_check.log */
SetColumns(0);
SetMemoryLimit(16*10^9);
Q := Rationals();
P<x> := PolynomialRing(Q);
Cv := func<f, h | HyperellipticCurve(P!f, P!h)>;
Tors := func<C | Invariants(TorsionSubgroup(Jacobian(SimplifiedModel(C))))>;
procedure Compare(a, an, b, bn)
  iso := IsIsomorphic(a, b);
  printf "%o vs %o: same G2-invariants %o, isomorphic over Q %o\n", an, bn, G2Invariants(a) eq G2Invariants(b), iso;
end procedure;

// census curves
c27b := Cv([0,4,20,17,9,0,1], [1,1]);              // 27.b (table2.txt row 22, label 7776.e.2)
c27c := Cv([9,-36,90,-192,249,-60,4], [0]);        // 27.c (issue #8)
c24b := Cv([375,0,97,0,8], [0,1,0,1]);             // 24.b (row 20, label 300.a.1)
c25b := Cv([1,8,13,9,14,4,2], [0,0,1,1]);          // 25.b (row 21, label 42336.dw.1)
c21c := Cv([1,3,3,1], [0,1,1,1]);                  // 21.c (issue #7, LMFDB 676.a.5408.1)
c39b := Cv([0,-2,0,5,5,3,1], [1,1,1]);             // 3.9.b (row 48, label 2916.e.2)
c66b := Cv([1,3,6,7,6,3,1], [0,1,1]);              // 6.6.b (row 57, LMFDB 196.a.21952.1)
c55a := Cv([-2,4,2,5,2,1], [1,1,1,1]);             // 5.5.a (row 55, label 121.a.1)
// the sources
C241 := HyperellipticCurve((2*x^2-2*x-1)*(2*x^4-10*x^3+7*x^2+4*x-4));
C242 := HyperellipticCurve((x^2-x+1)*(x^4-3*x^3+8*x^2-3*x+1));
C25 := HyperellipticCurve((-6*x^3+13*x^2-6*x+3)^2 - 12*x^2);
assert C25 eq HyperellipticCurve(36*x^6 - 156*x^5 + 241*x^4 - 192*x^3 + 102*x^2 - 36*x + 9);
A := 2*x^3-15*x^2-3; B := 2*x^3-15*x^2+12*x-3;
C27 := HyperellipticCurve(A*B);
assert A*B eq (2*x^3-15*x^2+6*x-3)^2 - 36*x^2;
C21 := HyperellipticCurve(4*x^6 - 12*x^5 + 13*x^4 - 6*x^3 + 3*x^2 - 2*x + 1);
H27 := [Cv([0,6], [-2,3,3,6]), Cv([0,0,0,1], [1,-2,0,1]), Cv([6,0,0,6], [2,-3,3,2]), Cv([0,0,-3], [-1,6,9,6])];
X22 := SmallModularCurve(22);
X28 := SmallModularCurve(28);
L39 := Cv([0,1,0,-6,5,-3,1], [1,1,1]);             // production LMFDB 2916.a.139968.1

print "== torsion of the source curves";
for pr in [<"C24,1", C241>, <"C24,2", C242>, <"C25", C25>, <"C27", C27>, <"C21", C21>] cat [<"Howe 27 #" cat IntegerToString(i), H27[i]> : i in [1..4]] do
  printf "J(Q)_tors(%o) = %o\n", pr[1], Tors(pr[2]);
end for;
print "== Leprevost's map from C27 onto E1 (p. 302)";
N := 648*(x-1)^6 + 297*(x-1)^4*A + 18*(x-1)^2*A^2 + A^3;
printf "648 p^3 + 297 p^2 q + 18 p q^2 + q^3 = B (2x^3 - 6x^2 + 30x + 6)^2 for p = (x-1)^2, q = A: %o\n", N eq B*(2*x^3-6*x^2+30*x+6)^2;
E1 := MinimalModel(EllipticCurve([0,297,0,18*648,648^2]));
printf "E1 = %o (conductor %o, torsion %o)\n", CremonaReference(E1), Conductor(E1), Invariants(TorsionSubgroup(E1));

print "== claim 1: census 27.b";
printf "G2(27.b) = %o\nG2(C27)  = %o\n", G2Invariants(c27b), G2Invariants(C27);
Compare(c27b, "27.b", C27, "Leprevost C27");
for i in [1..4] do Compare(c27b, "27.b", H27[i], "Howe 2015 Table 2, order 27, row " cat IntegerToString(i)); end for;
Compare(C27, "Leprevost C27", H27[1], "Howe 2015 Table 2, order 27, row 1");
Compare(c27c, "27.c", C27, "Leprevost C27");
print "== claim 2: census 24.b";
printf "G2(24.b)  = %o\nG2(C24,1) = %o\nG2(C24,2) = %o\n", G2Invariants(c24b), G2Invariants(C241), G2Invariants(C242);
Compare(c24b, "24.b", C241, "Leprevost C24,1");
Compare(c24b, "24.b", C242, "Leprevost C24,2");
print "== rows the review found correct, and the witnesses of 2026-10-14";
Compare(c25b, "25.b", C25, "Leprevost C25");
Compare(c55a, "5.5.a", X22, "X0(22) (HLP 2000, Cor. 5)");
Compare(c21c, "21.c", C21, "Leprevost C21");
print "== database rows of claim 3";
Compare(c39b, "3.9.b", L39, "LMFDB 2916.a.139968.1");
Compare(c66b, "6.6.b", X28, "X0(28)");
print "ATTRIBUTION_CHECK_DONE";
quit;
