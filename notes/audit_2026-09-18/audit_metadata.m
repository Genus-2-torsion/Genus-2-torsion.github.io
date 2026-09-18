SetColumns(0);SetMemoryLimit(8*10^9);R<x>:=PolynomialRing(Rationals());n:=0;errors:=0;
id:="1.a";f:=R![Rationals()|-14580,51520,-63684,31054,-3823,-637,-18];h:=R![Rationals()|1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-18,745,-7278,-9032,-100,1271,-168],R![Rationals()|0,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 1200662898753097226511006037;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 15 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="1.b";f:=R![Rationals()|-6,-33,-14,95,-50,3,-1];h:=R![Rationals()|1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-6,33,-14,-95,-50,-3,-1],R![Rationals()|1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -390963;
gg:=0;for p in [5, 7, 11, 13, 17, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 15 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.a";f:=R![Rationals()|-608,389,22,-40,0,1];h:=R![Rationals()|1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-18,94,-129,2,10,1],R![Rationals()|1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -295;
gg:=0;for p in [3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 14 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.b";f:=R![Rationals()|-56,0,-75,0,15,0,-1];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-56,0,-75,0,15,0,-1],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -172032;
gg:=0;for p in [5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 48 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="3.a";f:=R![Rationals()|0,0,-1,2,-2,1];h:=R![Rationals()|1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-1,-2,-2,-1],R![Rationals()|0,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 997;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 3 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="3.b";f:=R![Rationals()|-54,-318,-479,92,302,-5,-54];h:=R![Rationals()|1,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-54,-318,-479,92,302,-5,-54],R![Rationals()|1,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -405000000000000000000000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 21 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="4.a";f:=R![Rationals()|0,-1,0,1];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-1,0,1],R![Rationals()|1,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -2140;
gg:=0;for p in [3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 4 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="4.b";f:=R![Rationals()|-7536,-4824,-4239,-744,-304,76,-4];h:=R![Rationals()|0,0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
gg:=0;for p in [7, 11, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 16 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="5.a";f:=R![Rationals()|-6,11,-19,14,-9,1];h:=R![Rationals()|1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-6,11,-19,14,-9,1],R![Rationals()|1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 277;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 15 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="5.b";f:=R![Rationals()|-344,348,-334,17,14,-36,-8];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-344,348,-334,17,14,-36,-8],R![Rationals()|]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 28098410075794702336;
gg:=0;for p in [3, 5, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 25 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="6.a";f:=R![Rationals()|3,21,46,26,-12,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-1,-12,-27,46,-21,3],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 1038;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 12 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="6.b";f:=R![Rationals()|-92,42,3,37,-12,-3,-5];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-92,42,3,37,-12,-3,-5],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -765625000000;
gg:=0;for p in [3, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 36 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="7.a";f:=R![Rationals()|-2,3,0,-3,0,1];h:=R![Rationals()|0,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-3,-3,2,1,-2,-1],R![Rationals()|1,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 461;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 7 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="7.b";f:=R![Rationals()|-3,-4,2,7,2,-4,-3];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-3,12,-20,14,-10,3,-1],R![Rationals()|1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -5537792;
gg:=0;for p in [3, 5, 7, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 21 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="8.a";f:=R![Rationals()|0,0,0,-1,-2,-2,-1];h:=R![Rationals()|1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,1,1,2,1,1],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 464;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 16 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="8.b";f:=R![Rationals()|-15,0,7,0,-4];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-15,0,7,0,-4],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 58593750000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 48 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="9.a";f:=R![Rationals()|0,0,0,0,-1];h:=R![Rationals()|1,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,0,-1,-1,0,-1],R![Rationals()|1,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 713;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 29, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 9 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="9.b";f:=R![Rationals()|-3,1,-5,1,-2,-1];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-3,1,-5,1,-2,-1],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 376367048;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 9 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="10.a";f:=R![Rationals()|7,16,0,-8,-2,1];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-1,-8,4,3,-2],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 389;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 20 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="10.b";f:=R![Rationals()|-2,-1,10,7,-13,-11];h:=R![Rationals()|0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-2,11,-15,-3,6,2],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -43923;
gg:=0;for p in [5, 7, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 40 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="11.a";f:=R![Rationals()|0,0,1];h:=R![Rationals()|1,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,0,0,-1,1,-1],R![Rationals()|1,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -353;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 11 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="11.b";f:=R![Rationals()|9471400,94188680,99564424,-65405928,31412066,-5109634,1204142];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
gg:=0;for p in [3, 7, 13, 17, 19, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 11 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="12.a";f:=R![Rationals()|1,1,1];h:=R![Rationals()|0,1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,1,1],R![Rationals()|0,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -3048;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 24 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="12.b";f:=R![Rationals()|0,0,1,-1,1];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,0,1,-1,1],R![Rationals()|1,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -294;
gg:=0;for p in [5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 48 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="13.a";f:=R![Rationals()|0,0,-1,-1];h:=R![Rationals()|1,1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,0,-1,-1],R![Rationals()|1,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 349;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 13 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="14.a";f:=R![Rationals()|0,1,1];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,0,0,-1,1,-1],R![Rationals()|1,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 249;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 28 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="14.b";f:=R![Rationals()|0,-4,1,-3,1,-1];h:=R![Rationals()|0,0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-4,1,-3,1,-1],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 1124864;
gg:=0;for p in [3, 5, 7, 11, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 14 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="15.a";f:=R![Rationals()|0,-1,-1];h:=R![Rationals()|1,1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-1,-1],R![Rationals()|1,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 277;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 15 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="15.b";f:=R![Rationals()|1,0,2,0,1];h:=R![Rationals()|0,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,0,2,0,1],R![Rationals()|0,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -1936;
gg:=0;for p in [3, 5, 7, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 15 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="16.a";f:=R![Rationals()|1,1,-2,0,1,-1];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,-1,-1,-1,-1],R![Rationals()|0,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -6640;
gg:=0;for p in [3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 16 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="16.b";f:=R![Rationals()|-125,0,97,0,-26,0,2];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-125,0,97,0,-26,0,2],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 911250000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 48 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="17.a";f:=R![Rationals()|1,3,1,1];h:=R![Rationals()|0,1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,3,1,1],R![Rationals()|0,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 510976;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 17 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="18.a";f:=R![Rationals()|0,2,4,0,-2];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,0,-2,-1,4,-2],R![Rationals()|1,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -18880;
gg:=0;for p in [3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 18 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="18.b";f:=R![Rationals()|-6,5,-7,2,-1,-1];h:=R![Rationals()|0,0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-6,5,-7,2,-1,-1],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 2352980;
gg:=0;for p in [3, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 18 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="19.a";f:=R![Rationals()|-2,5,-3,-5,1,1,1];h:=R![Rationals()|1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-2,-5,-1,6,3,-2],R![Rationals()|1,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 6012928;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 19 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="19.b";f:=R![Rationals()|0,1,1];h:=R![Rationals()|1,0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,1,1],R![Rationals()|1,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -169;
gg:=0;for p in [3, 5, 7, 11, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 19 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="20.a";f:=R![Rationals()|0,0,0,0,0,-1];h:=R![Rationals()|1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,1,0,0,0,-1],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 3152;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 20 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="20.b";f:=R![Rationals()|-2,-7,0,3,-2,1];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-2,-7,0,3,-2,1],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -7073843073;
gg:=0;for p in [5, 7, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 40 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="21.a";f:=R![Rationals()|0,1,2,0,-1];h:=R![Rationals()|1,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,0,-2,1,0,-1],R![Rationals()|1,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 776;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 21 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="21.b";f:=R![Rationals()|0,0,1,2,2,1];h:=R![Rationals()|1,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,0,1,2,2,1],R![Rationals()|1,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -648;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 21 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="22.a";f:=R![Rationals()|1,-1,-2,1];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,-1,-2,1],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -19072;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 22 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="23.a";f:=R![Rationals()|1,-2,2,1,-2,-1,1];h:=R![Rationals()|0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,-2,2,1,-2,-1,1],R![Rationals()|0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -1736704;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 23 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="24.a";f:=R![Rationals()|0,2,4,3,2];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,0,2,-4,4,-2],R![Rationals()|1,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 183168;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 48 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="24.b";f:=R![Rationals()|375,0,97,0,8];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|375,0,97,0,8],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -3750000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 48 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="25.a";f:=R![Rationals()|2,-23,41,75,25,-9];h:=R![Rationals()|1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|2,-23,41,75,25,-9],R![Rationals()|1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -896806687500000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 25 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="25.b";f:=R![Rationals()|1,8,13,9,14,4,2];h:=R![Rationals()|0,0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|2,4,14,9,13,8,1],R![Rationals()|1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 3413194702848;
gg:=0;for p in [5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 25 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="26.a";f:=R![Rationals()|0,0,-3,-3,-3,-12];h:=R![Rationals()|1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-12,-3,-3,-3],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -55507949568;
gg:=0;for p in [5, 7, 11, 13, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 26 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="27.a";f:=R![Rationals()|0,-1,1,1,-1];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-1,1,1,-1],R![Rationals()|1,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 9664;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 27 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="27.b";f:=R![Rationals()|0,4,20,17,9,0,1];h:=R![Rationals()|1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,4,20,17,9,0,1],R![Rationals()|1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 278628139008;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 27 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="28.a";f:=R![Rationals()|2,3,1,1,0,-1];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|2,3,1,1,0,-1],R![Rationals()|1,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -6723;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 28 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="28.b";f:=R![Rationals()|4,10,-4,-1,5,-3,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|4,10,-4,-1,5,-3,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 129008402432;
gg:=0;for p in [3, 5, 11, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 28 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="29.a";f:=R![Rationals()|0,0,-1,2,0,-2,1];h:=R![Rationals()|1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,2,0,-2,-1,-1],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 999424;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 29 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="30.a";f:=R![Rationals()|0,-1,16,-12,-10,4,2];h:=R![Rationals()|0,0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|2,-5,-10,12,16,1],R![Rationals()|1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -1104207552;
gg:=0;for p in [5, 11, 13, 17, 19, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 30 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="30.b";f:=R![Rationals()|2,6,10,11,10,6,2];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|2,7,13,13,9,3,1],R![Rationals()|1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -67744512;
gg:=0;for p in [5, 7, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 60 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="31.a";f:=R![Rationals()|-126,816,-2466,4300,-4587,2841,-839];h:=R![Rationals()|0,-1,-1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-126,60,-276,77,-210,30,-60],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -5339228294700000000000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 31 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="32.a";f:=R![Rationals()|-47,-49,49,45,1,-1,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|9,-27,-55,10,21,-7,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -2759533088735232;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 32 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="33.a";f:=R![Rationals()|2,13,21,7,9,-7,1];h:=R![Rationals()|1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|2,13,21,7,9,-7,1],R![Rationals()|1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -622015552512;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 33 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="34.a";f:=R![Rationals()|15,12,7,3,0,-2];h:=R![Rationals()|0,0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|15,12,7,3,0,-2],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 89973669888;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 34 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="35.a";f:=R![Rationals()|160,80,166,182,116,30,10];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|160,-80,166,-183,116,-30,10],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -1534038085937500000000000;
gg:=0;for p in [3, 7, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 35 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="36.a";f:=R![Rationals()|9,6,-14,13,30,3];h:=R![Rationals()|0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|9,6,-14,13,30,3],R![Rationals()|0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -37186698240000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 36 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="36.b";f:=R![Rationals()|9,-9,-1,-5,5,-1,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|9,-9,-1,-5,5,-1,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 68797071360;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 108 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="39.a";f:=R![Rationals()|0,-1,1,2,1];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-1,1,2,1],R![Rationals()|1,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -214272;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 39 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="40.a";f:=R![Rationals()|0,0,-2,-2,4,3];h:=R![Rationals()|1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-3,4,2,-2,-1],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -1389312;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 40 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="40.b";f:=R![Rationals()|66,0,12,0,2];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|66,0,12,0,2],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -2909269592064;
gg:=0;for p in [5, 7, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 80 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="45.a";f:=R![Rationals()|168300000000,0,49996210000,0,29240200,0,13981];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|85536,351648,435234,181153,125529,41943,13981],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -829580939374927605395181061275152270815887000000000000;
gg:=0;for p in [7, 13, 19, 23, 29, 37, 43, 47, 53, 59, 61, 67] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 45 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="48.a";f:=R![Rationals()|3,0,9,0,-4];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|3,0,9,0,-4],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -1249949232;
gg:=0;for p in [5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 48 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="60.a";f:=R![Rationals()|1,-7,25,-24,25,-7,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,-7,25,-24,25,-7,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -25798901760000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 120 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="63.a";f:=R![Rationals()|-146398496,0,79136353,0,-197570,0,897];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-1054044,4921494,4872774,-96543,-46029,2691,897],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 4519153520051578364207900239722444433179693771085318478364672;
gg:=0;for p in [5, 7, 11, 17, 31, 37, 41, 43, 47, 53, 59, 61] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 63 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="70.a";f:=R![Rationals()|179,425,-51,1];h:=R![Rationals()|110,-41,-3,2];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1089,-231,309,-98,-9,9,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 293063822324235502645248;
gg:=0;for p in [5, 7, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 70 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.a";f:=R![Rationals()|0,1,16,72,33,4];h:=R![Rationals()|0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-4,33,-72,16,-1],R![Rationals()|0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -29696;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 16 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.b";f:=R![Rationals()|-60,-1,59,-1,59,0,-60];h:=R![Rationals()|1,1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-60,-1,59,-1,59,0,-60],R![Rationals()|1,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -11928619582710272;
gg:=0;for p in [3, 5, 7, 11, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 20 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.4.a";f:=R![Rationals()|0,-1,0,16,-8,1];h:=R![Rationals()|0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-1,0,16,-8,1],R![Rationals()|0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 997;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 16 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.4.b";f:=R![Rationals()|-16,6,-6,11,1,3];h:=R![Rationals()|1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-16,6,-6,11,1,3],R![Rationals()|1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 223948800;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.6.a";f:=R![Rationals()|0,0,-2,-1,4,4];h:=R![Rationals()|1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-1,0,-2,-3,1,1],R![Rationals()|1,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -45056;
gg:=0;for p in [3, 5, 7, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 12 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.6.b";f:=R![Rationals()|-23,13,-14,-3,1,-2];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-23,13,-14,-3,1,-2],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 106662334464;
gg:=0;for p in [5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 36 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.8.a";f:=R![Rationals()|0,0,-2,-4,3,8];h:=R![Rationals()|1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-3,2,3,-3,-1,1],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -29696;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 16 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.8.b";f:=R![Rationals()|-214,465,-279,447,-305,60,-175];h:=R![Rationals()|1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-214,465,-279,447,-305,60,-175],R![Rationals()|1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 267422836908796689945715200000;
gg:=0;for p in [11, 13, 17, 19, 23, 31, 37, 41, 43, 47, 53, 59] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.10.a";f:=R![Rationals()|0,1,1,-4,-2,3];h:=R![Rationals()|1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-3,-2,4,1,-2],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 8325;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 20 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.10.b";f:=R![Rationals()|0,0,-1,-1,-1,-1];h:=R![Rationals()|1,1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,0,-1,-1,-1,-1],R![Rationals()|1,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -512;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 20 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.12.a";f:=R![Rationals()|0,-1,2,14,-8,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-1,2,14,-8,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 82296;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 24 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.12.b";f:=R![Rationals()|-5,15,-8,-10,4,2];h:=R![Rationals()|0,1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,4,21,39,21,4],R![Rationals()|1,0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 36450;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 96 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.14.a";f:=R![Rationals()|1,1,0,-1,-2];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,1,0,-1,-2],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -135936;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 28 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.14.b";f:=R![Rationals()|10,-6,-24,-5,-9,-3,9];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|10,-6,-24,-5,-9,-3,9],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -10587979958080896;
gg:=0;for p in [5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 28 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.16.a";f:=R![Rationals()|8,4,-3,1,-1,-1];h:=R![Rationals()|0,0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|8,4,-3,1,-1,-1],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -182145024;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 32 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.16.b";f:=R![Rationals()|2,7,16,18,16,7,2];h:=R![Rationals()|1,1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|2,-8,16,-20,16,-8,2],R![Rationals()|1,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -1934917632;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 32 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.18.a";f:=R![Rationals()|7,2,1,3,-5,-1,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|7,2,1,3,-5,-1,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -18825765984;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 36 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.18.b";f:=R![Rationals()|36,66,12,-5,-3,-3,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|36,66,12,-5,-3,-3,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -104862024990720;
gg:=0;for p in [7, 11, 13, 17, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 36 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.20.a";f:=R![Rationals()|0,-180,-251,185,28,-30,4];h:=R![Rationals()|0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-244,-140,233,78,-62,-6,4],R![Rationals()|1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 17764050124800000000;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 40 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.20.b";f:=R![Rationals()|0,-18,42,-20,15,-3,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-18,42,-20,15,-3,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 41657590808064;
gg:=0;for p in [5, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 40 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.22.a";f:=R![Rationals()|0,-6,12,-5,9,-3,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-6,12,-5,9,-3,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 1151517855744;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 44 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.24.a";f:=R![Rationals()|15,0,-6,0,-1];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|15,0,-6,0,-1],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -25194240;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 96 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.26.a";f:=R![Rationals()|0,0,2,3,-12,-3,9];h:=R![Rationals()|1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|9,3,-12,-3,2],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -9110237184;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 52 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.28.a";f:=R![Rationals()|0,-3,4,34,-42,-3,9];h:=R![Rationals()|1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-20,61,-27,-7,1],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 19491840000000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 56 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.30.a";f:=R![Rationals()|12,-36,18,23,-13,-5,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,-1,-23,57,-23,-1,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 429229228032;
gg:=0;for p in [5, 7, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 60 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.48.a";f:=R![Rationals()|1,-13,19,111,19,-13,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,-13,19,111,19,-13,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 59521392000000000;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 192 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="3.3.a";f:=R![Rationals()|3,3,7,2,4,0,1];h:=R![Rationals()|1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|3,3,7,2,4,0,1],R![Rationals()|1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -123930000;
gg:=0;for p in [7, 11, 13, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 9 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="3.3.b";f:=R![Rationals()|-2,6,-8,6,-8,6,-2];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
gg:=0;for p in [3, 5, 7, 11, 13, 17, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 9 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="3.6.a";f:=R![Rationals()|21,-2,7,-5,0,-1];h:=R![Rationals()|1,1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|21,-2,7,-5,0,-1],R![Rationals()|1,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -6804000000;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 18 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="3.6.b";f:=R![Rationals()|27,-27,56,-13,21,8,6];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|27,-27,56,-13,21,8,6],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -5424906536727048;
gg:=0;for p in [5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 36 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="3.9.a";f:=R![Rationals()|326,-75,-230,-113,40,75,16];h:=R![Rationals()|1,0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|265,165,-159,152,-92,-23,16],R![Rationals()|0,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 16891051025390625000000000;
gg:=0;for p in [7, 11, 17, 23, 29, 31, 41, 43, 47, 53, 59, 61] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 27 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="3.9.b";f:=R![Rationals()|0,-2,0,5,5,3,1];h:=R![Rationals()|1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,1,0,-6,5,-3,1],R![Rationals()|1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -139968;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 27 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="3.12.a";f:=R![Rationals()|-8,12,6,-11,-3,3,1];h:=R![Rationals()|1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-8,12,6,-11,-3,3,1],R![Rationals()|1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 3359232000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 144 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="3.24.a";f:=R![Rationals()|4,6,-4,-6,5,-3,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|4,-6,-4,5,5,3,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -26873856000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 72 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="4.4.a";f:=R![Rationals()|20493,6534,10832,2601,1823,242,99];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|20493,6534,10832,2601,1823,242,99],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -33394123343324496616340163892734375;
gg:=0;for p in [13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 16 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="4.4.b";f:=R![Rationals()|448,0,-7,0,-2,0,7];h:=R![Rationals()|0,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|448,0,-7,0,-2,0,7],R![Rationals()|0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -3603291078374081785958400;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="4.8.a";f:=R![Rationals()|164,682,450,-435,272,-81];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|164,-682,450,434,272,81],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -786611447825421960000000000;
gg:=0;for p in [11, 13, 17, 19, 23, 31, 37, 41, 43, 47, 53, 59] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 32 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="4.8.b";f:=R![Rationals()|0,-15,-5,-5];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-15,-5,-5],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 345600000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="4.12.a";f:=R![Rationals()|9,-27,53,-62,53,-27,9];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|9,27,53,61,53,27,9],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -816293376000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 144 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="4.16.a";f:=R![Rationals()|252,0,-24,0,-4];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|252,0,-24,0,-4],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -2645395200000000;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 128 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="5.5.a";f:=R![Rationals()|-2,4,2,5,2,1];h:=R![Rationals()|1,1,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-2,4,2,5,2,1],R![Rationals()|1,1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 59969536;
gg:=0;for p in [3, 5, 7, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 25 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="5.10.a";f:=R![Rationals()|-135,135,11,-33,15,-3,9];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-135,135,11,-33,15,-3,9],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 62520309325076889600000;
gg:=0;for p in [7, 11, 13, 17, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 100 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="6.6.a";f:=R![Rationals()|0,4875,-1691,6969,-3000,1872];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,4875,-1691,6969,-3000,1872],R![Rationals()|]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 2608421204497982722931257344000000000000;
gg:=0;for p in [7, 11, 23, 31, 37, 41, 43, 47, 53, 59, 61, 67] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 36 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="6.6.b";f:=R![Rationals()|1,3,6,7,6,3,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,3,6,7,6,3,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -21952;
gg:=0;for p in [3, 5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 36 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="6.12.a";f:=R![Rationals()|88825,81950,75207,-13354,-6347,396,132];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|22206,20487,18801,-3339,-1587,99,33],R![Rationals()|1,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -1582766952336404027054696316525120000;
gg:=0;for p in [13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 144 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="7.7.a";f:=R![Rationals()|869675859,0,3232987,0,3025,0,1];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|13639248,202440,203197,1515,760,3,1],R![Rationals()|]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -131055146712245470001669515048389647007744;
gg:=0;for p in [5, 11, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 49 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="8.8.a";f:=R![Rationals()|84285504,-4535664,-4045487,1800118,88597,88596,836];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|21071376,-1133916,-1011372,450029,22149,22149,209],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 52561940892724075764583621835522035428873153802475520000000000;
gg:=0;for p in [13, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 256 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.a";f:=R![Rationals()|0,1,-3,-5,8,10];h:=R![Rationals()|0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-10,8,5,-3,-1],R![Rationals()|0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 338000;
gg:=0;for p in [3, 7, 11, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 8 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.b";f:=R![Rationals()|0,1,32,268,129,16];h:=R![Rationals()|0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-16,129,-268,32,-1],R![Rationals()|0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 243855360;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.4.a";f:=R![Rationals()|0,6,-8,-3,3,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,6,-8,-3,3,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 930852;
gg:=0;for p in [5, 7, 11, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 16 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.4.b";f:=R![Rationals()|0,-15,7,18,0,3,-2];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.6.a";f:=R![Rationals()|0,1,-1,-4,0,3];h:=R![Rationals()|1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-3,0,4,-1,-1],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 39168;
gg:=0;for p in [5, 7, 11, 13, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 24 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.6.b";f:=R![Rationals()|0,1,-4,-3,3,2];h:=R![Rationals()|0,1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,1,-4,-3,3,2],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 18000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 72 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.8.a";f:=R![Rationals()|0,-1,-6,-3,20,9];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|9,-21,3,10,-1,-1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 5586507792;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 32 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.8.b";f:=R![Rationals()|0,1,1,-5,-4,6];h:=R![Rationals()|1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-5,0,7,0,-3],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 6480;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.10.a";f:=R![Rationals()|0,0,-7,12,6,-12];h:=R![Rationals()|1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-12,7,9,-3,-2],R![Rationals()|0,0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 1961791488;
gg:=0;for p in [5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 40 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.10.b";f:=R![Rationals()|-20,-12,66,-2,-45,3,9];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|9,-57,105,-30,-42,-6],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 328628627712;
gg:=0;for p in [5, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 40 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.12.a";f:=R![Rationals()|0,-5,-5,16,8,-7,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|1,7,8,-17,-5,5],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 3321000000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 48 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.12.b";f:=R![Rationals()|1,-9,19,19,-55,-45];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-69,58,44,-27,-11,3,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 273375000;
gg:=0;for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 96 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.14.a";f:=R![Rationals()|-150965100,-169058940,-16976459,1218410,96973,-4020,36];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-37741275,-42264735,-4244115,304602,24243,-1005,9],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 59697750962036741304238387738855507968750000000000;
gg:=0;for p in [11, 17, 29, 31, 41, 43, 47, 53, 59, 61, 67, 71] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 56 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.16.a";f:=R![Rationals()|9,-33,3,61,3,-33,-11];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|9,-33,3,61,3,-33,-11],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 14549673600000;
gg:=0;for p in [13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 128 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.20.a";f:=R![Rationals()|1016064,10501344,15616273,-60534098,27387697,7579404,-1566684];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|254016,2625336,3904068,-15133525,6846924,1894851,-391671],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 21730082908258624452363375871188942218058008075322118963200000000;
gg:=0;for p in [19, 23, 29, 31, 37, 41, 43, 53, 59, 61, 67, 71] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 80 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.24.a";f:=R![Rationals()|0,-158457693276225,-25400003349776,19670508870560,-10404361990400,4237596000000];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-432494559,-2805860817,-3970444347,-2629108779,-143694570,1033698540,145362420],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -891234716919100996302158225554840110141749815144761640901908681598151333676225427520000000000000;
gg:=0;for p in [13, 17, 19, 23, 31, 37, 43, 47, 53, 59, 67, 71] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 192 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.4.4.a";f:=R![Rationals()|0,180,-210,-1850,2204,-116];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,116,2204,1849,-210,-180],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 1564333978866414231947865292800;
gg:=0;for p in [7, 13, 19, 23, 31, 37, 41, 43, 47, 53, 59, 61] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 32 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.4.4.b";f:=R![Rationals()|0,15,40,0,161,-240];h:=R![Rationals()|0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,240,161,0,40,-15],R![Rationals()|0,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -120812547111963240960000;
gg:=0;for p in [7, 11, 13, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.4.8.a";f:=R![Rationals()|0,-959040,522729,4402,15929,240];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-239760,130682,1100,3982,60],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -53071740815830485969233271053205739929600000;
gg:=0;for p in [17, 19, 23, 29, 31, 41, 47, 53, 59, 61, 67, 71] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.4.8.b";f:=R![Rationals()|7524,-1530,-10562,797,3630,140];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|7524,-1530,-10562,797,3630,140],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 490832899010034992656045955481600000;
gg:=0;for p in [11, 17, 19, 23, 29, 31, 41, 43, 47, 53, 59, 61] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 128 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.6.6.a";f:=R![Rationals()|18,-54,9,71,-40,-5,25];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|25,5,-40,-72,9,54,18],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq -18091291148437500000;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 144 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.2.a";f:=R![Rationals()|0,-1,-4,8,23,-27];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,27,23,-9,-4,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 125581640625;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 16 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.2.b";f:=R![Rationals()|26,41,-44,-30,25,-4];h:=R![Rationals()|1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|15,-56,-24,30,5,-4],R![Rationals()|0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 256096265048064;
gg:=0;for p in [5, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 16 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.4.a";f:=R![Rationals()|0,1,-6,-16,146,-225];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,18,-34,5,9,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 180837562500;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 32 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.4.b";f:=R![Rationals()|0,-21,-6,34,16,-3];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-21,-6,34,16,-3],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 189922851562500;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 128 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.6.a";f:=R![Rationals()|2,-12,11,18,-14,-6];h:=R![Rationals()|1,0,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,6,-14,-18,11,12,2],R![Rationals()|0,1,0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 508083840000;
gg:=0;for p in [7, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 48 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.6.b";f:=R![Rationals()|0,30,39,-72,-111,-30];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,30,39,-72,-111,-30],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 10090298369529000000;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 144 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.8.a";f:=R![Rationals()|0,463250390625,463480444900,230082726,28452,1];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|772392807,473592636,-27768796,-120020,14167,112],R![Rationals()|0,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 9376097829696673222178250408886962890625000000000000;
gg:=0;for p in [17, 19, 23, 29, 37, 41, 43, 47, 53, 59, 61, 67] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.8.b";f:=R![Rationals()|-45,42,30,-23,-9,3,1];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-45,-42,30,22,-9,-3,1],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 56710659600;
gg:=0;for p in [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 128 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.10.a";f:=R![Rationals()|0,-2275,191,3955,-767,-1680,576];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|0,-2275,191,3955,-767,-1680,576],R![Rationals()|]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 1004162066644028944180838400000000;
gg:=0;for p in [11, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 80 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.2.12.a";f:=R![Rationals()|0,-57510204598294806758400,9819234969339333136,-603825803380840,16672898913,-211210,1];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|-1282993930035013443975,-397058962729817115,1518598238654317,301623595822,-462983772,-36750,36],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 80555305668714798821273453211978701165384881364366063777249899756695218220841945965848560018846993868427264000000000000;
gg:=0;for p in [29, 31, 37, 43, 47, 53, 59, 61, 67, 71, 79, 83] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 96 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.4.4.a";f:=R![Rationals()|0,300512487407616,414985109760,173387808,26065,1];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|14400,3720,-18600,-3409,3636,252],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 4054938578206555962132640681574400000000;
gg:=0;for p in [13, 23, 29, 37, 41, 43, 47, 53, 59, 61, 67, 71] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.4.4.b";f:=R![Rationals()|4913,867,-5657,-671,1000,60];h:=R![Rationals()|0,1,1];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|4913,867,-5657,-671,1000,60],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 12135124134840201007303198679040000;
gg:=0;for p in [13, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 64 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
id:="2.2.4.8.a";f:=R![Rationals()|0,42938243694863192676456480983270400000000000000,558098950720994323170545736900000000,1821369625722709924500000,2288161182841,1];h:=R![Rationals()|0];C:=HyperellipticCurve(f,h);g:=4*f+h^2;
Cm:=HyperellipticCurve(R![Rationals()|12400419254552180004153,175916808980224671443,666572690814667403,438870061440610,-280020979290,32602500],R![Rationals()|0,1,1]);assert IsIsomorphic(C,Cm);assert Discriminant(Cm) eq 35953487634322382327337810591083846243130877337481149140388002527820026302467504496133458929508822585464045721626841852038845825024000000000000;
gg:=0;for p in [37, 41, 43, 47, 59, 61, 67, 71, 73, 79, 83, 89] do
gp:=PolynomialRing(GF(p))!g;assert IsSquarefree(gp);chi:=R!Reverse(Coefficients(LPolynomial(HyperellipticCurve(gp):Al:="Naive")));gg:=GCD(gg,Integers()!Evaluate(chi,1));end for;
if gg ne 256 then printf "FAIL GCD %o got %o\n",id,gg;errors+:=1;end if;n+:=1;
printf "METADATA_AUDIT_COMPLETE %o curves, %o errors\n",n,errors;quit;