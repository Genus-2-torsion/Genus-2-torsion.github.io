"""Curated, cited facts used by the census: the bibliography, and for each torsion group and each
side of the census (geometrically simple / geometrically split) what is known about *infinitely
many* realisations.

Nothing here is asserted from memory.  Every entry cites a source in SOURCES, and the sources that
are tables or records are copied verbatim into data/knowledge/sources/ (see PROVENANCE.md there).

Grades of the "infinitely many?" answer (there is no "finitely many" grade: no group is known to be
realised only finitely often, and no group is known not to be realised at all):

  exact   -- certified: infinitely many pairwise non-isomorphic curves whose Jacobian is in the
             class and has J(Q)_tors isomorphic to G *exactly*.  On the simple side this is the
             L1-L3 standard of the lab record (a marked family over a rational base with an exact
             fiber, one geometrically simple fiber, pairwise distinct G2-invariants); the
             certificates are in that record.
  family  -- a positive-dimensional family of curves with J(Q)_tors *containing* G is proven in
             the cited source (infinitely many members in moduli); exactness of the torsion, or
             the class, for infinitely many members is not certified.
  open    -- no positive-dimensional construction is recorded.
"""

from __future__ import annotations

# --------------------------------------------------------------------------- bibliography

SOURCES = {
    "BNSS2026": {
        "short": "Balakrishnan–Najman–Shnidman–Sutherland 2026",
        "cite": "J. S. Balakrishnan, F. Najman, A. Shnidman, A. V. Sutherland, Rational torsion on simple genus two Jacobians, arXiv:2608.28543 (2026).",
        "url": "https://arxiv.org/abs/2608.28543",
        "used_for": "Tables 1 and 2 (the seed curves and their sources); Theorem 3.3 ([2,2,2,12] family), Corollary 4.3 and 4.6 ([2,2,4,4]), Theorem 4.7 ([2,2,2,8]); the uniqueness of the [2,2,20] curve in its family",
    },
    "Genus2Torsion": {
        "short": "Genus2Torsion repository",
        "cite": "J. S. Balakrishnan, F. Najman, A. Shnidman, A. V. Sutherland, Genus2Torsion GitHub repository (verification scripts and machine-readable tables), https://github.com/AndrewVSutherland/Genus2Torsion, 2026.",
        "url": "https://github.com/AndrewVSutherland/Genus2Torsion",
        "used_for": "table1.txt, table2.txt (copied to data/knowledge/sources/), the geometric-simplicity criterion and the explicit covers used as splitness certificates",
    },
    "Najman2026": {
        "short": "Najman 2026 (unpublished)",
        "cite": "F. Najman, Which torsion groups are realized by infinitely many simple Jacobians? The record (torsion_jac/notes/census_infinity_record.md, 2026-08-21, unpublished; Magma scripts and logs in the private repository F-Najman/torsion_jac). Copied to data/knowledge/sources/census_infinity_record.md.",
        "url": "data/knowledge/sources/census_infinity_record.md",
        "used_for": "the 'exact' grade of the ∞ column on the simple side (39 groups), and the notes on the open ones",
    },
    "ThisSite": {
        "short": "this site (Magma log)",
        "cite": "Computation for this census on the Mordell workstation (Magma V2.29-4); script and log in the repository.",
        "url": "data/knowledge/sources/trivial_torsion_family.log",
        "used_for": "the trivial group: fibers of y² = x⁵ + x + t",
    },
    "HLP2000": {
        "short": "Howe–Leprévost–Poonen 2000",
        "cite": "E. W. Howe, F. Leprévost, B. Poonen, Large torsion subgroups of split Jacobians of curves of genus two or three, Forum Math. 12 (2000), 315–364 (arXiv:math/9809210).",
        "url": "https://arxiv.org/abs/math/9809210",
        "used_for": "Theorem 1 / Table 1: families of genus-2 curves over ℚ with split Jacobians whose torsion contains G, for 20 groups (excerpt in data/knowledge/sources/HLP2000_Theorem1_Table1.txt); seed curves of Table 2",
    },
    "BSSVY": {
        "short": "LMFDB (BSSVY 2016)",
        "cite": "A. R. Booker, J. Sijsling, A. V. Sutherland, J. Voight, D. Yasaki, A database of genus-2 curves over the rational numbers, LMS J. Comput. Math. 19(A) (2016), 235–254.",
        "url": "https://www.lmfdb.org/Genus2Curve/Q/",
        "used_for": "seed curves taken from the LMFDB",
    },
    "BookerSutherland": {
        "short": "Booker–Sutherland",
        "cite": "A. R. Booker and A. V. Sutherland, Genus 2 curves of small conductor, in preparation; preliminary data in the alpha version of the LMFDB.",
        "url": "https://alpha.lmfdb.org/Genus2Curve/Q/",
        "used_for": "seed curves taken from the extended database",
    },
    "LMFDB": {
        "short": "LMFDB",
        "cite": "The LMFDB Collaboration, The L-functions and modular forms database, https://www.lmfdb.org, 2026.",
        "url": "https://www.lmfdb.org",
        "used_for": "links from the curves to their LMFDB home pages",
    },
    "Ogg1973": {"short": "Ogg 1973", "cite": "A. Ogg, Rational points on certain elliptic modular curves, in Proc. Sympos. Pure Math. XXIV, AMS, 1973, 221–231.", "url": "", "used_for": "seed curves"},
    "MazurTate1973": {"short": "Mazur–Tate 1973", "cite": "B. Mazur and J. Tate, Points of order 13 on elliptic curves, Invent. Math. 22 (1973), 41–49.", "url": "https://doi.org/10.1007/BF01425572", "used_for": "the [19] split curve J₁(13)"},
    "Flynn1990": {"short": "Flynn 1990", "cite": "E. V. Flynn, Large rational torsion on abelian varieties, J. Number Theory 36 (1990), 257–265.", "url": "https://doi.org/10.1016/0022-314X(90)90036-Q", "used_for": "seed curves"},
    "Flynn1991": {"short": "Flynn 1991", "cite": "E. V. Flynn, Sequences of rational torsions on abelian varieties, Invent. Math. 106 (1991), 433–442.", "url": "https://doi.org/10.1007/BF01243920", "used_for": "seed curves"},
    "Leprevost1991a": {"short": "Leprévost 1991 (order 13)", "cite": "F. Leprévost, Famille de courbes de genre 2 munies d'une classe de diviseurs rationnels d'ordre 13, C. R. Acad. Sci. Paris Sér. I Math. 313 (1991), 451–454.", "url": "", "used_for": "seed curves; the family with a class of order 13"},
    "Leprevost1991": {"short": "Leprévost 1991 (orders 15, 17, 19, 21)", "cite": "F. Leprévost, Familles de courbes de genre 2 munies d'une classe de diviseurs rationnels d'ordre 15, 17, 19 ou 21, C. R. Acad. Sci. Paris Sér. I Math. 313 (1991), 771–774.", "url": "", "used_for": "seed curves; the families with classes of order 15, 17, 19, 21"},
    "Leprevost1993": {"short": "Leprévost 1993", "cite": "F. Leprévost, Points rationnels de torsion de jacobiennes de certaines courbes de genre 2, C. R. Acad. Sci. Paris Sér. I 316 (1993), 819–821.", "url": "", "used_for": "seed curves"},
    "Leprevost1995": {"short": "Leprévost 1995", "cite": "F. Leprévost, Jacobiennes de certaines courbes de genre 2: torsion et simplicité, J. Théor. Nombres Bordeaux 7 (1995), 283–306.", "url": "https://doi.org/10.5802/jtnb.145", "used_for": "seed curves"},
    "Leprevost1997": {"short": "Leprévost 1997", "cite": "F. Leprévost, Sur certains sous-groupes de torsion de jacobiennes de courbes hyperelliptiques de genre g ≥ 1, Manuscripta Math. 92 (1997), 47–63.", "url": "https://doi.org/10.1007/BF02678180", "used_for": "seed curves"},
    "Elkies2002": {"short": "Elkies 2002", "cite": "N. D. Elkies, Curves of genus 2 over ℚ whose Jacobians are absolutely simple abelian surfaces with torsion points of high order, https://people.math.harvard.edu/~elkies/g2_tors.html, 2001–2002, updated 2010.", "url": "https://people.math.harvard.edu/~elkies/g2_tors.html", "used_for": "seed curves; the order-32 component"},
    "Elkies2024": {"short": "Elkies 2024", "cite": "N. D. Elkies, Families of genus-2 curves with 5-torsion, in LuCaNT: LMFDB, Computation, and Number Theory, Contemp. Math. 796, AMS, 2024, 165–185.", "url": "https://doi.org/10.1090/conm/796", "used_for": "the two-parameter [2,2,2,10] family and the 5-torsion chart behind the [2,10], [2,2,10] families"},
    "Howe2015": {"short": "Howe 2015", "cite": "E. W. Howe, Genus-2 Jacobians with torsion points of large order, Bull. London Math. Soc. 47 (2015), 127–135.", "url": "https://doi.org/10.1112/blms/bdu110", "used_for": "seed curves"},
    "PP2012": {"short": "Platonov–Petrunin 2012a", "cite": "V. P. Platonov and M. M. Petrunin, New orders of torsion points in Jacobians of curves of genus 2 over the rational number field, Dokl. Math. 85 (2012), 286–288.", "url": "https://doi.org/10.1134/S1064562412020366", "used_for": "seed curves"},
    "PP2012b": {"short": "Platonov–Petrunin 2012b", "cite": "V. P. Platonov and M. M. Petrunin, On the torsion problem in Jacobians of curves of genus 2 over the rational number field, Dokl. Math. 86 (2012), 642–643.", "url": "https://doi.org/10.1134/S1064562412050146", "used_for": "seed curves"},
    "PP2015": {"short": "Platonov–Petrunin 2015", "cite": "V. P. Platonov and M. M. Petrunin, New curves of genus 2 over the field of rational numbers whose Jacobians contain torsion points of high order, Dokl. Math. 91 (2015), 220–221.", "url": "https://doi.org/10.1134/S1064562415020222", "used_for": "seed curves"},
    "PZP2013": {"short": "Platonov–Zhgun–Petrunin 2013", "cite": "V. P. Platonov, V. S. Zhgun, M. M. Petrunin, On the simplicity of Jacobians for hyperelliptic curves of genus 2 over the field of rational numbers with torsion points of high order, Dokl. Math. 87 (2013), 318–321.", "url": "https://doi.org/10.1134/S1064562413030216", "used_for": "seed curves"},
    "Platonov2014": {"short": "Platonov 2014", "cite": "V. P. Platonov, Number-theoretic properties of hyperelliptic fields and the torsion problem in Jacobians of hyperelliptic curves over the rational number field, Russian Math. Surveys 69 (2014), no. 1, 1–34.", "url": "https://doi.org/10.1070/RM2014v069n01ABEH004877", "used_for": "seed curves"},
    "BFT2014": {"short": "Bruin–Flynn–Testa 2014", "cite": "N. Bruin, E. V. Flynn, D. Testa, Descent via (3,3)-isogeny on Jacobians of genus 2 curves, Acta Arith. 165 (2014), 201–223.", "url": "https://doi.org/10.4064/aa165-3-1", "used_for": "seed curve; Theorem 6, the rational parametrization behind the [3,3] family"},
    "Nicholls2018": {"short": "Nicholls 2018", "cite": "C. Nicholls, Descent methods and torsion on Jacobians of higher genus curves, DPhil thesis, University of Oxford, 2018.", "url": "https://ora.ox.ac.uk/objects/uuid:0a2a3c5f-6c7f-4d0d-9b21-4c9a2d1e6d64", "used_for": "seed curve ([25])"},
    "costa": {"short": "Costa et al. (ModularAbelianSurfaces)", "cite": "E. Costa, N. D. Elkies, S. Hashimoto, A. Jha, K. Martin, B. Poonen, J. Voight, ModularAbelianSurfaces GitHub repository, https://github.com/edgarcosta/ModularAbelianSurfaces, 2022.", "url": "https://github.com/edgarcosta/ModularAbelianSurfaces", "used_for": "the [31] curve (attached to the newform 1830.2.a.q)"},
    "EpochAI": {"short": "Epoch AI 2026", "cite": "Epoch AI, A genus 2 curve over the rationals with a rational torsion point of prime order at least 31, https://epoch.ai/frontiermath/open-problems/genus-2-jacobian-torsion, 2026.", "url": "https://epoch.ai/frontiermath/open-problems/genus-2-jacobian-torsion", "used_for": "identification of the [31] curve"},
    "AlessandriCoppola": {"short": "Alessandrì–Coppola 2026", "cite": "J. Alessandrì and N. Coppola, Torsion points on GL₂-type abelian varieties, arXiv:2602.21047 (2026).", "url": "https://arxiv.org/abs/2602.21047", "used_for": "the 31-torsion candidate 1830.2.a.q"},
    "DaowsudSchmidt": {"short": "Daowsud–Schmidt 2018", "cite": "K. Daowsud and T. A. Schmidt, Continued fractions for rational torsion, J. Number Theory 189 (2018), 115–130; corrigendum J. Number Theory 246 (2023), 326–327.", "url": "https://doi.org/10.1016/j.jnt.2017.11.001", "used_for": "the order-11 family (and Flynn's older order-11 family quoted there)"},
    "KuruSadek": {"short": "Kuru–Sadek 2024", "cite": "H. Kuru and M. Sadek, Quadratic torsion orders on Jacobian varieties, arXiv:2410.14455 (2024).", "url": "https://arxiv.org/abs/2410.14455", "used_for": "the order-23 family (genus-2 specialization)"},
    "Choudhry": {"short": "Choudhry 2016", "cite": "A. Choudhry, Equal sums of like powers with minimum number of terms, Integers 16 (2016), Paper No. A77.", "url": "", "used_for": "the genus-1 curve on the [2,2,2,12] surface (Theorem 3.3 of BNSS2026)"},
    "Stoll1999": {"short": "Stoll 1999", "cite": "M. Stoll, On the height constant for curves of genus two, Acta Arith. 90 (1999), 183–201.", "url": "https://doi.org/10.4064/aa-90-2-183-201", "used_for": "the algorithm behind Magma's TorsionSubgroup for genus-2 Jacobians over ℚ"},
    "MullerStoll2016": {"short": "Müller–Stoll 2016", "cite": "J. S. Müller and M. Stoll, Canonical heights on genus-2 Jacobians, Algebra Number Theory 10 (2016), 2153–2234; errata, published online 2023.", "url": "https://doi.org/10.2140/ant.2016.10.2153", "used_for": "the height bounds used by Magma's TorsionSubgroup"},
    "Zywina2022": {"short": "Zywina 2022", "cite": "D. Zywina, Determining monodromy groups of abelian varieties, Res. Number Theory 8 (2022), article 89.", "url": "https://doi.org/10.1007/s40993-022-00389-x", "used_for": "background for the Frobenius criteria (End = ℤ certificates in the Genus2Torsion scripts; not used here)"},
    "Lombardo2019": {"short": "Lombardo 2019", "cite": "D. Lombardo, Computing the geometric endomorphism ring of a genus-2 Jacobian, Math. Comp. 88 (2019), 889–929.", "url": "https://doi.org/10.1090/mcom/3358", "used_for": "the endomorphism-algebra criterion cited by the paper for geometric simplicity"},
    "Magma": {"short": "Magma", "cite": "W. Bosma, J. J. Cannon, C. Fieker, A. Steel (eds.), Handbook of Magma functions, Edition 2.29 (2026).", "url": "http://magma.maths.usyd.edu.au/magma/", "used_for": "all verifications (TorsionSubgroup, RichelotIsogenousSurfaces, AutomorphismGroup, EulerFactor, Conductor)"},
}

# citation keys as they appear in the Source(s) column of table1.txt / table2.txt -> SOURCES keys
CITE_KEYS = {k: k for k in SOURCES}
CITE_KEYS["new"] = "BNSS2026"


def parse_source_column(s: str) -> tuple[list[str], bool]:
    """'\\cite{BookerSutherland}, new' -> (['BookerSutherland', 'BNSS2026'], rm_flag)."""
    import re
    rm = "RM" in s
    keys = []
    for m in re.finditer(r"\\cite\{([^}]*)\}", s):
        for k in m.group(1).split(","):
            k = k.strip()
            if k:
                keys.append(CITE_KEYS.get(k, k))
    rest = re.sub(r"\\cite\{[^}]*\}", "", s)
    rest = re.sub(r"\$\^\{?\\?(text|rm)?\{?RM\}?\}?\$", "", rest)
    if "new" in rest:
        keys.append("BNSS2026")
    out = []
    for k in keys:
        if k not in out:
            out.append(k)
    return out, rm


# --------------------------------------------------------------------------- infinitely many?

RECORD = "Najman2026"

# simple side: certified 'exact' families (census_infinity_record.md, sections 2 and 3)
_SIMPLE_EXACT = {
    "1": ("y² = x⁵ + x + t (t = 3, 4, 5, 7, 8: trivial torsion, strict primes 7, 7, 7, 13, 13, distinct G2-invariants)", ["ThisSite"]),
    "2": ("y² = x(x⁵ + x + t); exact and strict fiber t = 2 (p = 17)", [RECORD]),
    "2.2": ("y² = x(x−1)(x−2)(x³ + x + t); t = 3 (p = 11)", [RECORD]),
    "2.2.2": ("y² = x(x−1)(x−2)(x² + x + t); t = 3 (p = 29)", [RECORD]),
    "2.2.2.2": ("y² = x(x−1)(x−2)(x−3)(x−t); t = 5 (p = 11)", [RECORD]),
    "3": ("y² = (x³+x+1)² − t(x²+1)³ (contact form f = h² − λQ³); t = 2 (p = 13)", [RECORD]),
    "10": ("contact-5 chart at b = 1: y² = (1 + ax + x²)² − (a+2)²x⁵; eight exact fibers", [RECORD]),
    "8": ("the A(8) chart, slice p = 2, r = 3 with free parameter t; t = 2 (p = 13)", [RECORD]),
    "4": ("y² = x·Q₄(x), Q₄ the norm from ℚ(√2) of the automatic-halving quadratic q_{α,β}, α = 1 + s√2, β = 1; eight exact fibers", [RECORD]),
    "2.4": ("y² = x(x² + (2t+1)x + t²)(x² + 3x + 1); eight exact fibers", [RECORD]),
    "2.2.4": ("y² = x(x+4)(x+9)(x² + (2t+1)x + t²); seven exact fibers", [RECORD]),
    "2.2.2.4": ("y² = x(x+1)(x+4)(x+9)(x+t²) (the square-branch chart of BNSS2026, Lemma 2.3); t = 5 (p = 59)", [RECORD, "BNSS2026"]),
    "2.2.4.4": ("the rational curve in A(2,2,4,4) of BNSS2026, Proposition 4.4; s = 3 (p = 37)", [RECORD, "BNSS2026"]),
    "2.2.2.8": ("the rational curve on the K3 surface M(2,2,2,8) of BNSS2026, Theorem 4.7: (a,b,c,d) = (4t²(t+1)/(t²+t+1)², t/(t+1), −1, −t) on y² = x(x+a²)(x+b²)(x+c²)(x+d²); t = 2 (p = 31)", [RECORD, "BNSS2026"]),
    "2.2.2.6": ("the M(2,2,2,6) chart: y² = x(x+2s²−sn)(x+2s²+sm−2sn−mn)(x+2s²+sm−sn−mn)(2x−mn)(2x+4s²−4sn−mn); five exact fibers", [RECORD]),
    "4.4": ("Richelot partner #5 of the [2,2,4,4] family (fibers s = 3, 5, 7, 9, 11)", [RECORD]),
    "2.8": ("Richelot partner #4 of the [2,2,2,8] family (fibers t = 2, 3, 5, 7)", [RECORD]),
    "2.10": ("Elkies' 5-torsion chart with three rational Weierstrass points, slice (t, 2, 3); five exact fibers", [RECORD, "Elkies2024"]),
    "2.2.10": ("Elkies' 5-torsion chart with four rational Weierstrass points on the nodal cubic a³ + 5a² − (5/9)ab² + (1/3)b² = 0 (parametrised by τ); five exact fibers", [RECORD, "Elkies2024"]),
    "7": ("contact-7 chart: h = 1 − (7/2)x + ax² + bx³, f = (h² + (x−1)⁷)/x²", [RECORD]),
    "9": ("contact-9 chart: h = 1 − (9/2)x + (63/8)x² − (105/16)x³ + ax⁴, f = (h² + (x−1)⁹)/x⁴", [RECORD]),
    "11": ("Flynn's family y² = x⁶ + 2x⁵ + (2t+3)x⁴ + 2x³ + (t²+1)x² + 2t(1−t)x + t² and the Daowsud–Schmidt continued-fraction family", [RECORD, "DaowsudSchmidt"]),
    "6": ("the plain M(12) chart T₁₂ = ax² − x + r, h = (x−r)(T₁₂+1), f = ax²T₁₂(T₁₂+1)", [RECORD]),
    "12": ("the rational-Weierstrass subchart a = (1−z²)/(4(r+1)) of the M(12) chart", [RECORD]),
    "2.12": ("the line a = (1−r)/4 in the M(12) chart", [RECORD]),
    "2.6": ("contact-6 chart: f = h² − (x−1)⁶, h cubic", [RECORD]),
    "14": ("rational-root subfamily of the contact-7 chart", [RECORD]),
    "18": ("rational-root subfamily of the contact-9 chart", [RECORD]),
    "5": ("five-torsion families (see the record)", [RECORD]),
    "15": ("inherited from the cyclic [30] family (see the record for the exactness discussion)", [RECORD]),
    "20": ("the contact-5 + 4-torsion family", [RECORD]),
    "2.20": ("the extra-2 loci of the contact-5 + 4-torsion family", [RECORD]),
    "21": ("Leprévost's 1991 one-parameter degree-5 family f₂₁ = A₂₁² − k₂₁x³(x−1)²; marked class of order 21 verified symbolically over ℚ(t)", [RECORD, "Leprevost1991"]),
    "22": ("rational-branch subfamilies of the order-11 families", [RECORD]),
    "23": ("the Kuru–Sadek quadratic-order construction, genus-2 specialization", [RECORD, "KuruSadek"]),
    "30": ("the simultaneous contact-5 / contact-6 family", [RECORD]),
    "32": ("the reconstructed Elkies order-32 component", [RECORD, "Elkies2002"]),
    "3.3": ("Bruin–Flynn–Testa 2014, Theorem 6: rational parametrization of the pointwise A(3,3) chart", [RECORD, "BFT2014"]),
    "4.8": ("the tangent-cover family", [RECORD]),
    "2.2.2.10": ("Elkies 2024: the Clebsch–Klein full-level-2 plus 5-torsion two-parameter family", [RECORD, "Elkies2024"]),
}

# simple side: a family with J(Q)_tors ⊇ G is proven, but not exactness/simplicity of infinitely many members
_SIMPLE_FAMILY = {
    "13": ("Leprévost 1991: a family of genus-2 curves with a rational divisor class of order 13 (formulas not transcribed into the record; exactness and geometric simplicity of infinitely many members not certified)", ["Leprevost1991a"]),
    "17": ("Leprévost 1991: a family with a rational divisor class of order 17 (not transcribed; open legs as for [13])", ["Leprevost1991"]),
    "19": ("Leprévost 1991: a family with a rational divisor class of order 19 (not transcribed; open legs as for [13])", ["Leprevost1991"]),
    "2.2.2.12": ("BNSS2026, Theorem 3.3: infinitely many genus-2 Jacobians over ℚ with J(ℚ)tors ⊇ [2,2,2,12], from a genus-1 curve of rank 1 on the surface a²+b²+c² = u²+v², a⁴+b⁴+c⁴ = u⁴+v⁴; whether infinitely many of them have torsion exactly [2,2,2,12] and are geometrically simple is open", ["BNSS2026", "Choudhry"]),
}

# simple side: open, with a note from the record
_SIMPLE_OPEN_NOTES = {
    "16": "the record traces no positive-dimensional construction with generic torsion exactly [16] (the reconstructed Elkies [32] component contains [16] as 2D, but its generic torsion is [32])",
    "26": "the record mentions Platonov–Petrunin (2012) constructions, not transcribed or certified",
    "28": "the record mentions Platonov–Petrunin (2012) constructions, not transcribed or certified",
    "3.9": "only isolated examples; (2,2)-isogenies cannot create 3-power structure",
    "6.6": "the contact-6 × contact-6 locus is positive-dimensional with one certified point, but no rational curve on it is known",
    "2.2.8": "no positive-dimensional construction found (a second-level halving criterion was designed but not derived)",
    "2.4.4": "no positive-dimensional construction found; provably not a Richelot image of the [2,2,4,4] family",
    "2.4.8": "one verified exact hit on the one-split tangent subcover; the genus/rank of its base curve was not computed",
    "31": "every known example has real multiplication (by ℤ[√2]); a single curve is known",
    "2.22": "every known example has real multiplication (by ℤ[φ], φ² = φ + 1); two curves are known",
    "2.2.20": "the displayed curve is the unique one in its one-parameter family (BNSS2026, §'A second order-80 group', by Chabauty); Elkies' family gives infinitely many Jacobians of order 80, but with group [2,2,2,10]",
}

# split side: HLP 2000, Theorem 1 / Table 1 (families with J(Q)_tors ⊇ G; the Jacobians are
# isogenous to products of elliptic curves by construction)
_SPLIT_FAMILY = {
    "20": "ℙ²", "21": "ℙ²", "3.9": "ℙ²", "30": "ℙ²", "35": "positive rank elliptic curve",
    "6.6": "ℙ²", "3.12": "ℙ²", "40": "positive rank elliptic surface", "45": "positive rank elliptic curve",
    "2.24": "ℙ²", "5.10": "positive rank elliptic surface", "60": "positive rank elliptic curve",
    "8.8": "ℙ²", "2.4.8": "ℙ²", "6.12": "positive rank elliptic surface", "2.6.6": "positive rank elliptic surface",
    "2.2.24": "positive rank elliptic curve", "2.2.4.8": "positive rank elliptic surface",
}
_SPLIT_OPEN_NOTES = {
    "7.7": "HLP 2000, Table 1: a single curve (parameter space ℙ⁰)",
    "63": "HLP 2000, Table 1: a single curve (parameter space ℙ⁰)",
}


def infinite_record(side: str, key: str) -> dict:
    """The ∞ answer for a group key on a side ('simple' or 'split')."""
    if side == "simple":
        if key in _SIMPLE_EXACT:
            fam, src = _SIMPLE_EXACT[key]
            return {"grade": "exact", "sources": src, "note": fam}
        if key in _SIMPLE_FAMILY:
            fam, src = _SIMPLE_FAMILY[key]
            return {"grade": "family", "sources": src, "note": fam}
        if key in _SIMPLE_OPEN_NOTES:
            src = ["BNSS2026", RECORD] if key in ("31", "2.22", "2.2.20") else [RECORD]
            return {"grade": "open", "sources": src, "note": _SIMPLE_OPEN_NOTES[key]}
        return {"grade": "open", "sources": [], "note": ""}
    if side == "split":
        if key in _SPLIT_FAMILY:
            return {"grade": "family", "sources": ["HLP2000"],
                    "note": f"HLP 2000, Theorem 1: a family of genus-2 curves over ℚ with split Jacobians whose torsion contains the group, parametrised by the rational points of a non-empty open subset of a {_SPLIT_FAMILY[key]}"}
        return {"grade": "open", "sources": ["HLP2000"] if key in _SPLIT_OPEN_NOTES else [], "note": _SPLIT_OPEN_NOTES.get(key, "")}
    raise ValueError(side)


GRADE_LABEL = {"exact": "∞", "family": "∞ ⊇", "open": "?"}
GRADE_TEXT = {
    "exact": "infinitely many, with this exact torsion group (certified)",
    "family": "a positive-dimensional family with torsion containing the group is proven; exactness for infinitely many members is open",
    "open": "no positive-dimensional family is recorded",
}
