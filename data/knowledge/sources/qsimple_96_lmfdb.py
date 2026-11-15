# Validation queries for the census proof qsimple-96-obstruction (not used in the proof).
# Run: python3 qsimple_96_lmfdb.py > qsimple_96_lmfdb.log      (needs psycopg2; public LMFDB mirror)
# Tables: g2c_curves (the LMFDB genus-2 curves) and g2c_curves_new (the extended table).
import datetime
import psycopg2

conn = psycopg2.connect(host="devmirror.lmfdb.xyz", port=5432, dbname="lmfdb", user="lmfdb",
                        password="lmfdb", connect_timeout=60)
cur = conn.cursor()


def two_part(n):
    k = 1
    while n % 2 == 0:
        n //= 2
        k *= 2
    return k


ROT = ["E_3", "E_4", "E_6", "J(E_3)", "J(E_4)", "J(E_6)"]
ok = True
print("date", datetime.date.today().isoformat())
for tab in ["g2c_curves", "g2c_curves_new"]:
    cur.execute("select count(*) from %s" % tab)
    print("TABLE", tab, "rows", cur.fetchone()[0])

    # (1) geometrically split, simple over Q, with a rational point of order 3
    cur.execute("select torsion_order, torsion_subgroup, count(*) from %s "
                "where is_simple_geom = false and is_simple_base = true and torsion_order %% 3 = 0 "
                "group by torsion_order, torsion_subgroup order by torsion_order, torsion_subgroup" % tab)
    rows = cur.fetchall()
    n = sum(r[2] for r in rows)
    hist = {}
    for o, t, c in rows:
        hist[two_part(o)] = hist.get(two_part(o), 0) + c
    m = max(hist) if hist else 0
    print("  (1) geometrically split, simple over Q, 3 | #torsion: %d curves; 2-part of the torsion order -> number of curves: %s"
          % (n, sorted(hist.items())))
    print("      torsion subgroups with 2-part >= 8:", [(t, c) for o, t, c in rows if two_part(o) >= 8])
    print("      largest 2-part: %d (the Theorem predicts <= 16)" % m)
    ok = ok and m <= 16

    # (2) type (S) with an element of order 3, 4 or 6 in Gamma (Sato-Tate groups E_n, J(E_n), n = 3, 4, 6)
    for st in ROT:
        cur.execute("select torsion_order, count(*) from %s where st_group = %%s group by torsion_order" % tab, (st,))
        rows = cur.fetchall()
        n = sum(r[1] for r in rows)
        m2 = max([two_part(r[0]) for r in rows] or [0])
        mo = max([r[0] for r in rows] or [0])
        print("  (2) st_group %-7s %7d curves; largest torsion order %d; largest 2-part %d (Proposition 8 predicts <= 16)"
              % (st, n, mo, m2))
        ok = ok and m2 <= 16

    # (3) geometrically split, simple over Q, with at least 32 rational points of 2-power order, by Sato-Tate group
    cur.execute("select st_group, torsion_order, count(*) from %s "
                "where is_simple_geom = false and is_simple_base = true and torsion_order %% 32 = 0 "
                "group by st_group, torsion_order order by st_group, torsion_order" % tab)
    rows = cur.fetchall()
    print("  (3) geometrically split, simple over Q, 32 | #torsion: (st_group, torsion order, curves):", rows)
    bad = [r for r in rows if r[0] in ROT]
    ok = ok and not bad

    # (4) all curves with 96 | #torsion, by simplicity
    cur.execute("select is_simple_geom, is_simple_base, torsion_subgroup, count(*) from %s "
                "where torsion_order %% 96 = 0 group by is_simple_geom, is_simple_base, torsion_subgroup" % tab)
    rows = cur.fetchall()
    print("  (4) curves with 96 | #torsion (is_simple_geom, is_simple_base, torsion, curves):", rows)
    ok = ok and all(r[0] or not r[1] for r in rows)
print("LMFDB_CHECK_OK" if ok else "LMFDB_CHECK_FAILED")
