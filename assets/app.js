/* genus-2-torsion: static front end.  Reads data/groups.json, data/curves.json and data/sources.json
   (built by pipeline/build.py on Mordell) and renders the census.  No framework, no build step. */
(function () {
  "use strict";

  const REPO = "https://github.com/Genus-2-torsion/Genus-2-torsion.github.io";
  const ISSUE_FORM = REPO + "/issues/new?template=submit-curve.yml";
  const pageName = document.body.dataset.page || "";
  const CLASSES = [["simple", "geometrically simple"], ["qsplit", "split over ℚ"], ["gsplit", "geometrically split, simple over ℚ"]];
  const CLASS_TEXT = Object.fromEntries(CLASSES);
  const GRADE = { exact: "∞", family: "∞ ⊇", open: "?" };
  const GRADE_TEXT = {
    exact: "infinitely many, with this exact torsion group (certified)",
    family: "a positive-dimensional family with torsion containing the group is proven; exactness for infinitely many members is open",
    open: "no positive-dimensional family is recorded",
  };

  // ---------------------------------------------------------------- utilities
  const $ = (sel, root) => (root || document).querySelector(sel);
  const el = (tag, attrs, ...children) => {
    const node = document.createElement(tag);
    if (attrs) for (const [k, v] of Object.entries(attrs)) {
      if (k === "class") node.className = v;
      else if (k === "html") node.innerHTML = v;
      else if (k.startsWith("on")) node.addEventListener(k.slice(2), v);
      else if (v !== null && v !== undefined && v !== false) node.setAttribute(k, v);
    }
    for (const c of children.flat(Infinity)) {
      if (c === null || c === undefined || c === false) continue;
      node.append(c.nodeType ? c : document.createTextNode(String(c)));
    }
    return node;
  };
  const esc = (s) => String(s).replace(/[&<>"]/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;" }[c]));
  const params = new URLSearchParams(location.search);
  const cache = {};
  async function load(name) {
    if (!cache[name]) {
      cache[name] = fetch("data/" + name + ".json", { cache: "no-cache" }).then((r) => {
        if (!r.ok) throw new Error("could not load data/" + name + ".json (" + r.status + ")");
        return r.json();
      });
    }
    return cache[name];
  }

  // Polynomials as Magma prints them -> readable HTML (x^3 - 3*x^2 + 3 -> x³ − 3x² + 3).
  const SUP = { "0": "⁰", "1": "¹", "2": "²", "3": "³", "4": "⁴", "5": "⁵", "6": "⁶", "7": "⁷", "8": "⁸", "9": "⁹", "-": "⁻" };
  function math(s) {
    if (s === null || s === undefined || s === "") return "";
    let t = esc(String(s));
    t = t.replace(/\^(-?\d+)/g, (_, d) => d.split("").map((c) => SUP[c]).join(""));
    t = t.replace(/(\d|\)|[a-zA-Z])\*([a-zA-Zζ(])/g, "$1$2");   // 3*x -> 3x, (..)*y -> (..)y
    t = t.replace(/\*/g, "·");
    t = t.replace(/(^|[\s(=])-/g, "$1−").replace(/ - /g, " − ");
    t = t.replace(/\bQ\b/g, "ℚ").replace(/Qbar/g, "ℚ̄").replace(/sqrt\((-?\d+)\)/g, "√$1");
    t = t.replace(/[−-]&gt;/g, "→").replace(/E x E\^sigma/g, "E × E<sup>σ</sup>").replace(/\bchi\b/g, "χ").replace(/\bpi\b/g, "π");
    return '<span class="m">' + t + "</span>";
  }
  const eqHtml = (eq) => math(eq.text);
  // Magma-style text of y^2 + h y = f from the coefficient strings / polynomial strings
  const eqText = (fs, hs) => (hs === "0" ? "y^2 = " + fs : hs === "1" ? "y^2 + y = " + fs
    : (/ [+-] /.test(hs) ? "y^2 + (" + hs + ")*y = " : "y^2 + " + hs + "*y = ") + fs);
  const bracket = (inv) => "[" + inv.join(",") + "]";
  const chip = (status, cls) => {
    const text = status === "certified" ? "certified" : status === "verified" ? "verified · class not certified" : status;
    return el("span", { class: "chip " + (status === "verified" ? "undecided" : status) }, text);
  };
  const infBadge = (rec) => el("span", { class: "inf-badge " + rec.grade + (rec.inherited_from ? " inherited" : ""),
    title: GRADE_TEXT[rec.grade] + (rec.note ? " — " + rec.note : "") }, GRADE[rec.grade]);
  const discovery = (by, year) => (by ? by + (year ? " (" + year + ")" : "") : (year ? String(year) : "—"));
  const srcYear = (src) => (src.short.includes(String(src.year)) ? "" : " (" + src.year + ")");

  function lmfdbLink(l, long) {
    if (!l || !l.url) return null;
    if (l.kind === "production") return el("a", { class: "lmfdb", href: l.url, target: "_blank", rel: "noopener", title: "LMFDB home page of the curve" }, l.label);
    if (l.kind === "alpha-snapshot") return el("a", { class: "lmfdb alpha", href: l.url, target: "_blank", rel: "noopener",
      title: "in the extended (alpha) LMFDB; label " + l.label + " is a 2026-08 snapshot, the link looks the curve up by its equation" }, long ? "alpha LMFDB " + l.label : "alpha " + l.label);
    return el("a", { class: "lmfdb jump", href: l.url, target: "_blank", rel: "noopener", title: "not found in the LMFDB at verification time; this link searches the LMFDB by equation" }, "LMFDB?");
  }

  function prefillIssue(g, cls) {
    const u = new URL(ISSUE_FORM);
    if (g) u.searchParams.set("group", bracket(g));
    if (cls) u.searchParams.set("class", cls);
    u.searchParams.set("title", "Genus 2 curve with torsion " + (g ? bracket(g) : "[...]"));
    return u.toString();
  }

  // ---------------------------------------------------------------- header / footer
  function chrome() {
    const nav = [["index.html", "groups", "index"], ["submit.html", "submit a curve", "submit"],
                 ["about.html", "about & sources", "about"], [REPO, "github", ""]];
    const header = el("header", { class: "site" },
      el("div", { class: "inner" },
        el("h1", null, el("a", { href: "index.html" }, "Torsion of genus 2 Jacobians over ℚ")),
        el("nav", null, nav.map(([href, label, key]) =>
          el("a", { href, "aria-current": key && key === pageName ? "page" : null,
                    target: href.startsWith("http") ? "_blank" : null, rel: href.startsWith("http") ? "noopener" : null }, label)))));
    document.body.prepend(header);
    const footer = el("footer", { class: "site" },
      el("div", { class: "links" },
        el("a", { href: "index.html" }, "groups"), el("span", { class: "sep" }, "·"),
        el("a", { href: "submit.html" }, "submit"), el("span", { class: "sep" }, "·"),
        el("a", { href: "about.html" }, "about & sources"), el("span", { class: "sep" }, "·"),
        el("a", { href: "data/groups.json" }, "all groups (JSON)"), el("span", { class: "sep" }, "·"),
        el("a", { href: "data/curves.json" }, "all curves (JSON)"), el("span", { class: "sep" }, "·"),
        el("a", { href: REPO, target: "_blank", rel: "noopener" }, "source & data on GitHub")),
      el("div", null,
        "Submissions are verified with Magma on the Mordell workstation of the Department of Mathematics, ",
        "University of Zagreb, and stored with their verification logs in the GitHub repository. ",
        "Please cite the original references of the curves you use."));
    document.body.append(footer);
  }

  // ---------------------------------------------------------------- sortable tables
  function makeSortable(table) {
    const ths = table.querySelectorAll("thead th.sortable");
    ths.forEach((th) => th.addEventListener("click", () => {
      const idx = th.cellIndex;
      const asc = !(th.classList.contains("sorted") && th.classList.contains("asc"));
      ths.forEach((t) => t.classList.remove("sorted", "asc"));
      th.classList.add("sorted"); if (asc) th.classList.add("asc");
      const tbody = table.tBodies[0];
      const rows = Array.from(tbody.rows);
      const key = (r) => { const c = r.cells[idx]; const v = c.dataset.sort !== undefined ? c.dataset.sort : c.textContent.trim(); const f = parseFloat(v); return isNaN(f) ? v : f; };
      rows.sort((a, b) => { const x = key(a), y = key(b); return (x < y ? -1 : x > y ? 1 : 0) * (asc ? 1 : -1); });
      rows.forEach((r) => tbody.append(r));
    }));
  }

  // ---------------------------------------------------------------- index page
  const shortSource = (sources, k) => (sources[k] && sources[k].short) || k;
  function sourceLinks(sources, keys) {
    const out = [];
    (keys || []).forEach((k, i) => {
      if (i) out.push(", ");
      const s = sources[k];
      out.push(s ? el("a", { href: "about.html#src-" + k, title: s.cite }, s.short || k) : k);
    });
    return out;
  }

  function exampleCell(list, sources, g, cls, mode) {
    const side = cls === "simple" || cls === "qsplit" || (cls === "gsplit" && !$("#cls-qsplit").checked) ? " side-start" : "";
    if (!list.length) return el("td", { class: "ex none" + side }, "—");
    let c = list[0], dim = false;
    const meta = [];
    if (mode === "first") {
      const i = g.first_dated[cls];
      if (i === null || i === undefined) {
        // no curve of this class carries a discovery date: show the smallest-conductor curve, dimmed,
        // with the source the paper credits for the first realisation of the group
        dim = true;
        const fs = g.first_source[cls], src = fs && sources[fs];
        if (src) meta.push(el("span", { title: "the source credited for the first realisation of this group in this class; the curve shown is the smallest-conductor example, whose own discovery is not dated" },
          "first known: ", el("a", { href: "about.html#src-" + fs, title: src.cite }, src.short), srcYear(src)));
      } else {
        c = list[i];
        meta.push(el("span", null, discovery(c.discovered_by, c.year)));
      }
    }
    const lm = lmfdbLink(c.lmfdb);
    if (lm) meta.push(lm);
    if (c.conductor) meta.push(el("span", { title: c.conductor_note || "conductor" }, "N = " + c.conductor));
    if (mode !== "first") {
      if (c.sources && c.sources.length) meta.push(el("span", null, sourceLinks(sources, c.sources)));
      else if (c.discovered_by) meta.push(el("span", null, discovery(c.discovered_by, c.year)));
    }
    if (c.rm) meta.push(el("span", { title: "the Jacobian has real multiplication" }, "RM"));
    if (c.status !== "certified") meta.push(chip(c.status));
    if (list.length > 1) meta.push(el("a", { href: "group.html?g=" + g.key }, "+" + (list.length - 1) + " more"));
    const metaNodes = [];
    meta.forEach((m, i) => { if (i) metaNodes.push(el("span", { class: "sep" }, "·")); metaNodes.push(m); });
    return el("td", { class: "ex" + side },
      el("p", { class: "eq" + (dim ? " dim" : "") }, el("a", { href: "curve.html?id=" + encodeURIComponent(c.id), html: eqHtml(c.equation) })),
      el("span", { class: "meta" }, metaNodes));
  }

  async function renderIndex() {
    const [data, srcData] = await Promise.all([load("groups"), load("sources")]);
    const sources = srcData.sources;
    $("#stat-groups").textContent = data.n_groups;
    $("#stat-curves").textContent = data.n_curves;
    $("#stat-simple").textContent = data.n_groups_by_class.simple;
    $("#stat-qsplit").textContent = data.n_groups_by_class.qsplit;
    $("#stat-gsplit").textContent = data.n_groups_by_class.gsplit;
    $("#stat-inf").textContent = data.n_infinite.simple.exact;
    $("#stat-generated").textContent = data.generated;
    const tbody = $("#groups tbody");
    const boxes = { simple: $("#cls-simple"), qsplit: $("#cls-qsplit"), gsplit: $("#cls-gsplit") };
    const onlyAll = $("#only-all"), showInf = $("#show-inf");
    const sortSel = $("#sort"), modeSel = $("#example-mode");
    function draw() {
      const mode = modeSel.value;
      const shown = CLASSES.map(([k]) => k).filter((k) => boxes[k].checked);
      const showSimple = shown.includes("simple"), showSplit = shown.includes("qsplit") || shown.includes("gsplit");
      // header
      const thead = $("#groups thead");
      thead.innerHTML = "";
      const sup = el("tr", { class: "super" }, el("th", { colspan: 2 }, ""));
      const cols = el("tr", null,
        el("th", { class: "sortable", title: "J(ℚ)tors as invariant factors" }, "group"),
        el("th", { class: "num sortable" }, "order"));
      if (showSimple) {
        sup.append(el("th", { class: "side-simple", colspan: showInf.checked ? 2 : 1 }, "geometrically simple"));
        cols.append(el("th", { class: "side-start", title: mode === "first" ? "earliest known curve with geometrically simple Jacobian and this torsion" : "smallest-conductor known curve with geometrically simple Jacobian and this torsion" }, mode === "first" ? "first example" : "smallest conductor"));
        if (showInf.checked) cols.append(el("th", { class: "sortable", title: "infinitely many geometrically simple Jacobians with this torsion? ∞ certified exact · ∞ ⊇ family with torsion containing the group · ? open" }, "∞?"));
      }
      if (showSplit) {
        const n = (shown.includes("qsplit") ? 1 : 0) + (shown.includes("gsplit") ? 1 : 0) + (showInf.checked ? 1 : 0);
        sup.append(el("th", { class: "side-split", colspan: n }, "geometrically split"));
        if (shown.includes("qsplit")) cols.append(el("th", { class: "side-start", title: "smallest-conductor known curve whose Jacobian is isogenous over ℚ to a product of elliptic curves" }, "split over ℚ"));
        if (shown.includes("gsplit")) cols.append(el("th", { class: shown.includes("qsplit") ? "" : "side-start", title: "smallest-conductor known curve whose Jacobian is simple over ℚ but splits over ℚ̄" }, "simple over ℚ, split over ℚ̄"));
        if (showInf.checked) cols.append(el("th", { class: "sortable", title: "infinitely many geometrically split Jacobians with this torsion? ∞ ⊇ family with torsion containing the group (Howe–Leprévost–Poonen) · ? open" }, "∞?"));
      }
      thead.append(sup, cols);
      tbody.innerHTML = "";
      let groups = data.groups.slice();
      if (sortSel.value === "order") groups.sort((a, b) => a.order - b.order || a.rank - b.rank);
      for (const g of groups) {
        const has = shown.map((k) => g.known[k]);
        if (!shown.length || (onlyAll.checked ? !has.every(Boolean) : !has.some(Boolean))) continue;
        const undec = g.classes.undecided.length;
        const tr = el("tr", null,
          el("td", { class: "grp", "data-sort": String(g.rank).padStart(2, "0") + g.group.map((n) => String(n).padStart(4, "0")).join("") },
            el("a", { class: "g", href: "group.html?g=" + g.key }, g.bracket),
            el("span", { class: "lbl" }, g.label, undec ? [" · ", el("a", { href: "group.html?g=" + g.key, class: "muted", title: "certified torsion, but the class (or the splitting field) could not be certified" }, undec + " curve" + (undec > 1 ? "s" : "") + " with undecided class")] : null)),
          el("td", { class: "num" }, g.order));
        if (showSimple) {
          tr.append(exampleCell(g.classes.simple, sources, g, "simple", mode));
          if (showInf.checked) tr.append(el("td", { class: "inf", "data-sort": { exact: 0, family: 1, open: 2 }[g.infinite.simple.grade] }, infBadge(g.infinite.simple)));
        }
        if (showSplit) {
          if (shown.includes("qsplit")) tr.append(exampleCell(g.classes.qsplit, sources, g, "qsplit", mode));
          if (shown.includes("gsplit")) tr.append(exampleCell(g.classes.gsplit, sources, g, "gsplit", mode));
          if (showInf.checked) tr.append(el("td", { class: "inf", "data-sort": { exact: 0, family: 1, open: 2 }[g.infinite.split.grade] }, infBadge(g.infinite.split)));
        }
        tbody.append(tr);
      }
      makeSortable($("#groups"));
    }
    Object.values(boxes).forEach((b) => b.addEventListener("change", draw));
    onlyAll.addEventListener("change", draw); showInf.addEventListener("change", draw); sortSel.addEventListener("change", draw); modeSel.addEventListener("change", draw);
    draw();
    loadPending();
  }

  async function loadPending() {
    const box = $("#pending");
    if (!box) return;
    try {
      const r = await fetch("https://api.github.com/repos/Genus-2-torsion/Genus-2-torsion.github.io/issues?labels=submission&state=open&per_page=100");
      if (!r.ok) return;
      const issues = (await r.json()).filter((i) => !i.pull_request);
      if (issues.length) {
        box.innerHTML = "";
        box.append(el("a", { href: REPO + "/issues?q=is%3Aissue+is%3Aopen+label%3Asubmission" },
          issues.length + " submission" + (issues.length > 1 ? "s" : "") + " awaiting verification"));
      }
    } catch (e) { /* offline or rate-limited: say nothing */ }
  }

  // ---------------------------------------------------------------- group page
  function curvesTable(list, sources) {
    const table = el("table", { class: "data" },
      el("thead", null, el("tr", null,
        el("th", { class: "sortable" }, "id"), el("th", null, "curve"),
        el("th", { class: "num sortable", title: "conductor of the Jacobian (blank when the discriminant could not be factored)" }, "N"),
        el("th", null, "LMFDB"), el("th", { class: "sortable" }, "source / discovered by"), el("th", { class: "num sortable", title: "year of discovery, when recorded" }, "year"),
        el("th", null, "certificate"), el("th", null, "status"))),
      el("tbody", null, list.map((c) => el("tr", { class: "row-link", onclick: (e) => { if (e.target.tagName !== "A") location.href = "curve.html?id=" + encodeURIComponent(c.id); } },
        el("td", null, el("a", { class: "id", href: "curve.html?id=" + encodeURIComponent(c.id) }, c.id)),
        el("td", { class: "poly", html: eqHtml(c.equation) }),
        el("td", { class: "num", "data-sort": c.conductor_sort === null ? 1e300 : c.conductor_sort, title: c.conductor_note || "" }, c.conductor || "—"),
        el("td", null, lmfdbLink(c.lmfdb) || el("span", { class: "empty" }, "—")),
        el("td", null, c.sources && c.sources.length ? sourceLinks(sources, c.sources) : (c.discovered_by || "—"), c.rm ? el("span", { class: "muted" }, " · RM") : null, c.historical ? el("span", { class: "muted" }, " · historical example") : null),
        el("td", { class: "num", "data-sort": c.year === null ? 9999 : c.year }, c.year === null ? "" : c.year),
        el("td", { class: "muted", style: "max-width: 22rem; font-size: 0.85rem" }, c.certificate),
        el("td", null, chip(c.status), c.new_group ? [" ", el("span", { class: "chip new", title: "first curve in the census with this group" }, "new group")] : null)))));
    makeSortable(table);
    return el("div", { class: "table-wrap" }, table);
  }

  async function renderGroup() {
    const key = params.get("g");
    const [data, srcData] = await Promise.all([load("groups"), load("sources")]);
    const sources = srcData.sources;
    const g = data.groups.find((x) => x.key === key);
    const title = $("#title");
    if (!g) { title.textContent = "Unknown group"; $("#body").append(el("p", { class: "error" }, "No curve with torsion group " + esc(key) + " is in the census. "), el("a", { class: "btn", href: prefillIssue(key ? key.split(".").map(Number) : null) }, "Submit one")); return; }
    document.title = g.bracket + " — Torsion of genus 2 Jacobians";
    title.innerHTML = math("J(Q)_tors ≅ " + g.bracket).replace("_tors", "<sub>tors</sub>");
    $("#subtitle").innerHTML = (g.rank ? esc(g.label) + ", of order " + g.order : "the trivial group") + ". " +
      g.n_curves + " curve" + (g.n_curves > 1 ? "s" : "") + " in the census.";
    $("#submit-link").href = prefillIssue(g.group);
    const body = $("#body");
    const sideBox = (name, known, inf, extra) => el("div", { class: "panel" },
      el("h3", null, name),
      el("div", { class: "count" }, known),
      el("p", null, "infinitely many? ", infBadge(inf), " ", el("span", null, GRADE_TEXT[inf.grade]),
        inf.note ? ". " + inf.note : null,
        inf.sources && inf.sources.length ? [" (", sourceLinks(sources, inf.sources), ")"] : null),
      extra || null);
    const nS = g.classes.simple.length, nQ = g.classes.qsplit.length, nG = g.classes.gsplit.length;
    const undecSplit = g.classes.undecided.filter((c) => c.class === "split_undecided_over_Q"), undecAll = g.classes.undecided.filter((c) => c.class !== "split_undecided_over_Q");
    const nSp = nQ + nG + undecSplit.length;
    body.append(el("div", { class: "class-cols" },
      sideBox("geometrically simple", nS ? nS + " curve" + (nS > 1 ? "s" : "") : "no example known", g.infinite.simple),
      sideBox("geometrically split", nSp ? [nSp + " curve" + (nSp > 1 ? "s" : ""), el("span", { class: "muted", style: "font-size:0.9rem; font-family: var(--sans)" },
          " (" + nQ + " split over ℚ, " + nG + " simple over ℚ" + (undecSplit.length ? ", " + undecSplit.length + " undecided over ℚ" : "") + ")")] : "no example known", g.infinite.split)));
    if (undecSplit.length) body.append(el("p", { class: "notice" }, undecSplit.length + " curve" + (undecSplit.length > 1 ? "s" : "") + " with this torsion group " + (undecSplit.length > 1 ? "are" : "is") +
      " certified geometrically split, but neither a certificate of splitness over ℚ nor one of simplicity over ℚ was found (see below)."));
    if (undecAll.length) body.append(el("p", { class: "notice" }, undecAll.length + " curve" + (undecAll.length > 1 ? "s" : "") + " with this torsion group " + (undecAll.length > 1 ? "have" : "has") + " a verified torsion subgroup but no certificate of simplicity or splitness (see below)."));
    for (const [cls, name] of CLASSES) {
      if (!g.classes[cls].length) continue;
      const fs = g.first_source[cls], src = fs && sources[fs];
      body.append(el("h3", null, name),
        el("p", null, "Sorted by conductor. ", src ? ["First realisation credited to ", el("a", { href: "about.html#src-" + fs, title: src.cite }, src.short), srcYear(src),
          g.first_dated[cls] === null ? "; the curves below are the smallest-conductor examples, not necessarily the curve of that source." : "."] : null),
        curvesTable(g.classes[cls], sources));
    }
    if (undecSplit.length) body.append(el("h3", null, "geometrically split, undecided over ℚ"), curvesTable(undecSplit, sources));
    if (undecAll.length) body.append(el("h3", null, "class not certified"), curvesTable(undecAll, sources));
  }

  // ---------------------------------------------------------------- curve page
  async function renderCurve() {
    const id = params.get("id");
    const [cv, srcData, groups] = await Promise.all([load("curves"), load("sources"), load("groups")]);
    const sources = srcData.sources;
    const c = cv.curves.find((x) => x.id === id);
    if (!c) { $("#title").textContent = "Unknown curve"; $("#body").append(el("p", { class: "error" }, "No curve with id " + esc(id) + " in the census.")); return; }
    const g = groups.groups.find((x) => x.key === c.group_key) || {};
    document.title = c.id + " — Torsion of genus 2 Jacobians";
    $("#title").textContent = "Curve " + c.id;
    const cu = c.curve;
    const eq = { text: eqText(cu.f_str, cu.h_str) };
    $("#subtitle").innerHTML = "A genus 2 curve over ℚ whose Jacobian has rational torsion subgroup <b>" + esc(c.group_bracket) + "</b>" +
      (c.class_certified ? " and is " + esc(c.class_text) : "") + ".";
    $("#crumb-group").textContent = c.group_bracket;
    $("#crumb-group").href = "group.html?g=" + c.group_key;
    $("#json-link").href = "data/curves/" + c.id + ".json";
    $("#log-link").href = c.verification.log;
    const body = $("#body");

    // status panel
    const cert = c.split.certificate;
    body.append(el("div", { class: "panel " + (c.status === "certified" ? "accepted" : "open") },
      el("h3", null, chip(c.status), " ", c.new_group ? el("span", { class: "chip new" }, "first curve with this group") : c.new_for_class ? el("span", { class: "chip new" }, "first curve with this group in this class") : null),
      el("dl", { class: "facts" },
        el("dt", null, "J(ℚ)tors"), el("dd", null, el("b", null, c.group_bracket), " = " + (g.label || ""), ", order " + c.order,
          " ", el("span", { class: "src" }, "(Magma TorsionSubgroup; #J(ℚ)tors divides gcd #J(𝔽p) = " + c.torsion.point_count_gcd + " over p = " + c.torsion.point_count_primes.join(", ") + ")")),
        el("dt", null, "class"), el("dd", null, el("b", null, c.class_text), c.class_certified ? "" : " — not certified"),
        el("dt", null, "certificate"), el("dd", null,
          c.simplicity.geometrically_simple
            ? ["the Frobenius polynomial at p = " + c.simplicity.strict_prime + " is ", el("span", { html: math(c.simplicity.chi) }), ": irreducible, and no power πⁿ (n ≤ 12) of a root generates a proper subfield, so the Jacobian is geometrically simple"]
            : cert ? [el("span", { html: math(cert.detail) }), c.q_simple.certified ? ["; simple over ℚ because χ at p = " + c.q_simple.prime + " is irreducible: ", el("span", { html: math(c.q_simple.chi) })] : (c.split.split_over_Q ? "" : "; whether it splits over ℚ is undecided (χ is reducible at every good prime tried, and no certificate over ℚ was found)")]
            : ["none found: no strict prime below " + c.simplicity.primes_tried_below + " (" + c.simplicity.admissible_primes + " good primes, χ reducible at " + c.simplicity.reducible_chi + " of them)",
               c.split.signature.all_good_primes_fail_strictness ? "; every good prime below 200 fails the strictness test, the signature of a geometrically split Jacobian (evidence, not a proof)" : "",
               c.q_simple.certified ? ["; simple over ℚ (χ at p = " + c.q_simple.prime + " is irreducible)"] : ""]),
        c.split.cover_submitted ? el("dt", null, "submitted cover") : null,
        c.split.cover_submitted ? el("dd", null, c.split.cover_accepted ? "accepted (the polynomial identity holds)" : "rejected: " + c.split.cover_error) : null,
        el("dt", null, "LMFDB"), el("dd", null, lmfdbLink(c.lmfdb, true) || el("span", { class: "empty" }, "not in the LMFDB"),
          c.lmfdb.kind === "alpha-snapshot" ? el("span", { class: "muted" }, " (extended database; the label is a 2026-08 snapshot, not a permalink)") : c.lmfdb.kind === "jump" ? el("span", { class: "muted" }, " (not found at verification time; the link searches by equation)") : null),
        el("dt", null, "source"), el("dd", null, c.sources && c.sources.length ? sourceLinks(sources, c.sources) : discovery(c.discovery.by, c.discovery.year), c.rm ? " · the Jacobian has real multiplication" : ""),
        c.reference ? el("dt", null, "reference") : null, c.reference ? el("dd", null, c.reference) : null,
        c.notes ? el("dt", null, "notes") : null, c.notes ? el("dd", null, c.notes) : null,
        c.geometrically_isomorphic_to && c.geometrically_isomorphic_to.length ? el("dt", null, "same G2-invariants as") : null,
        c.geometrically_isomorphic_to && c.geometrically_isomorphic_to.length ? el("dd", null, c.geometrically_isomorphic_to.map((x, i) => [i ? ", " : "", el("a", { href: "curve.html?id=" + x }, x)]), " (geometrically isomorphic; a twist)") : null)));

    // the curve
    body.append(el("h3", null, "The curve"), el("p", { class: "eqbig", html: eqHtml(eq) }),
      el("dl", { class: "facts" },
        el("dt", null, "f, h"), el("dd", null, "[" + cu.f.join(", ") + "], [" + cu.h.join(", ") + "] (ascending coefficients, LMFDB convention)"),
        el("dt", null, "even model"), el("dd", { class: "poly", html: "y² = " + math(cu.g_str) + (cu.g_scale !== "1" ? " (4f + h² scaled by " + esc(cu.g_scale) + "²)" : " (= 4f + h²)") }),
        el("dt", null, "minimal model"), el("dd", { class: "poly" }, cu.minimal_model ? el("span", { html: eqHtml({ text: eqText(polyStr(cu.minimal_model.f), polyStr(cu.minimal_model.h)) }) }) : el("span", { class: "empty" }, "not computed")),
        el("dt", null, "conductor"), el("dd", null, cu.conductor ? [el("b", null, cu.conductor), cu.conductor_factored && cu.conductor_factored !== cu.conductor ? " = " + cu.conductor_factored : "", cu.conductor_note ? el("span", { class: "muted" }, " — " + cu.conductor_note) : null] : el("span", { class: "empty" }, cu.conductor_note || "not computed")),
        el("dt", null, "minimal discriminant"), el("dd", null, cu.minimal_discriminant ? [cu.minimal_discriminant, cu.minimal_discriminant_factored ? " = " + cu.minimal_discriminant_factored : ""] : el("span", { class: "empty" }, "not computed")),
        el("dt", null, "discriminant of the even model"), el("dd", null, cu.discriminant_digits > 40 ? cu.discriminant_digits + " digits" : cu.discriminant),
        el("dt", null, "Igusa–Clebsch invariants"), el("dd", null, "[" + cu.igusa_clebsch.join(", ") + "]"),
        el("dt", null, "G2-invariants"), el("dd", null, "[" + cu.g2_invariants.join(", ") + "]"),
        el("dt", null, "rational Weierstrass points"), el("dd", null, cu.rational_weierstrass_points)));

    // torsion
    body.append(el("h3", null, "The torsion subgroup"),
      el("p", { class: "wide" }, "Generators of J(ℚ)tors in Mumford representation ⟨a(x), b(x), d⟩ on the even model y² = 4f + h² (a divisor of degree d with x-coordinates the roots of a and y = b(x)); computed by Magma's TorsionSubgroup (Stoll's algorithm with the Müller–Stoll height bounds)."),
      c.torsion.generators.length ? el("dl", { class: "facts" }, c.torsion.generators.map((gen) => [
        el("dt", null, "order " + gen.order), el("dd", { class: "gen", html: "⟨" + math(gen.mumford[0]) + ", " + math(gen.mumford[1]) + ", " + esc(gen.mumford[2]) + "⟩" })])) : el("p", null, "The torsion subgroup is trivial."));

    // Magma snippet
    const magma = [
      "R<x> := PolynomialRing(Rationals());",
      "C := HyperellipticCurve(R![" + cu.f.join(", ") + "], R![" + cu.h.join(", ") + "]);",
      "J := Jacobian(SimplifiedModel(C));",
      "Invariants(TorsionSubgroup(J));   // " + c.group_bracket,
      c.simplicity.geometrically_simple ? "// geometric simplicity: EulerFactor(Jacobian(ChangeRing(SimplifiedModel(C), GF(" + c.simplicity.strict_prime + ")))) is irreducible with no degree drop of pi^n, n <= 12" : null,
    ].filter(Boolean).join("\n");
    body.append(el("h3", null, "Reproduce in Magma"), el("pre", null, el("code", null, magma)));

    // submission & verification
    const s = c.submitted, src = c.source || {};
    body.append(el("h3", null, "Submission and verification"),
      el("dl", { class: "facts" },
        el("dt", null, "submitted by"), el("dd", null, [c.submitter.name, c.submitter.affiliation, c.submitter.github ? "@" + c.submitter.github : ""].filter(Boolean).join(", ") || "—"),
        el("dt", null, "submitted"), el("dd", null, c.dates.submitted || "—", s.claimed_group ? ", claimed torsion " + bracket(s.claimed_group) : "", s.claimed_class ? ", claimed class: " + s.claimed_class : ""),
        el("dt", null, "source"), el("dd", null, src.kind === "issue" ? el("a", { href: REPO + "/issues/" + src.number }, "issue #" + src.number)
          : src.kind === "import" ? [src.file + ", line " + src.line + " ", el("a", { href: src.url }, "(Genus2Torsion repository)")] : (src.path || "file")),
        el("dt", null, "verified"), el("dd", null, c.dates.verified + " on " + c.verification.host),
        el("dt", null, "software"), el("dd", null, "Magma " + c.verification.magma_version + ", " + c.verification.cputime_seconds + " s CPU" +
          (c.verification.conductor_cputime_seconds !== null && c.verification.conductor_cputime_seconds !== undefined ? " + " + c.verification.conductor_cputime_seconds + " s for the conductor" : "") +
          "; verify_lib.m sha256 " + c.verification.verify_lib_sha256.slice(0, 12) + "…"),
        el("dt", null, "logs"), el("dd", null, el("a", { href: c.verification.log }, "verification log"), c.verification.conductor_log ? [", ", el("a", { href: c.verification.conductor_log }, "conductor log")] : null),
        el("dt", null, "checks"), el("dd", null, "the model defines a smooth curve of genus 2; J(ℚ)tors was computed exactly and its order divides the gcd of #J(𝔽p) over good primes; the class certificate above was recomputed from scratch" +
          (c.split.cover_submitted ? "; the submitted map to a genus-1 curve was checked by the polynomial identity (e₃p³ + e₂p²q + e₁pq² + e₀q³)·q = g·hh² with nonzero Wronskian" : "") + ".")));
  }
  function polyStr(coeffs) {
    // ascending coefficient strings -> Magma-style string
    const terms = [];
    for (let i = coeffs.length - 1; i >= 0; i--) {
      const c = coeffs[i];
      if (c === "0") continue;
      const neg = c.startsWith("-"); const mag = neg ? c.slice(1) : c;
      const xp = i === 0 ? "" : i === 1 ? "x" : "x^" + i;
      const body = i === 0 ? mag : (mag === "1" ? xp : mag + "*" + xp);
      terms.push([neg, body]);
    }
    if (!terms.length) return "0";
    let s = (terms[0][0] ? "-" : "") + terms[0][1];
    for (const [neg, body] of terms.slice(1)) s += (neg ? " - " : " + ") + body;
    return s;
  }

  // ---------------------------------------------------------------- about page: sources
  async function renderAbout() {
    const srcData = await load("sources");
    const ul = $("#sources");
    for (const [k, s] of Object.entries(srcData.sources)) {
      ul.append(el("li", { id: "src-" + k }, el("b", null, s.short || k), " — ", s.url ? el("a", { href: s.url, target: "_blank", rel: "noopener" }, s.cite) : s.cite,
        s.used_for ? el("span", { class: "muted" }, " Used for: " + s.used_for + ".") : null));
    }
    const groups = await load("groups");
    $("#generated").textContent = groups.generated;
  }

  // ---------------------------------------------------------------- boot
  chrome();
  const sl = $("#submit-any");
  if (sl) sl.href = prefillIssue();
  const run = { index: renderIndex, group: renderGroup, curve: renderCurve, about: renderAbout }[pageName];
  if (run) run().catch((e) => {
    const main = $("main");
    main.append(el("p", { class: "error" }, "Could not load the census data: " + e.message));
  });
})();
