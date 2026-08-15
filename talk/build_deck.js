// Builds the 12-slide JSDP talk deck from talk/slide_blueprint.md.
// Palette and motif mirror poster_final.pdf so the deck and the board read as one project.
const pptxgen = require("pptxgenjs");
const path = require("path");

const A = (f) => path.join(__dirname, "assets", f);

// ---- palette (from poster_template.html :root) --------------------------
const NAVY = "10243C";
const TEAL = "0B6F6A";
const EOSIN = "B8465E";
const HEMA = "3D2F6B";
const HEMA_LT = "6B5AA8";
const INK = "16202C";
const INK2 = "4A5768";
const MUTED = "7A8798";
const RULE = "D5DEE7";
const WASH = "F2F5F9";
const HEMA_WASH = "F0EDF9";
const EOSIN_WASH = "FDF0F3";
const WHITE = "FFFFFF";

const HEAD = "Cambria";
const BODY = "Calibri";

const W = 13.333, H = 7.5;
const ML = 0.6, MR = 0.6;
const CW = W - ML - MR;

const pres = new pptxgen();
pres.layout = "LAYOUT_WIDE";
pres.author = "Kaisar Dauyey";
pres.title = "Foundation-model perturbation analysis in SCLC";

// ---- helpers ------------------------------------------------------------
function title(slide, text, opts = {}) {
  slide.addText(text, {
    x: ML, y: 0.36, w: CW, h: opts.h || 1.15,
    fontFace: HEAD, fontSize: opts.size || 27, bold: true,
    color: opts.color || NAVY, align: "left", valign: "top", margin: 0,
  });
}

function kicker(slide, text, color) {
  slide.addText(text.toUpperCase(), {
    x: ML, y: 0.10, w: CW, h: 0.22,
    fontFace: BODY, fontSize: 11, bold: true, charSpacing: 1.6,
    color: color || TEAL, margin: 0, valign: "middle",
  });
}

function card(slide, x, y, w, h, fill) {
  slide.addShape(pres.ShapeType.roundRect, {
    x, y, w, h, rectRadius: 0.06,
    fill: { color: fill || WASH }, line: { color: RULE, width: 0.75 },
  });
}

function bullets(slide, items, x, y, w, h, opts = {}) {
  slide.addText(
    items.map((t, i) => ({
      text: t,
      options: { bullet: { indent: 14 }, breakLine: i !== items.length - 1 },
    })),
    {
      x, y, w, h, fontFace: BODY, fontSize: opts.size || 14,
      color: opts.color || INK, lineSpacing: opts.lineSpacing || 20,
      paraSpaceAfter: opts.after === undefined ? 7 : opts.after, margin: 0, valign: "top",
    }
  );
}

// numbered circle + label, the deck's repeated motif
function chip(slide, x, y, n, color) {
  slide.addShape(pres.ShapeType.ellipse, {
    x, y, w: 0.33, h: 0.33, fill: { color: color || TEAL }, line: { color: color || TEAL, width: 0 },
  });
  slide.addText(String(n), {
    x, y, w: 0.33, h: 0.33, fontFace: BODY, fontSize: 13, bold: true,
    color: WHITE, align: "center", valign: "middle", margin: 0,
  });
}

function stat(slide, x, y, w, value, label, color) {
  card(slide, x, y, w, 0.95, WHITE);
  slide.addText(value, {
    x: x + 0.04, y: y + 0.08, w: w - 0.08, h: 0.5,
    fontFace: HEAD, fontSize: 25, bold: true, color: color || NAVY,
    align: "center", valign: "middle", margin: 0,
  });
  slide.addText(label, {
    x: x + 0.04, y: y + 0.58, w: w - 0.08, h: 0.3,
    fontFace: BODY, fontSize: 10.5, color: INK2, align: "center", valign: "middle", margin: 0,
  });
}

function footer(slide, n, tag) {
  slide.addText(tag, {
    x: ML, y: 7.02, w: 8, h: 0.28, fontFace: BODY, fontSize: 9.5,
    color: MUTED, margin: 0, valign: "middle",
  });
  slide.addText(String(n), {
    x: W - MR - 0.7, y: 7.02, w: 0.7, h: 0.28, fontFace: BODY, fontSize: 9.5,
    color: MUTED, align: "right", margin: 0, valign: "middle",
  });
}

// caption under a figure
function cap(slide, text, x, y, w) {
  slide.addText(text, {
    x, y, w, h: 0.26, fontFace: BODY, fontSize: 9.5, italic: true,
    color: MUTED, margin: 0, valign: "top",
  });
}

/* ======================= SLIDE 1 — title ============================== */
{
  const s = pres.addSlide();
  s.background = { color: NAVY };

  s.addText("24JSDP  ·  ABSTRACT P25", {
    x: ML, y: 0.55, w: CW, h: 0.3, fontFace: BODY, fontSize: 12, bold: true,
    charSpacing: 2, color: "8FB3C9", margin: 0,
  });

  s.addText(
    "A foundation model points to antigen presentation in SCLC T-cell dysfunction — and independent tissue agrees",
    { x: ML, y: 1.0, w: CW - 0.2, h: 1.8, fontFace: HEAD, fontSize: 30, bold: true, color: WHITE, margin: 0, valign: "top" }
  );

  s.addText(
    [
      { text: "Kaisar Dauyey", options: { bold: true, color: WHITE } },
      { text: "  ·  ", options: { color: "6E8AA6" } },
      { text: "Shinji Nakaoka", options: { bold: true, color: WHITE } },
    ],
    { x: ML, y: 2.92, w: CW, h: 0.3, fontFace: BODY, fontSize: 15, margin: 0 }
  );
  s.addText("Laboratory of Mathematical Biology, Faculty of Advanced Life Science, Hokkaido University", {
    x: ML, y: 3.22, w: CW, h: 0.3, fontFace: BODY, fontSize: 12, color: "9DB4C8", margin: 0,
  });

  // three-node take-home flow
  const ny = 4.15, nh = 1.18, gap = 0.42;
  const nw = (CW - 2 * gap) / 3;
  const nodes = [
    ["Model-derived dysfunction\ncandidates in SCLC", TEAL],
    ["Antigen-presentation\nprogramme", HEMA_LT],
    ["Associated with T-cell abundance\nin independent tumour tissue", EOSIN],
  ];
  nodes.forEach(([txt, col], i) => {
    const x = ML + i * (nw + gap);
    s.addShape(pres.ShapeType.roundRect, {
      x, y: ny, w: nw, h: nh, rectRadius: 0.08,
      fill: { color: "1B3454" }, line: { color: col, width: 1.5 },
    });
    s.addText(txt, {
      x: x + 0.12, y: ny, w: nw - 0.24, h: nh, fontFace: BODY, fontSize: 12.5,
      color: WHITE, align: "center", valign: "middle", margin: 0,
    });
    if (i < 2) {
      s.addText("→", {
        x: x + nw, y: ny, w: gap, h: nh, fontFace: BODY, fontSize: 20,
        color: "6E8AA6", align: "center", valign: "middle", margin: 0,
      });
    }
  });

  s.addText("Hypothesis generator, tissue-validated — not a knockout screen", {
    x: ML, y: 5.62, w: CW - 1.45, h: 0.34, fontFace: BODY, fontSize: 13.5, bold: true, italic: true,
    color: "F0C9D3", margin: 0,
  });

  s.addImage({ path: A("qr.png"), x: W - MR - 1.18, y: 5.62, w: 1.18, h: 1.18 });
  s.addText("github.com/Kays3/geneformer-lung-tcell", {
    x: ML, y: 6.42, w: 7.5, h: 0.3, fontFace: BODY, fontSize: 11, color: "8FB3C9", margin: 0,
  });

  s.addNotes(
    "Open with the clinical stake: SCLC is the most aggressive thoracic malignancy and barely responds to checkpoint blockade, despite a high mutational burden that should make it visible to the immune system.\n\n" +
    "State the result up front: a foundation model, perturbed in silico, nominates dysfunction candidates in SCLC T cells; the programme it points to is associated with T-cell abundance in independent tumour tissue.\n\n" +
    "Name the honest framing now — this is a hypothesis generator validated against tissue, not a knockout experiment. Saying it early buys credibility for everything that follows."
  );
}

/* ====== SLIDE 2 — what an in silico perturbation is (poster panel) ===== */
{
  const s = pres.addSlide();
  kicker(s, "Method");
  title(s, "What is an in silico perturbation?", { h: 0.72 });
  footer(s, 2, "Reading a cell as a ranked gene list");

  s.addText(
    [
      { text: "Geneformer reads a cell as an ", options: { color: INK } },
      { text: "ordered list of its genes", options: { bold: true, color: INK } },
      { text: ", most-expressed first. Deleting or inserting one gene in that list and re-reading the cell asks a counterfactual question: ", options: { color: INK } },
      { text: "if this gene changed, would the cell look more like a different disease state?", options: { italic: true, color: TEAL } },
    ],
    { x: ML, y: 1.20, w: CW, h: 0.72, fontFace: BODY, fontSize: 14.5, lineSpacing: 21, margin: 0, valign: "top" }
  );

  const gap = 0.28;
  const cw = (CW - 3 * gap) / 4;
  const colX = (i) => ML + i * (cw + gap);
  const heads = ["One T cell", "Rank its genes", "Edit one gene", "Measure the move"];
  heads.forEach((h, i) => {
    chip(s, colX(i), 2.32, i + 1, i === 2 ? EOSIN : (i === 3 ? "1A6B3C" : TEAL));
    s.addText(h, {
      x: colX(i) + 0.42, y: 2.32, w: cw - 0.42, h: 0.33,
      fontFace: HEAD, fontSize: 15, bold: true, color: NAVY, margin: 0, valign: "middle",
    });
  });

  // ---- stage 1: the cell
  s.addShape(pres.ShapeType.ellipse, {
    x: colX(0) + 0.52, y: 3.25, w: 1.75, h: 1.75,
    fill: { color: "EAF1F8" }, line: { color: "2C6E9B", width: 1.5 },
  });
  s.addShape(pres.ShapeType.ellipse, {
    x: colX(0) + 1.12, y: 3.85, w: 0.6, h: 0.55,
    fill: { color: "CFE0EC" }, line: { color: "2C6E9B", width: 1 },
  });
  s.addText("nucleus", {
    x: colX(0) + 1.12, y: 3.85, w: 0.6, h: 0.55, fontFace: BODY, fontSize: 7.5,
    color: "2C6E9B", align: "center", valign: "middle", margin: 0,
  });
  [[0.28, 0.34], [1.18, 0.2], [0.34, 1.2], [1.3, 1.28], [0.72, 0.14], [1.42, 0.75]].forEach(([dx, dy]) => {
    s.addShape(pres.ShapeType.ellipse, {
      x: colX(0) + 0.52 + dx, y: 3.25 + dy, w: 0.09, h: 0.09,
      fill: { color: "2C6E9B" }, line: { color: "2C6E9B", width: 0 },
    });
  });

  // ---- stages 2 and 3: the ranked gene list
  const genes = [
    ["CD3D", 0.78, false],
    ["IL7R", 0.60, false],
    ["TIGIT", 0.70, true],
    ["GZMB", 0.30, false],
    ["TOX", 0.20, false],
  ];
  const rowH = 0.38, rowY0 = 3.12;

  s.addText("most expressed", {
    x: colX(1), y: 2.78, w: cw, h: 0.26, fontFace: BODY, fontSize: 9.5, color: MUTED, margin: 0,
  });
  genes.forEach(([g, frac, hl], i) => {
    const y = rowY0 + i * rowH;
    s.addShape(pres.ShapeType.roundRect, {
      x: colX(1), y, w: cw, h: 0.32, rectRadius: 0.04,
      fill: { color: hl ? EOSIN_WASH : WHITE }, line: { color: hl ? EOSIN : RULE, width: hl ? 1.25 : 0.75 },
    });
    s.addText(g, {
      x: colX(1) + 0.12, y, w: 1.0, h: 0.32, fontFace: BODY, fontSize: 11.5,
      bold: hl, color: hl ? EOSIN : INK, valign: "middle", margin: 0,
    });
    const barMax = cw - 1.3;
    s.addShape(pres.ShapeType.roundRect, {
      x: colX(1) + 1.18, y: y + 0.09, w: barMax * frac, h: 0.14, rectRadius: 0.05,
      fill: { color: hl ? EOSIN : "2C6E9B" }, line: { color: hl ? EOSIN : "2C6E9B", width: 0 },
    });
  });
  s.addText("least expressed", {
    x: colX(1), y: 5.08, w: cw, h: 0.26, fontFace: BODY, fontSize: 9.5, color: MUTED, align: "right", margin: 0,
  });

  s.addText("DELETE — drop it out", {
    x: colX(2), y: 2.78, w: cw, h: 0.26, fontFace: BODY, fontSize: 10.5, bold: true, color: EOSIN, margin: 0,
  });
  genes.forEach(([g, , hl], i) => {
    const y = rowY0 + i * rowH;
    s.addShape(pres.ShapeType.roundRect, {
      x: colX(2), y, w: cw, h: 0.32, rectRadius: 0.04,
      fill: { color: hl ? "F7F7F8" : WHITE }, line: { color: RULE, width: 0.75, dashType: hl ? "dash" : "solid" },
    });
    s.addText(g, {
      x: colX(2) + 0.12, y, w: 1.2, h: 0.32, fontFace: BODY, fontSize: 11.5,
      color: hl ? MUTED : INK, strike: hl ? true : false, valign: "middle", margin: 0,
    });
    if (hl) {
      s.addText("✕", {
        x: colX(2) + cw - 0.42, y, w: 0.3, h: 0.32, fontFace: BODY, fontSize: 12,
        color: EOSIN, align: "center", valign: "middle", margin: 0,
      });
    } else if (i > 2) {
      s.addText("▲", {
        x: colX(2) + cw - 0.5, y, w: 0.35, h: 0.32, fontFace: BODY, fontSize: 9,
        color: "1A6B3C", align: "center", valign: "middle", margin: 0,
      });
    }
  });
  s.addText("the rest move up a rank", {
    x: colX(2), y: 5.08, w: cw, h: 0.26, fontFace: BODY, fontSize: 9.5, color: MUTED, margin: 0,
  });
  s.addText(
    [
      { text: "OVEREXPRESS", options: { bold: true, color: TEAL } },
      { text: " — move it to the top", options: { color: INK2 } },
    ],
    { x: colX(2), y: 5.32, w: cw, h: 0.26, fontFace: BODY, fontSize: 10.5, margin: 0 }
  );

  // ---- stage 4: the state map
  card(s, colX(3), 2.9, cw, 2.48, WHITE);
  s.addText("model's map of cell states", {
    x: colX(3) + 0.12, y: 2.98, w: cw - 0.24, h: 0.24, fontFace: BODY, fontSize: 9, color: MUTED, align: "center", margin: 0,
  });
  s.addShape(pres.ShapeType.ellipse, {
    x: colX(3) + 0.18, y: 4.12, w: 1.35, h: 0.88,
    fill: { color: HEMA_WASH }, line: { color: HEMA, width: 1.25, dashType: "dash" },
  });
  s.addText("SCLC", {
    x: colX(3) + 0.18, y: 4.12, w: 1.35, h: 0.88, fontFace: BODY, fontSize: 12, bold: true,
    color: HEMA, align: "center", valign: "middle", margin: 0,
  });
  s.addShape(pres.ShapeType.ellipse, {
    x: colX(3) + cw - 1.55, y: 3.20, w: 1.35, h: 0.86,
    fill: { color: "E8F2EC" }, line: { color: "1A6B3C", width: 1.25, dashType: "dash" },
  });
  s.addText("Normal", {
    x: colX(3) + cw - 1.55, y: 3.20, w: 1.35, h: 0.86, fontFace: BODY, fontSize: 12, bold: true,
    color: "1A6B3C", align: "center", valign: "middle", margin: 0,
  });
  s.addShape(pres.ShapeType.line, {
    x: colX(3) + 1.02, y: 3.80, w: cw - 2.12, h: 0.60,
    line: { color: EOSIN, width: 2, endArrowType: "triangle" }, flipV: true,
  });
  s.addShape(pres.ShapeType.ellipse, {
    x: colX(3) + 0.96, y: 4.38, w: 0.14, h: 0.14, fill: { color: HEMA }, line: { color: HEMA, width: 0 },
  });
  s.addText("shift toward Normal", {
    x: colX(3) + 0.12, y: 5.05, w: cw - 0.24, h: 0.26, fontFace: BODY, fontSize: 10, italic: true,
    color: EOSIN, align: "right", margin: 0,
  });

  // ---- caveat
  card(s, ML, 5.72, CW, 1.02, WASH);
  s.addText(
    [
      { text: "Nothing is edited in a laboratory here. ", options: { bold: true, color: INK } },
      { text: "The edit is made to the model's input and the result is the model's ", options: { color: INK2 } },
      { text: "predicted", options: { italic: true, color: INK2 } },
      { text: " change in cell state, averaged over thousands of real cells — a hypothesis generator, not a knockout experiment.", options: { color: INK2 } },
    ],
    { x: ML + 0.3, y: 5.72, w: CW - 0.6, h: 1.02, fontFace: BODY, fontSize: 12.5, lineSpacing: 19, valign: "middle", margin: 0 }
  );

  s.addNotes(
    "Walk the four stages in plain language, pausing on stage 2 — the key mental model is that Geneformer reads a cell as an ordered list of its genes, most-expressed first. It is a rank order, not expression values.\n\n" +
    "Deleting a gene drops it out of the list and everything below moves up a rank; overexpressing moves it to the top. Re-reading the modified cell gives a predicted change in cell state, averaged over thousands of real cells.\n\n" +
    "Stage 4 is the measurement: the model places cells on a map of disease states, and the perturbation shifts a cell along that map — here, toward the normal T-cell state.\n\n" +
    "Say the caveat out loud: nothing is edited in a laboratory here. This is the conceptual keystone for a non-modelling audience, so do not rush it."
  );
}

/* ======== SLIDE 3 — the hit rule and what a shift means =============== */
{
  const s = pres.addSlide();
  kicker(s, "How a hit is called");
  title(s, "A gene counts only when deleting and overexpressing disagree");
  footer(s, 3, "Bidirectional concordance");

  s.addText(
    "Every number in the next four slides is a perturbation shift. Two rules decide which ones are allowed to count.",
    { x: ML, y: 1.45, w: CW, h: 0.35, fontFace: BODY, fontSize: 14, color: INK2, margin: 0, valign: "top" }
  );

  const lw = 7.15;

  // concordant case
  card(s, ML, 1.95, lw, 1.75, "EAF4F3");
  s.addText("CONCORDANT — counts as a hit", {
    x: ML + 0.28, y: 2.08, w: lw - 0.56, h: 0.3, fontFace: BODY, fontSize: 11.5, bold: true,
    charSpacing: 1.2, color: TEAL, margin: 0,
  });
  [["Delete the gene", "toward Normal", "1A6B3C", 2.48], ["Overexpress it", "away from Normal", EOSIN, 3.05]].forEach(([a, b, col, y]) => {
    s.addText(a, {
      x: ML + 0.28, y, w: 2.0, h: 0.4, fontFace: BODY, fontSize: 12.5, bold: true, color: INK, valign: "middle", margin: 0,
    });
    s.addShape(pres.ShapeType.line, {
      x: ML + 2.4, y: y + 0.2, w: 1.5, h: 0,
      line: { color: col, width: 2, endArrowType: "triangle" },
    });
    s.addText(b, {
      x: ML + 4.05, y, w: 3.0, h: 0.4, fontFace: BODY, fontSize: 12.5, color: col, valign: "middle", margin: 0,
    });
  });

  // discordant case
  card(s, ML, 3.85, lw, 1.6, WASH);
  s.addText("NOT CONCORDANT — discarded", {
    x: ML + 0.28, y: 3.98, w: lw - 0.56, h: 0.3, fontFace: BODY, fontSize: 11.5, bold: true,
    charSpacing: 1.2, color: MUTED, margin: 0,
  });
  s.addText("Both edits move the cell the same way — the model is not responding to the direction of the change, so the gene is not treated as a candidate.", {
    x: ML + 0.28, y: 4.35, w: lw - 0.56, h: 0.9, fontFace: BODY, fontSize: 12.5, color: INK2, margin: 0, valign: "top",
  });

  card(s, ML, 5.6, lw, 1.15, HEMA_WASH);
  s.addText(
    [
      { text: "And it must replicate. ", options: { bold: true, color: HEMA } },
      { text: "The four headline edits hold the same sign in all three SCLC donors, at FDR < 0.05 in both arms, in genes detected in at least 100 source cells.", options: { color: INK } },
    ],
    { x: ML + 0.3, y: 5.6, w: lw - 0.6, h: 1.15, fontFace: BODY, fontSize: 12.5, lineSpacing: 19, valign: "middle", margin: 0 }
  );

  // right: what the numbers are
  const rx = ML + lw + 0.45;
  const rw = CW - lw - 0.45;
  s.addText("What the numbers are", {
    x: rx, y: 1.95, w: rw, h: 0.3, fontFace: BODY, fontSize: 14, bold: true, color: NAVY, margin: 0,
  });

  const facts = [
    ["A change in similarity", "How much closer the perturbed cell sits to a target disease state, on the model's own scale."],
    ["Small by construction", "Shifts run ~0.001–0.01. Size is not significance — slide 7 shows why."],
    ["Averaged over cells", "Each value pools thousands of real T cells, not one simulated cell."],
  ];
  let fy = 2.4;
  facts.forEach(([h, d]) => {
    card(s, rx, fy, rw, 1.32, WHITE);
    s.addText(h, {
      x: rx + 0.24, y: fy + 0.14, w: rw - 0.48, h: 0.3, fontFace: BODY, fontSize: 12.5, bold: true, color: TEAL, margin: 0,
    });
    s.addText(d, {
      x: rx + 0.24, y: fy + 0.46, w: rw - 0.48, h: 0.75, fontFace: BODY, fontSize: 11, color: INK2, margin: 0, valign: "top",
    });
    fy += 1.45;
  });

  s.addNotes(
    "This slide is the bridge between the method and every number that follows. Two rules decide what counts.\n\n" +
    "Rule one is bidirectional concordance: a gene only counts when deleting it and overexpressing it move cells in OPPOSITE directions. If both edits push the cell the same way, the model is not responding to the direction of the change, and the gene is discarded. This is the internal consistency check, and it is stricter than a significance threshold.\n\n" +
    "Rule two is replication: the four headline edits hold the same sign in all three SCLC donors, at FDR < 0.05 in both arms, with at least 100 detected source cells.\n\n" +
    "On the right, be explicit about what a 'shift' is: a change in similarity on the model's own scale, running about 0.001 to 0.01, pooled over thousands of real cells. Small numbers are expected — and slide 7 explains why the big ones are the suspicious ones."
  );
}

/* ================= SLIDE 4 — cohort & donor guard ===================== */
{
  const s = pres.addSlide();
  kicker(s, "Data provenance");
  title(s, "42 individuals, zero donor leakage — so accuracy is measured on people the model never saw");
  footer(s, 4, "Cohort architecture");

  const lw = 6.9;
  s.addTable(
    [
      [
        { text: "Disease", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
        { text: "Train", options: { bold: true, color: WHITE, fill: { color: NAVY }, align: "right" } },
        { text: "Eval", options: { bold: true, color: WHITE, fill: { color: NAVY }, align: "right" } },
        { text: "Test", options: { bold: true, color: WHITE, fill: { color: NAVY }, align: "right" } },
        { text: "Donors", options: { bold: true, color: WHITE, fill: { color: NAVY }, align: "right" } },
      ],
      [
        { text: "SCLC", options: { bold: true } },
        { text: "7,037", options: { align: "right" } },
        { text: "2,330", options: { align: "right" } },
        { text: "2,424", options: { align: "right" } },
        { text: "19", options: { align: "right", bold: true } },
      ],
      [
        { text: "LUAD" },
        { text: "17,831", options: { align: "right" } },
        { text: "5,611", options: { align: "right" } },
        { text: "6,387", options: { align: "right" } },
        { text: "22", options: { align: "right", bold: true } },
      ],
      [
        { text: "Normal" },
        { text: "2,334", options: { align: "right" } },
        { text: "1,620", options: { align: "right" } },
        { text: "566", options: { align: "right" } },
        { text: "4", options: { align: "right", bold: true } },
      ],
    ],
    {
      x: ML, y: 1.75, w: lw, colW: [1.9, 1.3, 1.2, 1.2, 1.3],
      fontFace: BODY, fontSize: 13, color: INK, border: { type: "solid", color: RULE, pt: 0.75 },
      rowH: 0.42, valign: "middle",
    }
  );

  s.addText(
    [
      { text: "45 donor×disease groups = 42 people + 3 paired tumour–normal donors", options: { bold: true, color: INK, breakLine: true } },
      { text: "A donor contributing both tumour and normal tissue appears twice. The three extra rows are exactly RU675, RU682 and RU684.", options: { color: INK2 } },
    ],
    { x: ML, y: 3.6, w: lw, h: 1.0, fontFace: BODY, fontSize: 12.5, margin: 0, valign: "top", lineSpacing: 18 }
  );

  s.addText("46,140 cells  ·  24,540 features  ·  median 390–1,172 cells per donor", {
    x: ML, y: 4.7, w: lw, h: 0.3, fontFace: BODY, fontSize: 11.5, color: MUTED, margin: 0,
  });

  // right: donor split schematic
  const rx = ML + lw + 0.45;
  const rw = CW - lw - 0.45;
  card(s, rx, 1.75, rw, 3.85, WASH);
  s.addText("Donor-disjoint splits", {
    x: rx + 0.3, y: 1.95, w: rw - 0.6, h: 0.32, fontFace: BODY, fontSize: 14, bold: true, color: NAVY, margin: 0,
  });

  const blocks = [["Train", 27, TEAL], ["Eval", 7, HEMA_LT], ["Test", 8, EOSIN]];
  let by = 2.45;
  blocks.forEach(([lab, n, col]) => {
    s.addShape(pres.ShapeType.ellipse, { x: rx + 0.32, y: by + 0.05, w: 0.26, h: 0.26, fill: { color: col }, line: { color: col, width: 0 } });
    s.addText(`${lab} — ${n} donors`, {
      x: rx + 0.7, y: by, w: rw - 1.0, h: 0.36, fontFace: BODY, fontSize: 13, bold: true, color: INK, margin: 0, valign: "middle",
    });
    by += 0.46;
  });

  s.addText("No individual appears in two splits. Paired tumour–normal donors are confined to one split by a cross-disease guard.", {
    x: rx + 0.3, y: 4.0, w: rw - 0.6, h: 1.0, fontFace: BODY, fontSize: 12, color: INK2, margin: 0, valign: "top",
  });

  card(s, rx, 5.75, rw, 0.62, "E8F2EC");
  s.addText("Leakage audit: PASS  ·  0 / 42 donors", {
    x: rx, y: 5.75, w: rw, h: 0.62, fontFace: BODY, fontSize: 13, bold: true,
    color: "1A6B3C", align: "center", valign: "middle", margin: 0,
  });

  s.addText("Held-out test set:  9,377 cells  from  8 donors", {
    x: ML, y: 5.75, w: lw, h: 0.62, fontFace: HEAD, fontSize: 16, bold: true, color: NAVY, margin: 0, valign: "middle",
  });

  s.addNotes(
    "Provenance first: HTAN MSK collection via the CELLxGENE 'T cells' dataset — 46,140 CD4+CD8 T cells x 24,540 features from 42 individuals.\n\n" +
    "Explain the 45-vs-42 arithmetic before anyone asks: the split table has one row per donor x disease, so a donor contributing both tumour and normal tissue appears twice. The three-row difference is exactly the three paired tumour-normal donors, RU675, RU682, RU684.\n\n" +
    "The credibility point: splits are donor-disjoint and audited — zero leakage across all 42 individuals, with a cross-disease guard for the paired donors. Every performance number downstream is measured on people the model never saw. This is checked by assertion at build time, not by trusting a stored PASS string."
  );
}

/* ==================== SLIDE 5 — model performance ===================== */
{
  const s = pres.addSlide();
  kicker(s, "Model");
  title(s, "91.9% accuracy on held-out donors — with an SCLC→LUAD confusion worth noticing");
  footer(s, 5, "Fine-tuning & performance");

  const sw = (CW - 3 * 0.3) / 4;
  stat(s, ML + 0 * (sw + 0.3), 1.62, sw, "91.9%", "Test accuracy", TEAL);
  stat(s, ML + 1 * (sw + 0.3), 1.62, sw, "0.903", "Macro F1", TEAL);
  stat(s, ML + 2 * (sw + 0.3), 1.62, sw, "9,377", "Held-out cells", NAVY);
  stat(s, ML + 3 * (sw + 0.3), 1.62, sw, "8", "Test donors", NAVY);

  s.addImage({ path: A("confusion.png"), x: ML, y: 2.85, w: 5.9, h: 3.74 });
  cap(s, "Held-out test confusion matrix, row-normalised.", ML, 6.66, 5.9);

  const rx = ML + 6.55;
  const rw = CW - 6.55;

  card(s, rx, 2.85, rw, 1.75, EOSIN_WASH);
  s.addText("Look at the diagonal", {
    x: rx + 0.28, y: 3.0, w: rw - 0.56, h: 0.3, fontFace: BODY, fontSize: 13.5, bold: true, color: EOSIN, margin: 0,
  });
  bullets(s, [
    "Normal 98.6%  ·  LUAD 95.9%",
    "SCLC 80.0% — 484 missed cells go to LUAD",
    "Same axis as the open question on slide 7",
  ], rx + 0.28, 3.38, rw - 0.56, 1.1, { size: 12.5, color: INK });

  card(s, rx, 4.78, rw, 2.05, WHITE);
  s.addText("Fine-tuning setup", {
    x: rx + 0.28, y: 4.93, w: rw - 0.56, h: 0.3, fontFace: BODY, fontSize: 13.5, bold: true, color: NAVY, margin: 0,
  });
  bullets(s, [
    "Geneformer-V2-104M",
    "1 epoch · LR 5e-5 · 6 frozen layers",
    "NVIDIA GB10, 119 GiB unified memory",
    "~48 h GPU · ~5.88M cell-gene perturbations",
  ], rx + 0.28, 5.3, rw - 0.56, 1.45, { size: 12, color: INK2, after: 5 });

  s.addNotes(
    "Headline the metrics: 91.9% accuracy, macro F1 0.903 on 9,377 held-out cells from 8 donors. Donor-held-out performance separates SCLC, Normal and LUAD well enough to support the downstream perturbation screen.\n\n" +
    "Do not skip the diagonal. Normal recall is 98.6% and LUAD 95.9%, but SCLC is 80.0% — and the 20% that is missed goes almost entirely to LUAD, 484 cells. Flag it yourself before a reviewer does.\n\n" +
    "Preview the tension: that SCLC-to-LUAD confusion reappears on slide 7 as a directional result the model produces, and it is the one thing you have not reconciled. Planting it here makes slide 7 land as rigour rather than as a hole."
  );
}

/* ================= SLIDE 6 — applied checkpoint screen ================ */
{
  const s = pres.addSlide();
  kicker(s, "Screen");
  title(s, "Four edits replicate across every donor: TIM-3, TIGIT, CTLA-4, IL7R");
  footer(s, 6, "Applied checkpoint screen");

  const sw = (CW - 3 * 0.3) / 4;
  stat(s, ML + 0 * (sw + 0.3), 1.6, sw, "50", "Genes screened", NAVY);
  stat(s, ML + 1 * (sw + 0.3), 1.6, sw, "4", "Replicated edits", TEAL);
  stat(s, ML + 2 * (sw + 0.3), 1.6, sw, "3/3", "Donors agree", TEAL);
  stat(s, ML + 3 * (sw + 0.3), 1.6, sw, "8", "No deletion result", MUTED);

  const tw = 7.15;
  s.addTable(
    [
      [
        { text: "Gene", options: { bold: true, color: WHITE, fill: { color: NAVY } } },
        { text: "Delete → Normal", options: { bold: true, color: WHITE, fill: { color: NAVY }, align: "right" } },
        { text: "Overexpress", options: { bold: true, color: WHITE, fill: { color: NAVY }, align: "right" } },
        { text: "Detected", options: { bold: true, color: WHITE, fill: { color: NAVY }, align: "right" } },
      ],
      [
        { text: "HAVCR2 (TIM-3)", options: { bold: true } },
        { text: "+0.0019", options: { align: "right", color: TEAL, bold: true } },
        { text: "−0.0060", options: { align: "right", color: EOSIN } },
        { text: "202", options: { align: "right" } },
      ],
      [
        { text: "TIGIT", options: { bold: true } },
        { text: "+0.0018", options: { align: "right", color: TEAL, bold: true } },
        { text: "−0.0082", options: { align: "right", color: EOSIN } },
        { text: "438", options: { align: "right" } },
      ],
      [
        { text: "CTLA-4", options: { bold: true } },
        { text: "+0.0014", options: { align: "right", color: TEAL, bold: true } },
        { text: "−0.0014", options: { align: "right", color: EOSIN } },
        { text: "282", options: { align: "right" } },
      ],
      [
        { text: "IL7R (persistence)", options: { bold: true } },
        { text: "+0.0013", options: { align: "right", color: TEAL, bold: true } },
        { text: "−0.0054", options: { align: "right", color: EOSIN } },
        { text: "1,131", options: { align: "right" } },
      ],
    ],
    {
      x: ML, y: 2.85, w: tw, colW: [2.35, 1.75, 1.55, 1.5],
      fontFace: BODY, fontSize: 12.5, color: INK,
      border: { type: "solid", color: RULE, pt: 0.75 }, rowH: 0.4, valign: "middle",
    }
  );

  s.addText("All four: same sign in all 3 SCLC donors, FDR < 0.05 in both arms, ≥100 source cells detected.", {
    x: ML, y: 5.1, w: tw, h: 0.5, fontFace: BODY, fontSize: 12, color: INK2, margin: 0, valign: "top",
  });

  // filter funnel
  card(s, ML, 5.7, tw, 1.1, WASH);
  const funnel = ["50 genes", "concordant\nboth arms", "FDR < 0.05\nboth arms", "all 3\ndonors", "≥100\ndetected", "4"];
  const fw = (tw - 0.4) / funnel.length;
  funnel.forEach((t, i) => {
    const last = i === funnel.length - 1;
    s.addText(t, {
      x: ML + 0.2 + i * fw, y: 5.78, w: fw - 0.26, h: 0.9,
      fontFace: BODY, fontSize: last ? 17 : 9.5, bold: last,
      color: last ? TEAL : INK2, align: "center", valign: "middle", margin: 0,
    });
    if (!last) {
      s.addText("›", {
        x: ML + 0.2 + (i + 1) * fw - 0.24, y: 5.78, w: 0.22, h: 0.9,
        fontFace: BODY, fontSize: 13, color: "AEBCCB", align: "center", valign: "middle", margin: 0,
      });
    }
  });

  const rx = ML + tw + 0.4;
  const rw = CW - tw - 0.4;
  const bw = 3.72, bh = bw * 0.8834;
  s.addImage({ path: A("screen_b.png"), x: rx + (rw - bw) / 2, y: 2.85, w: bw, h: bh });
  cap(s, "SCLC → normal: four knockouts replicate across all 3 donors.", rx, 6.20, rw);

  s.addNotes(
    "Scope the screen: 50 checkpoint and T-cell engineering genes, tested on the SCLC to Normal transition. Of 123 concordant hits genome-wide, 43 are fully donor-consistent; within this targeted panel, four edits survive every filter.\n\n" +
    "State the conjunction, because it is stricter than 'biggest effect': concordant in both perturbation arms, FDR < 0.05 in both, same sign in all three SCLC donors, and detected in at least 100 source cells. These rows are selected in code by that conjunction, so the table cannot drift from the criteria printed beside it.\n\n" +
    "Read the biology plainly: deleting TIM-3, TIGIT, CTLA-4 or IL7R moves SCLC T cells toward a normal T-cell state; overexpressing moves them away. Then set up slide 7 — before you believe the ranking, one control changes how you read these numbers."
  );
}

/* ========= SLIDE 7 — detection confound + open question =============== */
{
  const s = pres.addSlide();
  kicker(s, "Quality control", EOSIN);
  title(s, "Sparse genes fake big effects (ρ = −0.60) — and one model result still contradicts the clinic");
  footer(s, 7, "Detection confound & the open axis");

  const lw = 7.0;
  s.addImage({ path: A("detection.png"), x: ML, y: 1.75, w: lw, h: lw * 0.5036 });
  cap(s, "Detection vs absolute deletion effect across the full 50-gene panel.", ML, 5.34, lw);

  bullets(s, [
    "Detection vs effect size: Spearman ρ = −0.60 (p = 2.7×10⁻⁵, n = 42)",
    "8 of 50 genes have no deletion result — undetected, so undefined, not zero",
    "18 of 50 fall below 100 detected SCLC source cells",
    "Effect size alone is not evidence — which is why slide 6 ranked on replication",
  ], ML, 5.62, lw, 1.32, { size: 12, color: INK, after: 3 });

  // open question
  const rx = ML + lw + 0.4;
  const rw = CW - lw - 0.4;
  card(s, rx, 1.75, rw, 5.05, EOSIN_WASH);
  s.addText("OPEN QUESTION", {
    x: rx + 0.28, y: 1.95, w: rw - 0.56, h: 0.32, fontFace: BODY, fontSize: 12, bold: true,
    charSpacing: 1.6, color: EOSIN, margin: 0,
  });
  s.addText("The model orders the exhaustion axis", {
    x: rx + 0.28, y: 2.32, w: rw - 0.56, h: 0.32, fontFace: BODY, fontSize: 13, color: INK, margin: 0,
  });

  // axis graphic
  const axY = 3.1;
  s.addShape(pres.ShapeType.line, {
    x: rx + 0.45, y: axY, w: rw - 0.9, h: 0, line: { color: NAVY, width: 1.75 },
  });
  ["Normal", "SCLC", "LUAD"].forEach((lab, i) => {
    const cx = rx + 0.45 + i * ((rw - 0.9) / 2);
    s.addShape(pres.ShapeType.ellipse, {
      x: cx - 0.09, y: axY - 0.09, w: 0.18, h: 0.18,
      fill: { color: i === 1 ? EOSIN : NAVY }, line: { color: WHITE, width: 1 },
    });
    s.addText(lab, {
      x: cx - 0.6, y: axY + 0.16, w: 1.2, h: 0.28, fontFace: BODY, fontSize: 11.5,
      bold: i === 1, color: i === 1 ? EOSIN : INK2, align: "center", margin: 0,
    });
  });

  s.addText("TIGIT overexpression moves SCLC cells away from Normal (−0.0082) but toward LUAD (+0.0256).", {
    x: rx + 0.28, y: 3.75, w: rw - 0.56, h: 0.95, fontFace: BODY, fontSize: 12.5, color: INK, margin: 0, valign: "top",
  });

  card(s, rx + 0.28, 4.75, rw - 0.56, 0.72, WHITE);
  s.addText("Model axis: T-cell intrinsic", {
    x: rx + 0.42, y: 4.75, w: rw - 0.84, h: 0.72, fontFace: BODY, fontSize: 12, bold: true, color: TEAL, valign: "middle", margin: 0,
  });
  card(s, rx + 0.28, 5.55, rw - 0.56, 0.72, WHITE);
  s.addText("Clinic: tumour level — SCLC cold, ICI-resistant", {
    x: rx + 0.42, y: 5.55, w: rw - 0.84, h: 0.72, fontFace: BODY, fontSize: 12, bold: true, color: EOSIN, valign: "middle", margin: 0,
  });
  s.addText("Not the same measurement. Presented as an open question, not a reconciled result.", {
    x: rx + 0.28, y: 6.28, w: rw - 0.56, h: 0.48, fontFace: BODY, fontSize: 11, italic: true, color: INK2, margin: 0, valign: "top",
  });

  s.addNotes(
    "The control that reframes the screen: across the 50-gene panel, detection count and nominal effect size are strongly anti-correlated, Spearman rho = -0.60 (p = 2.7e-5, n = 42). The fewer cells a gene is detected in, the larger its apparent shift. Effect magnitude alone is therefore not a ranking criterion — which is exactly why slide 6 ranked on replication and detection instead.\n\n" +
    "Reinforce with the gaps: 8 of 50 genes have no deletion result at all — undetected in SCLC source cells, so deletion is undefined, not zero — and 18 fall below 100 detected cells.\n\n" +
    "Then the honest part, and give it real time: on the exhaustion axis the model orders states Normal < SCLC < LUAD. TIGIT overexpression moves SCLC cells away from Normal but toward LUAD. That runs opposite to the clinical picture of SCLC as a cold, ICI-resistant tumour. Present it as an open question — a tumour-level clinical phenotype and a T-cell-intrinsic model axis are not the same measurement."
  );
}

/* ================== SLIDE 8 — CAR-T relevance ========================= */
{
  const s = pres.addSlide();
  kicker(s, "Translational scope");
  title(s, "Two of the four hits are live CAR-T engineering targets — the tumour antigens are out of reach here");
  footer(s, 8, "CAR-T relevance and design scope");

  const lw = 6.3;
  s.addText("How interpretable is each hit?", {
    x: ML, y: 1.65, w: lw, h: 0.32, fontFace: BODY, fontSize: 14, bold: true, color: NAVY, margin: 0,
  });

  const hits = [
    ["TIGIT", "Strongest external CAR-T support", TEAL, "Exhaustion checkpoint"],
    ["TIM-3 (HAVCR2)", "Direct KO precedent in solid tumours", TEAL, "Exhaustion checkpoint"],
    ["CTLA-4", "Real signal, less CAR-T-specific evidence", HEMA_LT, "Exhaustion checkpoint"],
    ["IL7R", "Persistence / fitness — not a checkpoint", EOSIN, "Persistence target"],
  ];
  let hy = 2.08;
  hits.forEach(([g, d, col, kind]) => {
    card(s, ML, hy, lw, 1.02, WHITE);
    s.addShape(pres.ShapeType.ellipse, {
      x: ML + 0.22, y: hy + 0.29, w: 0.44, h: 0.44, fill: { color: col }, line: { color: col, width: 0 },
    });
    s.addText(g, {
      x: ML + 0.82, y: hy + 0.12, w: lw - 3.15, h: 0.35, fontFace: HEAD, fontSize: 15, bold: true, color: NAVY, margin: 0, valign: "middle",
    });
    s.addText(d, {
      x: ML + 0.82, y: hy + 0.47, w: lw - 1.1, h: 0.42, fontFace: BODY, fontSize: 11.5, color: INK2, margin: 0, valign: "top",
    });
    s.addText(kind, {
      x: ML + lw - 2.05, y: hy + 0.12, w: 1.85, h: 0.32, fontFace: BODY, fontSize: 9.5,
      color: col, align: "right", margin: 0, valign: "middle",
    });
    hy += 1.12;
  });

  // can / cannot
  const rx = ML + lw + 0.45;
  const rw = CW - lw - 0.45;

  card(s, rx, 1.98, rw, 2.2, "EAF4F3");
  s.addText("This design CAN test", {
    x: rx + 0.28, y: 2.14, w: rw - 0.56, h: 0.3, fontFace: BODY, fontSize: 13.5, bold: true, color: TEAL, margin: 0,
  });
  bullets(s, [
    "T-cell-intrinsic checkpoint edits",
    "Persistence receptors (IL7R)",
    "Exhaustion-programme shifts in SCLC T cells",
  ], rx + 0.28, 2.52, rw - 0.56, 1.5, { size: 12.5, color: INK });

  card(s, rx, 4.35, rw, 2.45, WASH);
  s.addText("This design CANNOT test", {
    x: rx + 0.28, y: 4.51, w: rw - 0.56, h: 0.3, fontFace: BODY, fontSize: 13.5, bold: true, color: MUTED, margin: 0,
  });
  bullets(s, [
    "CAR-T tumour antigens: DLL3, SEZ6, NCAM1, CD276, CEACAM5",
    "PD-L1 (CD274) — absent from the T-cell atlas",
  ], rx + 0.28, 4.89, rw - 0.56, 1.2, { size: 12.5, color: INK2 });
  s.addText("Not expressed by T cells — this is scope, not failure. No PD-L1 claim is made.", {
    x: rx + 0.28, y: 6.05, w: rw - 0.56, h: 0.6, fontFace: BODY, fontSize: 11, italic: true, color: MUTED, margin: 0, valign: "top",
  });

  s.addNotes(
    "Separate the four hits by how interpretable they are, because they are not equivalent: TIGIT has the strongest external CAR-T support; TIM-3 has direct knockout precedent in solid-tumour models; CTLA-4 is a real signal with less CAR-T-specific evidence; IL7R is a persistence/fitness target, not an exhaustion checkpoint at all. Grouping them as 'four checkpoints' would be wrong.\n\n" +
    "The structural limit, stated as design not apology: this is a T-cell atlas, so it can only test T-cell-intrinsic edits. The CAR-T tumour antigens — DLL3, SEZ6, NCAM1, CD276, CEACAM5 — are not perturbable here because they are not expressed by T cells. Likewise PD-L1 is absent from the atlas, so this screen cannot speak to PD-L1 biology.\n\n" +
    "Land the translational read: the actionable axis is checkpoint knockout in the T-cell product (TIGIT, HAVCR2) plus persistence engineering (IL7R) — which is precisely the CRISPR experiment on slide 12."
  );
}

/* ================= SLIDE 9 — spatial validation ======================= */
{
  const s = pres.addSlide();
  kicker(s, "Tissue validation", EOSIN);
  title(s, "In intact tumour tissue, dysfunction tracks T-cell abundance — and antigen presentation tracks it hardest");
  footer(s, 9, "GSE263196 Visium validation");

  // The figure column is sized to the composite's own aspect ratio (2000x1352)
  // so the image fills its box exactly. Sizing it by height alone letterboxed
  // 2.26" of dead width and left the slide-mount strip far wider than the
  // picture under it. The separate forest plot is gone: panel C of this
  // composite is the same per-specimen forest with the same pooled estimate.
  const iw = 4.3;
  const ih = 2.78;               // 4.3 * (1287 / 1989), matches the trimmed composite

  s.addShape(pres.ShapeType.roundRect, {
    x: ML, y: 1.68, w: iw, h: 0.32, rectRadius: 0.03,
    fill: { color: NAVY }, line: { color: NAVY, width: 0 },
  });
  s.addText("GSE263196  ·  10x Visium  ·  5 specimens  ·  15,632 spots", {
    x: ML + 0.12, y: 1.68, w: iw - 0.24, h: 0.32, fontFace: BODY, fontSize: 8.5,
    color: "CFE0EC", valign: "middle", margin: 0,
  });
  s.addImage({ path: A("spatial.png"), x: ML, y: 2.0, w: iw, h: ih });
  cap(s, "A: T-cell identity.  B: dysfunction.  C: per-specimen and pooled ρ.", ML, 4.97, iw);

  card(s, ML, 5.45, iw, 1.32, WASH);
  s.addText("Why this is independent", {
    x: ML + 0.24, y: 5.58, w: iw - 0.48, h: 0.28, fontFace: BODY, fontSize: 12, bold: true, color: NAVY, margin: 0,
  });
  bullets(s, [
    "Different cohort, assay and lab",
    "15,632 / 15,774 in-tissue spots kept (99.1%)",
    "Marker panels share no genes — not a self-correlation",
  ], ML + 0.24, 5.9, iw - 0.48, 0.85, { size: 10, color: INK2, after: 2, lineSpacing: 14 });

  const rx = ML + iw + 0.45;
  const rw = CW - iw - 0.45;

  card(s, rx, 1.68, rw, 1.12, WASH);
  s.addText(
    [
      { text: "ρ = 0.161", options: { fontSize: 22, bold: true, color: NAVY, fontFace: HEAD, breakLine: true } },
      { text: "7-gene dysfunction benchmark  ·  [0.146, 0.176]  ·  4 of 5 specimens p < 0.001", options: { fontSize: 11, color: INK2 } },
    ],
    { x: rx + 0.26, y: 1.68, w: rw - 0.52, h: 1.12, valign: "middle", margin: 0, fontFace: BODY }
  );

  card(s, rx, 2.95, rw, 1.12, HEMA_WASH);
  s.addText(
    [
      { text: "ρ = 0.361 → 0.384", options: { fontSize: 22, bold: true, color: HEMA, fontFace: HEAD, breakLine: true } },
      { text: "13-gene antigen-presentation programme  ·  depth-controlled  ·  7.95 σ above null", options: { fontSize: 11, color: INK2 } },
    ],
    { x: rx + 0.26, y: 2.95, w: rw - 0.52, h: 1.12, valign: "middle", margin: 0, fontFace: BODY }
  );

  s.addChart(
    pres.ChartType.bar,
    [{
      name: "Raw rho",
      labels: ["Cytotoxic effector", "Antigen presentation", "Interferon", "Checkpoint/exhaustion", "Dysfunction benchmark"],
      values: [0.413, 0.361, 0.280, 0.218, 0.161],
    }],
    {
      x: rx, y: 4.22, w: rw, h: 2.35,
      barDir: "bar", chartColors: [TEAL, HEMA, "8FA6B8", "AEBCCB", "C9D3DE"],
      varyColors: true, showLegend: false,
      showTitle: true, title: "Spatial ρ by programme", titleFontSize: 12.5, titleColor: NAVY, titleFontFace: BODY,
      showValue: true, dataLabelPosition: "outEnd", dataLabelFormatCode: "0.000",
      dataLabelFontSize: 10, dataLabelColor: INK2, dataLabelFontFace: BODY,
      catAxisLabelColor: INK2, catAxisLabelFontSize: 10.5, catAxisLabelFontFace: BODY,
      valAxisHidden: true, valGridLine: { style: "none" }, catGridLine: { style: "none" },
      barGapWidthPct: 45, valAxisMaxVal: 0.5,
    }
  );

  s.addText("Modest effects, high between-specimen heterogeneity (I² ≈ 99%) — real differences in microenvironment composition, not measurement noise.", {
    x: rx, y: 6.50, w: rw, h: 0.48, fontFace: BODY, fontSize: 10.5, italic: true, color: MUTED, margin: 0, valign: "top",
  });

  s.addNotes(
    "Set up the independence: GSE263196 is a different cohort, a different assay and a different lab — five fresh-frozen SCLC sections, 15,632 of 15,774 in-tissue spots retained (99.1%), all 21 marker genes present in every specimen. No donor, cell or read is shared with the discovery atlas.\n\n" +
    "Explain the measurement in one breath: at every spatially indexed spot you score two disjoint marker panels — T-cell abundance and T-cell dysfunction — then correlate them within each specimen. The panels share no genes, so the correlation cannot be a self-correlation.\n\n" +
    "Walk the figure top to bottom: row A is the T-cell identity score, row B the dysfunction score, row C the per-specimen effects with the pooled estimate.\n\n" +
    "Give both numbers, correctly attributed. The benchmark 7-gene dysfunction panel gives pooled rho = 0.161 [0.146, 0.176], significant in 4 of 5 specimens. The 13-gene antigen-presentation programme is stronger at rho = 0.361, rising to 0.384 controlling for sequencing depth, and sits 7.95 sigma above an expression-matched random null. Depth control RAISING the estimate is the opposite of what a technical artifact does.\n\n" +
    "Be straight about effect size: modest correlations with high between-specimen heterogeneity, I-squared about 99%, consistent with real differences in tumour microenvironment composition rather than measurement noise."
  );
}

/* ============ SLIDE 10 — mechanism & literature overlay =============== */
{
  const s = pres.addSlide();
  kicker(s, "Mechanism");
  title(s, "The hits sit in one connected interaction neighbourhood — from prior knowledge, not inferred here");
  footer(s, 10, "STRING literature overlay");

  const iw = 6.5;
  s.addImage({ path: A("network_c.png"), x: ML, y: 1.7, w: iw, h: iw * 0.6726 });
  cap(s, "All 11 ICI/CAR-T candidates on the STRING map; colour = SCLC→Normal deletion shift.", ML, 6.12, iw);

  const rx = ML + iw + 0.45;
  const rw = CW - iw - 0.45;

  s.addText("Evidence behind each hit", {
    x: rx, y: 1.7, w: rw, h: 0.32, fontFace: BODY, fontSize: 14, bold: true, color: NAVY, margin: 0,
  });

  const ev = [
    ["TIGIT", "CAR-T support", 0.95, TEAL],
    ["TIM-3", "Solid-tumour KO precedent", 0.9, TEAL],
    ["CTLA-4", "Real signal, less CAR-T-specific", 0.55, HEMA_LT],
    ["IL7R", "Persistence, not exhaustion", 0.45, EOSIN],
  ];
  let ey = 2.12;
  ev.forEach(([g, d, frac, col]) => {
    s.addText(g, {
      x: rx, y: ey, w: 1.5, h: 0.28, fontFace: BODY, fontSize: 13, bold: true, color: NAVY, margin: 0, valign: "middle",
    });
    s.addShape(pres.ShapeType.roundRect, {
      x: rx + 1.5, y: ey + 0.07, w: rw - 1.5, h: 0.15, rectRadius: 0.07,
      fill: { color: "E4EAF0" }, line: { color: "E4EAF0", width: 0 },
    });
    s.addShape(pres.ShapeType.roundRect, {
      x: rx + 1.5, y: ey + 0.07, w: (rw - 1.5) * frac, h: 0.15, rectRadius: 0.07,
      fill: { color: col }, line: { color: col, width: 0 },
    });
    s.addText(d, {
      x: rx, y: ey + 0.3, w: rw, h: 0.3, fontFace: BODY, fontSize: 10.5, color: INK2, margin: 0, valign: "top",
    });
    ey += 0.78;
  });

  card(s, rx, 5.4, rw, 1.35, WASH);
  bullets(s, [
    "54 STRING edges among 16 context genes",
    "TOX, LAYN: no edge above threshold",
    "PD-L1 absent from the T-cell atlas",
  ], rx + 0.24, 5.56, rw - 0.48, 1.05, { size: 11, color: INK2, after: 3 });

  s.addText("Prior interaction evidence — not a network inferred from these cells.", {
    x: ML, y: 6.55, w: iw, h: 0.32, fontFace: BODY, fontSize: 12, bold: true, color: EOSIN, margin: 0, valign: "middle",
  });

  s.addNotes(
    "Show that the four replicated edits are not scattered: all 11 ICI/CAR-T candidates map onto a connected STRING neighbourhood, 54 filtered edges among 16 context genes, with TIM-3, TIGIT, CTLA-4 and IL7R sitting close to the cytotoxic and interferon machinery.\n\n" +
    "State the epistemics without hedging: this map is a literature overlay, prior interaction evidence — it is NOT a gene-to-gene network inferred from these cells. It provides context for the hits; it is not independent evidence for them.\n\n" +
    "Point out the informative absences: TOX and LAYN carry no STRING edge above threshold, and PD-L1 is not in the T-cell atlas at all. TIM-3 is the strongest individual hit but has no supported non-text-mining edge to another checkpoint here — an explicit asymmetry, not a plotting omission."
  );
}

/* ==================== SLIDE 11 — scope limits ========================= */
{
  const s = pres.addSlide();
  kicker(s, "Scope", MUTED);
  title(s, "What this study does not claim");
  footer(s, 11, "Limits");

  const limits = [
    ["Model-derived predictions, not experimental knockouts", "Perturbation shifts are interventions on a rank-value encoding. Everything upstream is a prioritised hypothesis."],
    ["Marker-score proxy, not deconvolution", "Spot scores reflect panel expression, not inferred cell-type composition."],
    ["SCLC-specificity untested", "The validation cohort has no Normal or LUAD sections to serve as controls."],
    ["Correlational, in archival tissue", "Detection is source-state specific — absence of a result is not evidence of absence."],
  ];
  let y = 2.05;
  limits.forEach(([h, d], i) => {
    s.addShape(pres.ShapeType.ellipse, {
      x: ML + 0.1, y: y + 0.14, w: 0.3, h: 0.3,
      fill: { color: WHITE }, line: { color: i === 0 ? EOSIN : "AEBCCB", width: 2 },
    });
    s.addText(h, {
      x: ML + 0.75, y: y, w: CW - 1.2, h: 0.42, fontFace: HEAD, fontSize: 18, bold: true,
      color: i === 0 ? EOSIN : NAVY, margin: 0, valign: "middle",
    });
    s.addText(d, {
      x: ML + 0.75, y: y + 0.44, w: CW - 1.6, h: 0.44, fontFace: BODY, fontSize: 13,
      color: INK2, margin: 0, valign: "top",
    });
    y += 1.12;
  });

  s.addText("Stating these plainly is what makes everything on the previous ten slides usable.", {
    x: ML + 0.75, y: 6.5, w: CW - 1.2, h: 0.4, fontFace: BODY, fontSize: 12.5, italic: true, color: MUTED, margin: 0, valign: "top",
  });

  s.addNotes(
    "Deliver these as deliberate scope, briskly and without apology — a clean limits slide pre-empts the hostile version of every Q&A question, and this audience will respect it.\n\n" +
    "The one that matters most: perturbation shifts are model-derived interventions on a rank-value encoding, not experimental knockouts. Everything upstream is a prioritised hypothesis. Say this sentence verbatim.\n\n" +
    "The measurement caveats, quickly: detection is source-state specific, so absence of a result is not evidence of absence; the spatial work uses marker-score proxies rather than formal deconvolution; there are no Normal or LUAD sections, so SCLC-specificity of the spatial association is untested; and the design is correlational, in archival tissue.\n\n" +
    "This is a spoken slide. Hold it briefly — the audience should be listening, not reading ahead."
  );
}

/* =============== SLIDE 12 — summary & next steps ====================== */
{
  const s = pres.addSlide();
  s.background = { color: NAVY };

  s.addText("A tissue-validated hypothesis, and the CRISPR experiment that would test it", {
    x: ML, y: 0.42, w: CW - 0.2, h: 1.15, fontFace: HEAD, fontSize: 28, bold: true, color: WHITE, margin: 0, valign: "top",
  });

  // closed loop
  const ny = 1.62, nh = 0.92, gap = 0.38;
  const nw = (CW - 2 * gap) / 3;
  const loop = [
    ["Donor-clean model", "91.9%  ·  F1 0.903", TEAL],
    ["Four replicated edits", "TIM-3 · TIGIT · CTLA-4 · IL7R", HEMA_LT],
    ["Tissue-validated programme", "ρ 0.361  ·  7.95 σ", EOSIN],
  ];
  loop.forEach(([h, d, col], i) => {
    const x = ML + i * (nw + gap);
    s.addShape(pres.ShapeType.roundRect, {
      x, y: ny, w: nw, h: nh, rectRadius: 0.08,
      fill: { color: "1B3454" }, line: { color: col, width: 1.5 },
    });
    s.addText("✓", {
      x: x + 0.16, y: ny, w: 0.4, h: nh, fontFace: BODY, fontSize: 15, bold: true,
      color: col, valign: "middle", margin: 0,
    });
    s.addText(
      [
        { text: h, options: { bold: true, color: WHITE, fontSize: 13, breakLine: true } },
        { text: d, options: { color: "9DB4C8", fontSize: 11 } },
      ],
      { x: x + 0.6, y: ny, w: nw - 0.75, h: nh, fontFace: BODY, valign: "middle", margin: 0 }
    );
    if (i < 2) {
      s.addText("→", {
        x: x + nw, y: ny, w: gap, h: nh, fontFace: BODY, fontSize: 17,
        color: "6E8AA6", align: "center", valign: "middle", margin: 0,
      });
    }
  });

  // next steps
  const cw2 = (CW - 0.45) / 2;
  [["Experimental", ["CRISPR knockout of TIGIT and HAVCR2 in SCLC co-culture", "Perturb-seq in primary T cells vs in silico predictions", "Clinical outcome correlation once metadata is available"], EOSIN],
   ["Computational", ["cell2location deconvolution in place of marker scores", "Normal / LUAD sections as specificity controls", "scGPT cross-model replication"], TEAL]]
    .forEach(([h, items, col], i) => {
      const x = ML + i * (cw2 + 0.45);
      s.addShape(pres.ShapeType.roundRect, {
        x, y: 2.85, w: cw2, h: 2.0, rectRadius: 0.08,
        fill: { color: "18304E" }, line: { color: "2B4A6E", width: 1 },
      });
      s.addText(h.toUpperCase(), {
        x: x + 0.28, y: 3.0, w: cw2 - 0.56, h: 0.3, fontFace: BODY, fontSize: 11.5, bold: true,
        charSpacing: 1.4, color: col, margin: 0,
      });
      s.addText(
        items.map((t, j) => ({ text: t, options: { bullet: { indent: 14 }, breakLine: j !== items.length - 1 } })),
        { x: x + 0.28, y: 3.38, w: cw2 - 0.56, h: 1.35, fontFace: BODY, fontSize: 12.5, color: "DCE7F0", lineSpacing: 19, paraSpaceAfter: 6, margin: 0, valign: "top" }
      );
    });

  // open question hook
  s.addShape(pres.ShapeType.roundRect, {
    x: ML, y: 5.1, w: CW - 1.85, h: 0.8, rectRadius: 0.08,
    fill: { color: "3A2230" }, line: { color: EOSIN, width: 1.25 },
  });
  s.addText("Still open: why does the model order Normal < SCLC < LUAD on the exhaustion axis?", {
    x: ML + 0.28, y: 5.1, w: CW - 2.4, h: 0.8, fontFace: BODY, fontSize: 14, bold: true,
    color: "F0C9D3", valign: "middle", margin: 0,
  });

  s.addImage({ path: A("qr.png"), x: W - MR - 1.5, y: 5.1, w: 1.5, h: 1.5 });

  s.addText("github.com/Kays3/geneformer-lung-tcell  ·  snakaoka@sci.hokudai.ac.jp", {
    x: ML, y: 6.1, w: CW - 1.85, h: 0.35, fontFace: BODY, fontSize: 11.5, color: "8FB3C9", margin: 0, valign: "middle",
  });
  s.addText("Kaisar Dauyey & Shinji Nakaoka  ·  Hokkaido University  ·  24JSDP P25", {
    x: ML, y: 6.5, w: CW - 1.85, h: 0.35, fontFace: BODY, fontSize: 10.5, color: "6E8AA6", margin: 0, valign: "middle",
  });

  s.addNotes(
    "Close the loop in three beats: a donor-clean fine-tuned model (91.9%, F1 0.903) leads to four replicated, detection-adequate candidate edits, which lead to an antigen-presentation programme that tracks T-cell abundance in independent tumour tissue (rho 0.361, 7.95 sigma).\n\n" +
    "Name the falsifying experiment concretely: CRISPR knockout of TIGIT and HAVCR2 in an SCLC co-culture, and Perturb-seq in primary T cells read directly against these in silico predictions. That is what turns predicted shifts into measured ones.\n\n" +
    "Leave the open question open — the Normal < SCLC < LUAD axis. Inviting the room to argue about it is a better close than a tidy conclusion, and it is the honest state of the work."
  );
}

pres.writeFile({ fileName: path.join(__dirname, "JSDP_P25_talk.pptx") })
  .then((f) => console.log("written", f));
