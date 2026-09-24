# Rebuilding the documentation

The checked-in PDF, PNG, SVG and GIF files are ready to read. The steps below are only needed when changing their sources.

## Figures and animations

From the repository root:

```bash
python -m pip install -r requirements-docs.txt
python scripts/generate_idea_visuals.py --out docs/figures --animations docs/animations
```

`scripts/geometry_model.py` defines the smooth demonstration road. The generator computes its tangent events and ray intersections, checks the region against an independent line-of-sight calculation, and then writes the graphics. Its numerical check summary is `docs/figures/geometry_checks.json`. The scientific scope and each figure's intended use are described in the [Chinese visual guide](visual_guide_zh.md).

To reproduce the two analytic examples included in the proof note:

```bash
python docs/proofs/verify_examples.py
```

This refreshes `proof_examples.png` and `verification_results.json` beside the proof source. These example checks do not validate the original sampled pipeline on all road geometries.

## Full-lap Monza replay

```bash
python -m pip install -r requirements-docs.txt
python scripts/generate_monza_replay.py --track Monza.csv --out docs/animations --figures docs/figures
```

This generates a 30-second MP4, a looping README GIF, a poster and a numerical check summary. The video encoder is provided by `imageio-ffmpeg`. The original root-level GIFs and JSON cache are preserved. See [Monza replay notes](monza_replay.md) for the periodic road construction and the distinction between sampled algorithm markers and reference visibility shading.

## Mathematical note PDF

The master source is `docs/proofs/geometric_visibility_proofs.md`. The main note and its three companion Markdown notes are in English. The companion notes contain extended geometry, event and control derivations. The main PDF contains the geometric proofs and the essential planning-interface results; it is not a concatenation of all companion notes.

With a current Node.js installation satisfying the dependencies in `package.json`:

```bash
cd docs
npm install
npx playwright install chromium
npm run build:proof
```

Alternatively, use an already installed Chromium-based browser. For example, on macOS:

```bash
CHROME_EXECUTABLE='/Applications/Google Chrome.app/Contents/MacOS/Google Chrome' npm run build:proof
```

The builder typesets the formulas with KaTeX, embeds the math fonts and example image into a standalone HTML file, and prints the PDF with Chromium. It checks formula parsing and display-equation overflow, and saves `build_check.json`. The HTML file also provides a browser-readable copy without needing a math-rendering service. HTML and `build_check.json` are local build outputs ignored by Git; the English Markdown source and PDF remain checked in. The English note uses standard serif and sans-serif fonts; mathematical fonts are embedded.

After changing the mathematical content, review the formulas and rendered PDF pages. Passing a build is a formatting check, not a proof of mathematical correctness.

## Verification of this documentation edition

Edition: 2026-09-24.

- The new single-viewpoint example was exercised against the original geometry module on Monza at `s0 = 930 m`.
- The demonstration road's three cuts were checked for alignment and forward station order; 4,000 randomly sampled points agreed with an independent line-of-sight calculation.
- Both animations and all three static figures were visually inspected. The figure generator was rerun from the packaged script paths.
- The English PDF was rendered and visually inspected, with formula parsing and display-equation overflow checks. Rebuilding produces a local `docs/proofs/build_check.json` log.

The original full OCP experiments were not rerun in this documentation edition. Their external data and implementation issues are listed in [implementation notes](implementation_notes.md).

## Local outputs removed from version control

The compact `fov_cache.npz` duplicates arrays in the retained JSON and has no current consumer. The existing cache generator still writes it when run. The old `monza_930.png` snapshot can be recreated from the JSON:

```bash
MPLBACKEND=Agg python plot_fov.py --csv Monza.csv --cache fov_cache.json \
  --s0 930 --mode segment --markers --fill --show-tangent-rays \
  --save monza_930.png
```

These two outputs, the proof HTML and its build log are ignored by Git. Presentation images, vector figures, all GIFs, the MP4 and the proof PDF remain checked in.
