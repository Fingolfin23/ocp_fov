/* Build the standalone mathematical note. Run npm install in docs first.
 * CHROME_EXECUTABLE optionally selects an installed Chromium/Chrome binary.
 * KATEX_DIST and DOCS_MODULE_ROOT are optional build-environment overrides.
 */
const fs = require('fs');
const path = require('path');
const { pathToFileURL } = require('url');
const base = __dirname;
const root = process.env.DOCS_MODULE_ROOT || path.join(base, 'node_modules');
const katexDist = process.env.KATEX_DIST || path.join(root, 'katex', 'dist');
const katex = require(path.join(katexDist, 'katex.min.js'));
const { chromium } = require(path.join(root, 'playwright'));

(async () => {
  const { marked } = await import(pathToFileURL(path.join(root, 'marked/lib/marked.esm.js')).href);
  let md = fs.readFileSync(path.join(base, 'proofs/geometric_visibility_proofs.md'), 'utf8');
  const sourceTitle = md.split('\n')[0].replace(/^# /, '');
  md = md.replace(/^# .*\n/, '');
  const mathErrors = [];
  let mathCount = 0;
  const formulas = [];
  function math(tex, displayMode) {
    mathCount++;
    try {
      return katex.renderToString(tex.trim(), {displayMode, throwOnError:true, strict:'ignore', output:'htmlAndMathml'});
    } catch (e) { mathErrors.push({tex,error:e.message}); return ''; }
  }
  md = md.replace(/\$\$([\s\S]*?)\$\$/g, (_,tex) => {
    const id=formulas.length;
    formulas.push(`<div class="equation">${math(tex,true)}</div>`);
    return `\n\nMATHBLOCKPLACEHOLDER${id}END\n\n`;
  });
  md = md.replace(/\$([^$\n]+)\$/g, (_,tex) => math(tex,false));
  let body = marked.parse(md);
  body = body.replace(/<p>MATHBLOCKPLACEHOLDER(\d+)END<\/p>/g,(_,i)=>formulas[+i]);
  if (mathErrors.length) throw new Error(JSON.stringify(mathErrors,null,2));
  body=body.replace(/(<h3>B\. Companion files<\/h3>\s*<ul>[\s\S]*?<\/ul>)/, '<div class="companions">$1</div>');
  let katexCss=fs.readFileSync(path.join(katexDist,'katex.min.css'),'utf8');
  katexCss=katexCss.replace(/src:url\(fonts\/([^)]+\.woff2)\)[^;]*;/g,(_,name)=>`src:url(data:font/woff2;base64,${fs.readFileSync(path.join(katexDist,'fonts',name)).toString('base64')}) format("woff2");`);
  body=body.replace(/src="proof_examples.png"/g,`src="data:image/png;base64,${fs.readFileSync(path.join(base,'proofs/proof_examples.png')).toString('base64')}"`);
  const css=`
  @page {size:A4; margin:21mm 20mm 21mm;}
  * {box-sizing:border-box;}
  body {margin:0;width:170mm;color:#172630;font-family:Georgia,"Times New Roman",serif;font-size:10.5pt;line-height:1.48;}
  .cover {break-after:page; padding-top:30mm;}
  .eyebrow {font-family:Arial,sans-serif;letter-spacing:2px;color:#167d8a;font-size:10pt;}
  .cover h1 {font-family:Arial,Helvetica,sans-serif;font-size:30pt;line-height:1.18;letter-spacing:0.3px;margin:15mm 0 8mm;}
  .cover .en {font-family:Georgia,serif;font-size:15pt;line-height:1.5;color:#486371;}
  .cover .rule {height:3px;background:#167d8a;margin:14mm 0 9mm;width:45mm;}
  .cover .meta {font-family:Arial,Helvetica,sans-serif;font-size:10pt;color:#50616a;}
  .cover .abstract {margin-top:14mm;padding:6mm;background:#eff7f6;border-left:3px solid #167d8a;}
  h2,h3,h4 {font-family:Arial,Helvetica,sans-serif;line-height:1.4;break-after:avoid;page-break-after:avoid;}
  h2 {font-size:15pt;color:#125e69;border-bottom:1px solid #bdd7d8;padding-bottom:3mm;margin:9mm 0 4mm;}
  h3 {font-size:11.5pt;color:#203d4a;margin:6mm 0 3mm;}
  p {margin:0 0 3.5mm;orphans:3;widows:3;}
  strong {font-family:Arial,Helvetica,sans-serif;font-weight:600;}
  ul,ol {padding-left:6mm;margin:2mm 0 4mm;}
  li {margin-bottom:2mm;break-inside:avoid;page-break-inside:avoid;}
  .equation {margin:3mm 0 4mm;break-inside:avoid;page-break-inside:avoid;}
  .katex {font-size:1.02em;}
  .katex-display {margin:0.7em 0;}
  pre {font-family:"Menlo",monospace;font-size:8.7pt;line-height:1.7;white-space:pre-wrap;padding:4mm;background:#f2f5f7;border-left:3px solid #167d8a;break-inside:avoid;}
  code {font-family:"Menlo",monospace;font-size:0.86em;overflow-wrap:anywhere;}
  table {border-collapse:collapse;width:100%;font-size:8.6pt;line-height:1.55;margin:4mm 0;}
  thead {display:table-header-group;}
  th {background:#e9f3f3;text-align:left;color:#125e69;font-family:Arial,Helvetica,sans-serif;}
  td,th {border-bottom:1px solid #d3dfe3;padding:2.2mm;vertical-align:top;}
  tr {break-inside:avoid;}
  img {max-width:100%;height:auto;}
  .companions {break-inside:avoid;page-break-inside:avoid;}
  a {color:#126777;text-decoration:none;}
  blockquote {margin:4mm 0;padding:2mm 5mm;border-left:2px solid #9ab5bd;}
  `;
  const html=`<!doctype html><html lang="en"><head><meta charset="utf-8"><title>${sourceTitle}</title><style>${katexCss}\n${css}</style></head><body>
  <section class="cover"><div class="eyebrow">OCP FOV / MATHEMATICAL NOTE</div>
  <h1>Road-Structured<br>Visibility</h1><div class="en">Formalization and Correctness Proofs<br>Common road stations and tangential events</div>
  <div class="rule"></div><p class="meta">Algorithm and source report: Yanxing Chen<br>English edition: 24 September 2026<br>Technical note for review and manuscript preparation</p>
  <div class="abstract"><strong>Two comparison stages connect local events to global visibility.</strong><br>This note formalizes tangential events, proves the classification of same-side and opposite-side intersections, characterizes the farthest visible road station, and establishes the conditions for complete visibility partitions and planning constraints.</div>
  <p style="margin-top:12mm;color:#566a73;font-size:9pt">The model is an embedded road strip with a common longitudinal coordinate and convex cross-sections. The scope of the geometric theorems, their sampled implementation, and vehicle-safety claims is stated separately.</p></section>
  <main>${body}</main></body></html>`;
  const htmlPath=path.join(base,'proofs/geometric_visibility_proofs.html');
  fs.writeFileSync(htmlPath,html);
  const launch={headless:true};
  if(process.env.CHROME_EXECUTABLE) launch.executablePath=process.env.CHROME_EXECUTABLE;
  const browser=await chromium.launch(launch);
  try {
    const page=await browser.newPage({viewport:{width:794,height:1123}});
    await page.goto(pathToFileURL(htmlPath).href);
    await page.evaluate(()=>document.fonts.ready);
    await page.emulateMedia({media:'print'});
    const overflow=await page.evaluate(()=>Array.from(document.querySelectorAll('.equation')).map((e,i)=>({i,width:e.scrollWidth,available:e.clientWidth})).filter(x=>x.width>x.available+2));
    if(overflow.length) throw new Error('Equation overflow: '+JSON.stringify(overflow));
    const pdfPath=path.join(base,'proofs/geometric_visibility_proofs.pdf');
    await page.pdf({path:pdfPath,format:'A4',printBackground:true,preferCSSPageSize:true,displayHeaderFooter:true,tagged:true,outline:true,
      headerTemplate:'<span></span>',footerTemplate:'<div style="font-family:Arial,sans-serif;font-size:8px;color:#69808c;width:100%;margin:0 20mm;display:flex;justify-content:space-between"><span>OCP FOV · MATHEMATICAL NOTE · 2026-09-24</span><span><span class="pageNumber"></span> / <span class="totalPages"></span></span></div>'});
    fs.writeFileSync(path.join(base,'proofs/build_check.json'),JSON.stringify({mathCount,mathErrors,overflow,source:'geometric_visibility_proofs.md'},null,2));
    console.log(`Rendered ${mathCount} formulas to ${pdfPath}`);
  } finally {await browser.close();}
})();
