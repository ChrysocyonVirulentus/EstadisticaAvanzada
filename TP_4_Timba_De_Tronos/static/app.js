'use strict';
/* =====================================================================
   Timba de Tronos -- frontend (JS puro, sin dependencias).
   Toda la estadística se pide al backend (hw_core.py vía Flask).
   ===================================================================== */
const $  = (s, r = document) => r.querySelector(s);
const $$ = (s, r = document) => [...r.querySelectorAll(s)];
const NS = 'http://www.w3.org/2000/svg';
const fmt  = (x, d = 3) => (x == null || Number.isNaN(x)) ? '—' : Number(x).toFixed(d);
const fmtP = p => p == null ? '—' : p < 0.0001 ? '< 0.0001' : p.toFixed(4);
const LOCI = ['Pelo', 'Ojos'];

async function api(url, body) {
  const opt = body === undefined ? {} :
    { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify(body) };
  const r = await fetch(url, opt);
  const j = await r.json().catch(() => ({}));
  if (!r.ok) throw new Error(j.error || r.statusText);
  return j;
}

const S = { casas: [], alelos: {}, niveles: {}, costos: {}, sel: null, nSim: 1000, total: 0, rondas: 0, record: 0 };
let G = null;                                    // partida en curso

/* ---------------------------- utilidades visuales ---------------------------- */
const lum = hex => { const n = parseInt(hex.slice(1), 16);
  return (0.299 * (n >> 16 & 255) + 0.587 * (n >> 8 & 255) + 0.114 * (n & 255)) / 255; };
const chip = a => { const m = S.alelos[a] || { color: '#777', nombre: a };
  return `<span class="chip ${lum(m.color) > .6 ? 'light' : ''}" style="--c:${m.color}" title="${m.nombre}">${a}</span>`; };
const genoChips = g => `<span class="geno">${chip(g[0])}${chip(g[1])}</span>`;
const casaDe = n => S.casas.find(c => c.nombre === n);
const etiq = (a, b) => a === b ? `Casa ${a}` : `${a} × ${b}`;

/* ---------------------------- gráficos SVG ---------------------------- */
const FONT = `font-family="EB Garamond, Georgia, serif"`;

function barChart(cats, series, { w = 330, h = 210 } = {}) {
  const m = { l: 36, r: 8, t: 36, b: 30 }, iw = w - m.l - m.r, ih = h - m.t - m.b;
  const maxv = Math.max(0.0001, ...series.flatMap(s => s.vals));
  const ymax = Math.min(1, Math.max(0.2, Math.ceil(maxv * 5) / 5));
  const gw = iw / cats.length, bw = Math.min(36, gw * 0.84 / series.length);
  let s = `<svg viewBox="0 0 ${w} ${h}" role="img">`;
  for (let i = 0; i <= 4; i++) {
    const y = m.t + ih - ih * i / 4;
    s += `<line x1="${m.l}" x2="${w - m.r}" y1="${y}" y2="${y}" stroke="#5c4318" stroke-opacity="${i ? .15 : .6}"/>` +
         `<text x="${m.l - 5}" y="${y + 3.5}" text-anchor="end" font-size="10" fill="#5a4630" ${FONT}>${(ymax * i / 4).toFixed(2)}</text>`;
  }
  cats.forEach((c, ci) => {
    const cx = m.l + gw * (ci + .5);
    series.forEach((se, si) => {
      const v = se.vals[ci], bh = ih * v / ymax, x = cx - bw * series.length / 2 + si * bw, y = m.t + ih - bh;
      s += `<rect x="${x}" y="${y}" width="${bw - 1.5}" height="${Math.max(bh, 0)}" rx="1.5" fill="${se.color}" stroke="#0005"/>` +
           `<text x="${x + (bw - 1.5) / 2}" y="${y - 3}" text-anchor="middle" font-size="9.5" fill="#3a2410" ${FONT}>${v.toFixed(2)}</text>`;
    });
    s += `<text x="${cx}" y="${h - 11}" text-anchor="middle" font-size="13" font-weight="600" fill="#3a2410" ${FONT}>${c}</text>`;
  });
  series.forEach((se, i) => {
    const x = m.l + i * 105;
    s += `<rect x="${x}" y="7" width="11" height="11" fill="${se.color}" stroke="#0005"/>` +
         `<text x="${x + 15}" y="16.5" font-size="11" fill="#3a2410" ${FONT}>${se.name}</text>`;
  });
  return s + '</svg>';
}

/* Curva clásica de Hardy-Weinberg: p², 2pq y q² en función de p, con la población marcada */
function curvaHW(alelos) {
  const [A, pA] = alelos[0], B = alelos[1] ? alelos[1][0] : null;
  const w = 330, h = 210, m = { l: 36, r: 12, t: 28, b: 34 }, iw = w - m.l - m.r, ih = h - m.t - m.b;
  const X = p => m.l + iw * p, Y = f => m.t + ih - ih * f;
  const curvas = [
    { f: p => p * p,         c: '#b8921f', n: `${A}${A}` },
    { f: p => 2 * p * (1 - p), c: '#8b1a1a', n: B ? `${A}${B}` : '2pq' },
    { f: p => (1 - p) ** 2,  c: '#4f7aa3', n: B ? `${B}${B}` : 'q²' },
  ];
  let s = `<svg viewBox="0 0 ${w} ${h}" role="img">`;
  for (let i = 0; i <= 4; i++)
    s += `<line x1="${m.l}" x2="${w - m.r}" y1="${Y(i / 4)}" y2="${Y(i / 4)}" stroke="#5c4318" stroke-opacity="${i ? .15 : .6}"/>` +
         `<text x="${m.l - 5}" y="${Y(i / 4) + 3.5}" text-anchor="end" font-size="10" fill="#5a4630" ${FONT}>${(i / 4).toFixed(2)}</text>`;
  for (let i = 0; i <= 4; i++)
    s += `<text x="${X(i / 4)}" y="${h - 17}" text-anchor="middle" font-size="10" fill="#5a4630" ${FONT}>${(i / 4).toFixed(2)}</text>`;
  s += `<text x="${m.l + iw / 2}" y="${h - 3}" text-anchor="middle" font-size="11" fill="#3a2410" ${FONT}>frecuencia del alelo ${A}  (p)</text>`;
  curvas.forEach(cv => {
    const pts = Array.from({ length: 41 }, (_, i) => `${X(i / 40).toFixed(1)},${Y(cv.f(i / 40)).toFixed(1)}`).join(' ');
    s += `<polyline points="${pts}" fill="none" stroke="${cv.c}" stroke-width="2.2"/>`;
  });
  s += `<line x1="${X(pA)}" x2="${X(pA)}" y1="${Y(0)}" y2="${Y(1)}" stroke="#2a1d10" stroke-dasharray="4 3"/>`;
  curvas.forEach(cv => { s += `<circle cx="${X(pA)}" cy="${Y(cv.f(pA))}" r="4.5" fill="${cv.c}" stroke="#fff" stroke-width="1.5"/>`; });
  curvas.forEach((cv, i) => {
    const x = m.l + i * 70;
    s += `<line x1="${x}" x2="${x + 14}" y1="12" y2="12" stroke="${cv.c}" stroke-width="3"/>` +
         `<text x="${x + 18}" y="16" font-size="11" fill="#3a2410" ${FONT}>${cv.n}</text>`;
  });
  s += `<text x="${w - m.r}" y="16" text-anchor="end" font-size="11" font-style="italic" fill="#3a2410" ${FONT}>p = ${pA.toFixed(2)}</text>`;
  return s + '</svg>';
}

/* ---------------------------- pestañas ---------------------------- */
function setupTabs() {
  $$('.tab').forEach(t => t.addEventListener('click', () => {
    $$('.tab').forEach(x => x.classList.toggle('active', x === t));
    $$('.view').forEach(v => v.classList.toggle('active', v.id === 'view-' + t.dataset.tab));
  }));
}
function setScore() {
  $('#score').textContent = S.total; $('#rondas').textContent = S.rondas; $('#record').textContent = S.record;
  try { localStorage.setItem('hw_record', S.record); } catch (e) { /* sin storage: no pasa nada */ }
}

/* =====================================================================
   ATLAS
   ===================================================================== */
function buildPins() {
  const cont = $('#pins');
  S.casas.forEach(c => {
    const g = document.createElementNS(NS, 'g');
    g.setAttribute('class', 'pin'); g.setAttribute('tabindex', '0'); g.setAttribute('role', 'button');
    g.setAttribute('data-n', c.nombre); g.setAttribute('aria-label', `${c.titulo}, ${c.lugar}`);
    g.setAttribute('transform', `translate(${c.pos.x},${c.pos.y})`);
    g.innerHTML = `<g class="in"><circle class="halo" r="28"/>
        <circle class="disc" r="21" fill="${c.color}" stroke="${c.color2}"/>
        <text text-anchor="middle" y="8" font-size="22">${c.emblema}</text></g>
      <text class="lbl" y="44" text-anchor="middle">${c.nombre}</text>
      <text class="sub" y="57" text-anchor="middle">${c.lugar}</text>`;
    g.addEventListener('click', () => mostrarCasa(c.nombre));
    g.addEventListener('keydown', e => { if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); mostrarCasa(c.nombre); } });
    cont.appendChild(g);
  });
}

function ecuacion(loc) {
  const sim = loc.genotipos.map(([g]) => g[0] === g[1] ? `${g[0]}²` : `2${g}`).join(' + ') + ' = 1';
  const num = loc.genotipos.map(([, p]) => p.toFixed(2)).join(' + ') + ' = 1';
  return `${sim}<br><span class="num">${num}</span>`;
}

function locusHTML(c, key, sim, n) {
  const loc = c.loci[key], sl = sim.loci[key], obsAl = Object.fromEntries(sl.alelos_obs);
  const cats = sl.genotipos.map(x => x.g);
  const bars = barChart(cats, [
    { name: 'Observado', color: '#7f9db9', vals: sl.genotipos.map(x => x.obs / n) },
    { name: 'Esperado',  color: '#c9a227', vals: sl.genotipos.map(x => x.p) }]);
  const fija = loc.alelos.length === 1;
  const t = sl.chi2;
  const chi = t.gl === 0
    ? `<span class="badge na">LOCUS FIJO</span> un solo alelo (frecuencia 1): no hay variación ni test posible.`
    : `χ² = <b>${fmt(t.chi2, 3)}</b> · gl = ${t.gl} · p = <b>${fmtP(t.p)}</b>
       <span class="badge ${t.p > .05 ? 'ok' : 'rej'}">${t.p > .05 ? 'AJUSTA A H-W' : 'NO AJUSTA'}</span>`;
  return `<section class="locus"><h3>Locus ${key}</h3>
    <div class="alelos">${loc.alelos.map(([a, f]) =>
      `<span class="alelo-box">${chip(a)}<span>${S.alelos[a] ? S.alelos[a].nombre : ''} · <b>p = ${f.toFixed(3)}</b>
       <small>(obs. ${fmt(obsAl[a])})</small></span></span>`).join('')}</div>
    <div class="eq">${ecuacion(loc)}</div>
    <div class="charts">
      <div class="chart">${bars}<h4>GENOTIPOS · N = ${n}</h4></div>
      <div class="chart">${curvaHW(loc.alelos)}<h4>${fija ? 'LOCUS FIJO (p = 1)' : 'CURVA DE HARDY-WEINBERG'}</h4></div>
    </div>
    <table class="tbl"><tr><th>Genotipo</th><th>Esperado</th><th>Esp. (n)</th><th>Obs. (n)</th><th>Obs.</th></tr>
      ${sl.genotipos.map(x => `<tr><td>${genoChips(x.g)}</td><td>${x.p.toFixed(3)}</td><td>${x.esp.toFixed(1)}</td>
        <td>${x.obs}</td><td>${(x.obs / n).toFixed(3)}</td></tr>`).join('')}</table>
    <p class="chi">${chi}</p></section>`;
}

async function mostrarCasa(nombre, n) {
  S.sel = nombre; n = n || S.nSim; S.nSim = n;
  $$('.pin').forEach(p => p.classList.toggle('sel', p.dataset.n === nombre));
  const c = casaDe(nombre), panel = $('#panel-casa');
  let sim;
  try { sim = await api(`/api/casa/${nombre}/simular?n=${n}`); }
  catch (e) { panel.innerHTML = `<div class="vacio">⚠ ${e.message}</div>`; return; }
  if (S.sel !== nombre) return;                      // el usuario ya clickeó otra casa
  const ind = sim.independencia;
  panel.innerHTML = `
    <div class="casa-head" style="--c1:${c.color};--c2:${c.color2}">
      <div class="escudo">${c.emblema}</div>
      <div><h2>${c.titulo}</h2><p class="lugar">📍 ${c.lugar}</p></div></div>
    <div class="sim-ctrl">Población simulada:
      <select id="n-sim">${[50, 200, 1000, 10000].map(v => `<option ${v === n ? 'selected' : ''}>${v}</option>`).join('')}</select> individuos
      <button class="btn" id="re-sim">🎲 Simular otra muestra</button></div>
    ${LOCI.map(k => locusHTML(c, k, sim, n)).join('')}
    <p class="chi"><b>Independencia entre loci:</b> ${ind.gl === 0
      ? 'no aplica (algún locus es fijo).'
      : `χ² = ${fmt(ind.chi2, 3)} · gl = ${ind.gl} · p = ${fmtP(ind.p)} — ${ind.veredicto}`}</p>`;
  $('#n-sim').addEventListener('change', e => mostrarCasa(nombre, +e.target.value));
  $('#re-sim').addEventListener('click', () => mostrarCasa(nombre, n));
}

/* =====================================================================
   JUEGO
   ===================================================================== */
function renderIntro() {
  const v = $('#view-juego');
  v.innerHTML = `<div class="card intro"><h2 style="font-size:1.7rem">⚔ Timba de Tronos</h2>
    <p>Un maestre te muestra algunos habitantes de Poniente. Todos descienden de un <b>linaje secreto</b>:
    o bien de <b>una sola casa</b> que se casó entre sus propios miembros (endogamia), o bien del
    <b>cruce de dos casas</b> distintas.</p>
    <p>Contá los alelos de cada locus, calculá las frecuencias, ponelas a prueba con el <b>χ²</b> contra cada hipótesis
    y decidí. Más evidencia ayuda… pero cuesta puntos.</p>
    <p class="hint">Consejo: estudiá primero el Atlas. Ahí está la composición de alelos de cada casa.</p>
    <div class="niveles">${Object.entries(S.niveles).map(([k, n]) => `
      <div class="nivel" data-nivel="${k}"><b>${n.titulo}</b>
        <span>${n.n_inicial} individuos iniciales${n.blackfyre ? '<br>Incluye la casa Blackfyre' : '<br>Sin Blackfyre'}</span></div>`).join('')}</div>
  </div>`;
  $$('.nivel', v).forEach(el => el.addEventListener('click', () => nuevaPartida(el.dataset.nivel)));
}

async function nuevaPartida(nivel) {
  try {
    const st = await api('/api/juego/nuevo', { nivel });
    G = { id: st.id, nivel, individuos: st.individuos, extras: 0, tests: 0, puntos: st.puntos, modo: 'casa',
          a: S.casas[0].nombre, b: S.casas[1].nombre, historial: [], chi: null, calc: {}, calcShow: false,
          nuevo: -1, over: false };
    renderJuego();
  } catch (e) { alert(e.message); }
}

function renderJuego() {
  const v = $('#view-juego'), c = S.costos;
  v.innerHTML = `
    <div class="card" style="max-width:1400px;margin:0 auto 1.2rem;padding:.6rem 1.2rem;display:flex;gap:1.2rem;align-items:center;flex-wrap:wrap">
      <b style="font-family:var(--title)">Rango: ${S.niveles[G.nivel].titulo}</b>
      <span>Puntos en juego: <b id="pts" style="font-size:1.3rem;color:var(--crimson)"></b></span>
      <span class="hint" style="margin:0">Evidencia extra −${c.evidencia} · consulta χ² nueva −${c.test} · acertar vale hasta ${c.base}</span>
      <button class="btn ghost" id="abandonar" style="margin-left:auto">↺ Nueva partida</button></div>
    <div class="game-grid">
      <div class="col">
        <section class="card"><h2>👥 Los sospechosos <small>descendientes de un linaje secreto</small></h2>
          <div class="personas" id="personas"></div>
          <div class="fila-btns"><button class="btn" id="btn-evi">🔍 Pedir más evidencia <small>(−${c.evidencia})</small></button>
            <span class="hint" id="evi-info"></span></div></section>
        <section class="card"><h2>📜 Pergamino de cálculo <small>frecuencias alélicas de la muestra</small></h2>
          <p class="hint" style="margin-top:0">p̂(A) = (2·#AA + #Aa) / (2n). Escribilo como decimal o fracción (3/8).</p>
          <div id="calc" class="calc"></div>
          <div class="fila-btns"><button class="btn" id="btn-check">✔ Comprobar</button>
            <button class="btn ghost" id="btn-show">👁 Mostrar solución</button></div></section>
      </div>
      <div class="col">
        <section class="card"><h2>⚖ Tu hipótesis <small>H₀ a poner a prueba</small></h2>
          <div class="seg"><button id="m-casa">Una sola casa (endogamia)</button><button id="m-cruce">Cruce de dos casas</button></div>
          <div class="selrow" id="selrow"></div>
          <div class="fila-btns">
            <button class="btn" id="btn-chi">⚖ Consultar χ² <small>(−${c.test})</small></button>
            <button class="btn red" id="btn-ver">👑 Proclamar veredicto</button></div>
          <p class="hint" id="hip-hint"></p>
          <div id="chi-res"></div><div class="hist" id="hist"></div></section>
        <section class="card libro"><details open><summary>📖 Libro de las Casas</summary><div id="libro"></div></details></section>
      </div>
    </div>`;
  $('#abandonar').addEventListener('click', renderIntro);
  $('#btn-evi').addEventListener('click', pedirEvidencia);
  $('#btn-check').addEventListener('click', () => comprobar(false));
  $('#btn-show').addEventListener('click', () => comprobar(true));
  $('#m-casa').addEventListener('click', () => { G.modo = 'casa'; renderHip(); });
  $('#m-cruce').addEventListener('click', () => { G.modo = 'cruce'; if (G.a === G.b) G.b = S.casas.find(x => x.nombre !== G.a).nombre; renderHip(); });
  $('#btn-chi').addEventListener('click', consultarChi);
  $('#btn-ver').addEventListener('click', proclamar);
  renderLibro(); renderTodo();
}

function renderTodo() { renderPersonas(); renderCalc(); renderHip(); renderChi(); renderHist(); renderPuntos(); }
function renderPuntos() { const e = $('#pts'); if (e) e.textContent = G.puntos; }

function renderPersonas() {
  $('#personas').innerHTML = G.individuos.map((p, i) => `
    <div class="persona ${i === G.nuevo ? 'nueva' : ''}"><div class="n">PERSONA ${i + 1}</div>
      <div class="row"><i>Pelo</i>${genoChips(p.pelo)}</div><div class="row"><i>Ojos</i>${genoChips(p.ojos)}</div></div>`).join('');
  const quedan = S.costos.max_extra - G.extras;
  $('#btn-evi').disabled = G.over || quedan <= 0;
  $('#evi-info').textContent = quedan > 0 ? `n = ${G.individuos.length} · quedan ${quedan} evidencias`
                                          : `n = ${G.individuos.length} · no hay más evidencia`;
}

async function pedirEvidencia() {
  try {
    const st = await api(`/api/juego/${G.id}/evidencia`, {});
    G.individuos = st.individuos; G.extras = st.extras; G.puntos = st.puntos; G.nuevo = G.individuos.length - 1;
    G.calc = {}; G.calcShow = false; G.chi = null;                 // los conteos cambiaron
    renderTodo();
  } catch (e) { alert(e.message); }
}

/* ---- pergamino de cálculo ---- */
function copiasAlelos(key) {
  const cnt = {};
  G.individuos.forEach(i => { for (const a of i[key]) cnt[a] = (cnt[a] || 0) + 1; });
  return cnt;
}
const parseNum = s => { s = String(s).trim().replace(',', '.'); if (!s) return NaN;
  if (s.includes('/')) { const [a, b] = s.split('/').map(Number); return b ? a / b : NaN; } return Number(s); };

function renderCalc() {
  const N2 = 2 * G.individuos.length;
  $('#calc').innerHTML = LOCI.map(loc => {
    const cnt = copiasAlelos(loc.toLowerCase());
    const filas = Object.entries(cnt).sort((x, y) => y[1] - x[1]).map(([a, k]) => {
      const st = G.calc[loc + a] || {};
      return `<tr><td>${chip(a)} <small>${S.alelos[a] ? S.alelos[a].nombre : ''}</small></td>
        <td>${G.calcShow ? `${k} / ${N2}` : '<span class="hint">?</span>'}</td>
        <td><input data-k="${loc}${a}" class="${st.est || ''}" value="${st.val ?? ''}" inputmode="decimal" placeholder="p̂"></td>
        <td class="est">${st.est === 'ok' ? '✔' : st.est === 'mal' ? '✘' : ''}</td></tr>`;
    }).join('');
    return `<div class="calc-locus"><h3>Locus ${loc}</h3><table class="tbl"><tr><th>Alelo</th><th>Copias</th><th>Tu p̂</th><th></th></tr>${filas}</table></div>`;
  }).join('');
  $$('#calc input').forEach(inp => inp.addEventListener('input', () => {
    G.calc[inp.dataset.k] = { val: inp.value }; inp.className = ''; inp.closest('tr').querySelector('.est').textContent = ''; }));
}

function comprobar(mostrar) {
  const N2 = 2 * G.individuos.length;
  if (mostrar) G.calcShow = true;
  LOCI.forEach(loc => {
    const cnt = copiasAlelos(loc.toLowerCase());
    Object.entries(cnt).forEach(([a, k]) => {
      const real = k / N2, key = loc + a, cur = G.calc[key] || {};
      if (mostrar) { G.calc[key] = { val: real.toFixed(3), est: 'ok' }; return; }
      if (!String(cur.val ?? '').trim()) { G.calc[key] = { val: cur.val }; return; }
      G.calc[key] = { val: cur.val, est: Math.abs(parseNum(cur.val) - real) <= 0.006 ? 'ok' : 'mal' };
    });
  });
  renderCalc();
}

/* ---- hipótesis ---- */
function renderHip() {
  const opts = sel => S.casas.map(c => `<option value="${c.nombre}" ${c.nombre === sel ? 'selected' : ''}>${c.emblema} ${c.nombre}</option>`).join('');
  $('#m-casa').classList.toggle('on', G.modo === 'casa');
  $('#m-cruce').classList.toggle('on', G.modo === 'cruce');
  $('#selrow').innerHTML = G.modo === 'casa'
    ? `<span>Todos descienden de la casa</span><select id="sel-a">${opts(G.a)}</select>`
    : `<select id="sel-a">${opts(G.a)}</select><span class="x">×</span><select id="sel-b">${opts(G.b)}</select>`;
  $('#sel-a').addEventListener('change', e => { G.a = e.target.value; validarHip(); });
  if ($('#sel-b')) $('#sel-b').addEventListener('change', e => { G.b = e.target.value; validarHip(); });
  validarHip();
}
const hipotesis = () => G.modo === 'casa' ? [G.a, G.a] : [G.a, G.b];
function validarHip() {
  const bad = G.modo === 'cruce' && G.a === G.b;
  $('#btn-chi').disabled = $('#btn-ver').disabled = bad || G.over;
  $('#hip-hint').textContent = bad ? 'Elegí dos casas distintas para un cruce.' :
    `H₀: la muestra proviene de la F1 de «${etiq(...hipotesis())}».`;
}

async function consultarChi() {
  const [a, b] = hipotesis();
  try {
    const t = await api(`/api/juego/${G.id}/chi2`, { a, b });
    G.puntos = t.puntos; G.chi = t;
    G.historial = G.historial.filter(h => !(h.a === t.a && h.b === t.b && h.n === t.n));
    G.historial.unshift({ a: t.a, b: t.b, n: t.n, p: t.p, rechaza: t.rechaza, imposible: t.imposible });
    renderChi(); renderHist(); renderPuntos();
  } catch (e) { alert(e.message); }
}

function renderChi() {
  const t = G.chi, box = $('#chi-res');
  if (!t) { box.innerHTML = ''; return; }
  const tablas = LOCI.map(loc => {
    const L = t.loci[loc];
    return `<h4>Locus ${loc} — χ² = ${L.imposibles.length ? '∞' : fmt(L.chi2, 2)} · gl = ${L.gl}</h4>
      <table class="tbl"><tr><th>Genotipo</th><th>Obs.</th><th>Esp.</th><th>(O−E)²/E</th></tr>
      ${L.filas.map(f => `<tr class="${f.contrib == null ? 'bad' : ''}"><td>${genoChips(f.g)}</td><td>${f.obs}</td>
        <td>${f.esp.toFixed(2)}</td><td>${f.contrib == null ? 'imposible' : f.contrib.toFixed(2)}</td></tr>`).join('')}</table>`;
  }).join('');
  const badge = t.imposible || t.rechaza ? 'rej' : (t.gl === 0 ? 'na' : 'ok');
  box.innerHTML = `<div class="chi-card ${t.rechaza ? 'rej' : 'ok'}">
    <h3>H₀: «${t.etiqueta}» · n = ${t.n}</h3>${tablas}
    <div class="resumen"><span>χ² = <b>${t.imposible ? '∞' : fmt(t.chi2, 2)}</b></span><span>gl = <b>${t.gl}</b></span>
      <span>p = <b>${t.imposible ? '0' : fmtP(t.p)}</b></span><span class="badge ${badge}">${t.veredicto.toUpperCase()}</span></div>
    ${t.esperados_bajos && !t.imposible ? `<div class="warn">⚠ Hay frecuencias esperadas &lt; 5 (regla de Cochran): con muestras tan chicas el χ²
      es solo aproximado y rara vez rechaza. Más evidencia lo vuelve más confiable.</div>` : ''}</div>`;
}

function renderHist() {
  const h = $('#hist');
  h.innerHTML = G.historial.length ? `<h4>CONSULTAS REALIZADAS</h4>` + G.historial.map((x, i) => `
    <button data-i="${i}"><span>${etiq(x.a, x.b)} <small>(n=${x.n})</small></span>
      <span>${x.imposible ? 'imposible' : 'p = ' + fmtP(x.p)} ${x.rechaza ? '✘' : '✔'}</span></button>`).join('') : '';
  $$('button', h).forEach(b => b.addEventListener('click', () => {
    const x = G.historial[+b.dataset.i]; G.modo = x.a === x.b ? 'casa' : 'cruce'; G.a = x.a; G.b = x.b; renderHip(); }));
}

function renderLibro() {
  $('#libro').innerHTML = `<table class="tbl"><tr><th>Casa</th><th>Pelo</th><th>Ojos</th></tr>` + S.casas.map(c => {
    const f = k => c.loci[k].alelos.map(([a, p]) => `${chip(a)}&thinsp;${p.toFixed(2)}`).join(' &nbsp;');
    return `<tr><td>${c.emblema} <b>${c.nombre}</b></td><td style="text-align:left">${f('Pelo')}</td><td style="text-align:left">${f('Ojos')}</td></tr>`;
  }).join('') + `</table><p class="hint">Frecuencias alélicas de cada casa (las mismas del Atlas).
    Una casa sola está en H-W (p², 2pq, q²); un cruce A×B da una F1 con el promedio de ambas frecuencias.</p>`;
}

/* ---- veredicto ---- */
async function proclamar() {
  const [a, b] = hipotesis();
  if (!confirm(`¿Proclamar que el linaje es «${etiq(a, b)}»?\nNo hay vuelta atrás.`)) return;
  try {
    const r = await api(`/api/juego/${G.id}/responder`, { a, b });
    G.over = true; S.rondas++; S.total += r.puntos; S.record = Math.max(S.record, S.total); setScore();
    $('#btn-evi').disabled = $('#btn-chi').disabled = $('#btn-ver').disabled = true;
    mostrarResultado(r);
  } catch (e) { alert(e.message); }
}

function mostrarResultado(r) {
  const win = r.acierto;
  const chart = loc => {
    const rows = r.f1[loc];
    return `<div class="chart">${barChart(rows.map(x => x.g), [
      { name: `F1 simulada (${r.pool})`, color: '#7f9db9', vals: rows.map(x => x.obs) },
      { name: 'Esperado', color: '#c9a227', vals: rows.map(x => x.esp) }], { w: 330, h: 200 })}<h4>F1 SECRETA · ${loc.toUpperCase()}</h4></div>`;
  };
  $('#modal-body').innerHTML = `
    <h2 class="veredicto ${win ? 'win' : 'lose'}">${win ? '👑 ¡Veredicto acertado!' : '💀 El linaje te engañó'}</h2>
    <div class="reveal-grid">
      <div><p>Linaje secreto:<br><b style="font-size:1.4rem">${r.secreto.tipo === 'casa' ? '🏰' : '💍'} ${r.secreto.etiqueta}</b></p>
        <p>Tu respuesta: <b>${r.respuesta.etiqueta}</b>${win && r.respuesta.etiqueta !== r.secreto.etiqueta
          ? '<br><i>(estadísticamente indistinguible del secreto: misma F1 ✔)</i>' : ''}</p>
        <p class="hint">${r.n} individuos observados · ${r.extras} evidencias extra · ${r.tests} consultas χ²</p></div>
      <div><div class="cifra">+${r.puntos}</div><p style="text-align:center" class="hint">puntos · total ${S.total}</p></div></div>
    <div class="charts">${chart('Pelo')}${chart('Ojos')}</div>
    <h3 style="margin:1rem 0 .3rem;font-size:.95rem">Ranking de hipótesis contra TODA la evidencia (por p-valor)</h3>
    <table class="tbl"><tr><th>Hipótesis</th><th>p-valor</th><th></th></tr>${r.ranking.map(x => `
      <tr class="${x.imposible ? 'bad' : ''}"><td>${x.es_verdad ? '👑 ' : ''}${x.etiqueta}${x.a === r.respuesta.a && x.b === r.respuesta.b ? ' <small>(tu respuesta)</small>' : ''}</td>
      <td>${x.imposible ? 'imposible' : fmtP(x.p)}</td><td>${x.rechaza ? '✘ rechazada' : '✔ compatible'}</td></tr>`).join('')}</table>
    <div class="fila-btns" style="margin-top:1rem;justify-content:center">
      <button class="btn" id="m-otra">⚔ Nueva ronda</button><button class="btn ghost" id="m-cerrar">Cerrar</button></div>`;
  $('#modal').classList.remove('hidden');
  $('#m-cerrar').addEventListener('click', () => $('#modal').classList.add('hidden'));
  $('#m-otra').addEventListener('click', () => { $('#modal').classList.add('hidden'); nuevaPartida(G.nivel); });
}

/* ---------------------------- arranque ---------------------------- */
(async function init() {
  setupTabs();
  try { S.record = +localStorage.getItem('hw_record') || 0; } catch (e) { S.record = 0; }
  try {
    const d = await api('/api/casas');
    Object.assign(S, { casas: d.casas, alelos: d.alelos, niveles: d.niveles, costos: d.costos });
    buildPins(); renderIntro(); setScore();
  } catch (e) {
    $('#panel-casa').innerHTML = `<div class="vacio">⚠ No se pudo contactar al servidor: ${e.message}</div>`;
  }
})();
