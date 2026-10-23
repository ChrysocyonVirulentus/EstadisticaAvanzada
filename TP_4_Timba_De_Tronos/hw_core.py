"""
hw_core.py -- Motor estadístico del simulador Hardy-Weinberg de Poniente.

Las clases y funciones de las secciones 2-10 vienen TAL CUAL del notebook
`Simulador_Hardy_Weinberg_Poniente.ipynb` (sin matplotlib, que acá no hace
falta: los gráficos los dibuja el frontend). Todo lo que NO está en el
notebook está en la sección "EXTENSIONES PARA EL JUEGO" al final del archivo.
"""
import math
import numpy as np
from abc import ABC, abstractmethod


class LocusHardyWeinberg(ABC):
    """
    Clase abstracta que representa UN locus génico en equilibrio
    Hardy-Weinberg. Define la estructura común: frecuencias alélicas,
    probabilidades genotípicas y simulación multinomial.

    Las subclases concretas deben implementar _calcular_probabilidades(),
    que define cómo se obtienen las probabilidades genotípicas a partir
    de las frecuencias alélicas (la ley de H-W: p², 2pq, q², ...).
    """

    def __init__(self, nombre: str, frecuencias: dict, seed: int = None):
        """
        Constructor base (lo heredan todas las subclases).

        Parámetros:
            nombre (str): nombre del locus, p. ej. "Pelo" u "Ojos".
            frecuencias (dict): {alelo: frecuencia}, debe sumar 1.
            seed (int, opcional): semilla para reproducibilidad.
        """
        # --- Validaciones ---
        if not frecuencias:
            raise ValueError("El locus necesita al menos un alelo.")
        if any(f < 0 for f in frecuencias.values()):
            raise ValueError("Las frecuencias no pueden ser negativas.")
        if abs(sum(frecuencias.values()) - 1.0) > 1e-9:
            raise ValueError("Las frecuencias alélicas deben sumar 1.")

        # --- Atributos (propiedades de la instancia) ---
        self.nombre = nombre
        self.frecuencias = dict(frecuencias)   # frecuencias alélicas (teóricas)
        self.seed = seed

        # Generador aleatorio propio de cada instancia (buena práctica OOP)
        self.rng = np.random.default_rng(seed)

        # Se calculan al construir la instancia (propiedad derivada)
        self.prob_genotipicas = self._calcular_probabilidades()

    @property
    def alelos(self) -> list:
        """Lista de alelos del locus."""
        return list(self.frecuencias)

    @abstractmethod
    def _calcular_probabilidades(self) -> dict:
        """
        Método abstracto: cada subclase define aquí la ley de H-W
        que le corresponde. Debe devolver {genotipo: probabilidad}.
        """

    def simular(self, n: int) -> dict:
        """
        Simula n individuos para ESTE locus con una multinomial
        sobre los genotipos posibles.
        """
        genotipos = list(self.prob_genotipicas)
        probs = list(self.prob_genotipicas.values())
        conteos = self.rng.multinomial(n, probs)

        return {
            "conteos": dict(zip(genotipos, conteos)),
            "frec_obs": dict(zip(genotipos, conteos / n)),
            "frec_esp": self.prob_genotipicas,
        }

    def frecuencias_alelicas_observadas(self, conteos: dict) -> dict:
        """
        Reestima las frecuencias alélicas desde los conteos genotípicos:
        cada homocigoto aporta 2 copias, cada heterocigoto 1 por alelo.
        """
        total = sum(conteos.values()) * 2          # copias totales de alelos
        obs = {a: 0.0 for a in self.alelos}
        for gen, c in conteos.items():
            obs[gen[0]] += c
            obs[gen[1]] += c
        return {a: obs[a] / total for a in self.alelos}

class LocusEquilibrio(LocusHardyWeinberg):
    """
    Locus concreto con la ley de Hardy-Weinberg generalizada a k alelos:

        P(AiAi)   = pi²          (homocigotos)
        P(AiAj)   = 2 * pi * pj  (heterocigotos, i != j)

    Con 2 alelos se reduce exactamente a p² + 2pq + q² = 1.
    """

    def _calcular_probabilidades(self) -> dict:
        probs = {}
        alelos = self.alelos
        for i, a1 in enumerate(alelos):
            for j in range(i, len(alelos)):
                a2 = alelos[j]
                p = self.frecuencias[a1] * self.frecuencias[a2]
                p *= 2 if i != j else 1          # factor 2 del heterocigoto
                # Etiqueta del genotipo: alelo dominante (mayúscula) primero
                gen = a1 + a2 if i == j else "".join(
                    sorted([a1, a2], key=lambda x: (not x.isupper(), x)))
                probs[gen] = p
        return probs

class DosLocus:
    """
    Combina dos loci en equilibrio Hardy-Weinberg bajo la hipótesis
    de independencia (no hay ligamiento génico).

    La probabilidad conjunta es el producto de las marginales:
        P(geno_pelo, geno_ojos) = P(geno_pelo) * P(geno_ojos)

    y la simulación de n individuos es una multinomial sobre las
    k1 * k2 combinaciones genotípicas posibles.
    """

    def __init__(self, locus_1: LocusHardyWeinberg, locus_2: LocusHardyWeinberg,
                 seed: int = None):
        self.locus_1 = locus_1
        self.locus_2 = locus_2
        self.rng = np.random.default_rng(seed)

        # Tabla de probabilidades conjuntas (producto de marginales)
        self.prob_conjunta = {
            f"{g1}-{g2}": p1 * p2
            for g1, p1 in locus_1.prob_genotipicas.items()
            for g2, p2 in locus_2.prob_genotipicas.items()
        }

    def simular(self, n: int) -> dict:
        """Simula n individuos a dos loci con una multinomial conjunta."""
        combos = list(self.prob_conjunta)
        probs = list(self.prob_conjunta.values())
        conteos_vec = self.rng.multinomial(n, probs)
        conteos = dict(zip(combos, conteos_vec))

        # Marginales por locus (sumando sobre el otro locus)
        def marginal(locus, idx):
            marg = {g: 0 for g in locus.prob_genotipicas}
            for combo, c in conteos.items():
                marg[combo.split("-")[idx]] += c
            return marg

        return {
            "conteos_conjunta": conteos,
            "marginal_locus_1": marginal(self.locus_1, 0),
            "marginal_locus_2": marginal(self.locus_2, 1),
        }

class Casa:
    """
    Casa de Poniente definida por las frecuencias alélicas de sus dos
    loci fenotípicos (pelo y ojos). Al instanciarla se generan los
    dos loci y el simulador de dos loci, listos para usar.

    Uso:
        stark = Casa("Stark", frec_pelo={"M": 0.7, "r": 0.3},
                     frec_ojos={"C": 1.0}, seed=42)
        stark.simular(10_000)
    """

    def __init__(self, nombre: str, frec_pelo: dict, frec_ojos: dict,
                 seed: int = None):
        self.nombre = nombre
        self.locus_pelo = LocusEquilibrio("Pelo", frec_pelo, seed=seed)
        self.locus_ojos = LocusEquilibrio("Ojos", frec_ojos, seed=seed)
        self.genoma = DosLocus(self.locus_pelo, self.locus_ojos, seed=seed)

    def simular(self, n: int = 10_000) -> dict:
        """Simula n personas de la casa y devuelve el resultado a dos loci."""
        return self.genoma.simular(n)

    def mostrar(self, n: int = 10_000, top: int = 6) -> None:
        """Muestra frecuencias alélicas y los genotipos más frecuentes."""
        res = self.simular(n)
        print("=" * 60)
        print(f"  CASA {self.nombre.upper()}  (N = {n})")
        print("=" * 60)
        for locus, marg in ((self.locus_pelo, res["marginal_locus_1"]),
                            (self.locus_ojos, res["marginal_locus_2"])):
            obs = locus.frecuencias_alelicas_observadas(marg)
            teor = ", ".join(f"{a}: {f:.3f}" for a, f in locus.frecuencias.items())
            obss = ", ".join(f"{a}: {f:.3f}" for a, f in obs.items())
            print(f"  Locus {locus.nombre}: teórico {{{teor}}} | observado {{{obss}}}")
        ranking = sorted(res["conteos_conjunta"].items(),
                         key=lambda kv: -kv[1])[:top]
        print("  Genotipos (pelo-ojos) más frecuentes:")
        for gen, c in ranking:
            print(f"    {gen:8s}: {c:>6} ({c/n:.2%})")
        print()

def crear_casas(seed: int = 42) -> dict:
    """
    Crea las instancias de cada Casa según las especificaciones.

    Nota: las ecuaciones tipo M² + 2Mr + r² = 1 describen el EQUILIBRIO
    pero no fijan las frecuencias: hay que elegir p y q (ajustables
    en cada entrada del diccionario).
    """
    specs = {
        "Stark": Casa("Stark", {"M": 0.7, "r": 0.3}, {"C": 1.0}, seed=seed),
        "Lannister": Casa("Lannister", {"n": 1.0}, {"a": 1.0}, seed=seed),
        "Targaryen": Casa("Targaryen", {"b": 0.6, "N": 0.4}, {"t": 1.0}, seed=seed),
        "Martell": Casa("Martell", {"v": 0.5, "M": 0.5}, {"C": 1.0}, seed=seed),
        "Blackfyre": Casa("Blackfyre", {"b": 0.5, "N": 0.5},
                          {"t": 0.5, "v": 0.5}, seed=seed),
        "Baratheon": Casa("Baratheon", {"N": 1.0}, {"v": 1.0}, seed=seed),
        "Tyrell": Casa("Tyrell", {"n": 0.5, "N": 0.5}, {"a": 1.0}, seed=seed),
        "Hightower": Casa("Hightower", {"r": 0.4, "M": 0.6}, {"C": 1.0}, seed=seed),
    }
    return specs

def _etiqueta_genotipo(a1: str, a2: str) -> str:
    """Etiqueta canónica: alelo dominante (mayúscula) primero."""
    if a1 == a2:
        return a1 + a2
    return "".join(sorted([a1, a2], key=lambda x: (not x.isupper(), x)))

class CruceEntreCasas:
    """
    Matrimonio entre dos casas. Calcula la distribución genotípica
    esperada en la descendencia (F1) y la simula muestreando gametos
    de cada progenitor.

    Genética de población (sin selección ni deriva, una generación):

      * Frecuencia alélica en F1 = promedio de las de los padres:
            p_F1(a) = (p_padre(a) + p_madre(a)) / 2

      * Pero la distribución genotípica de F1 NO es H-W cuando las
        frecuencias parentales difieren: hay exceso de heterocigotos
        en la F1 ("heterosis" de la primera generación). La
        probabilidad exacta, muestreando un gameto de cada padre, es:

            P(AiAi) = p_padre(Ai) * p_madre(Ai)
            P(AiAj) = p_padre(Ai) * p_madre(Aj)
                    + p_padre(Aj) * p_madre(Ai)

      * Recién a partir de la F2 (cruzando F1 x F1) vuelve al
        equilibrio H-W con p = promedio parental.
    """

    def __init__(self, padre: Casa, madre: Casa, seed: int = None):
        self.padre = padre
        self.madre = madre
        self.rng = np.random.default_rng(seed)

        self.loci = []           # [(nombre, locus_padre, locus_madre), ...]
        self.frec_f1 = {}        # {nombre_locus: {alelo: frec F1}}
        self.prob_f1 = {}        # {nombre_locus: {genotipo: prob exacta F1}}

        for nombre_l, l_p, l_m in (
                ("Pelo", padre.locus_pelo, madre.locus_pelo),
                ("Ojos", padre.locus_ojos, madre.locus_ojos)):
            self.loci.append((nombre_l, l_p, l_m))
            alelos = sorted(set(l_p.alelos) | set(l_m.alelos))
            self.frec_f1[nombre_l] = {
                a: (l_p.frecuencias.get(a, 0.0) + l_m.frecuencias.get(a, 0.0)) / 2
                for a in alelos}
            self.prob_f1[nombre_l] = self._distribucion_f1(l_p, l_m)

    @staticmethod
    def _distribucion_f1(l_p: LocusHardyWeinberg,
                         l_m: LocusHardyWeinberg) -> dict:
        """Distribución genotípica exacta de F1 (producto de gametos)."""
        alelos = sorted(set(l_p.alelos) | set(l_m.alelos))
        pf = {a: l_p.frecuencias.get(a, 0.0) for a in alelos}
        pm = {a: l_m.frecuencias.get(a, 0.0) for a in alelos}
        probs = {}
        for i, a1 in enumerate(alelos):
            for j in range(i, len(alelos)):
                a2 = alelos[j]
                if i == j:
                    p = pf[a1] * pm[a1]
                else:
                    p = pf[a1] * pm[a2] + pf[a2] * pm[a1]
                probs[_etiqueta_genotipo(a1, a2)] = p
        return probs

    def simular(self, n: int = 10_000) -> dict:
        """
        Simula n hijos muestreando un gameto de cada progenitor
        para cada locus.
        """
        hijos = {}
        for nombre_l, l_p, l_m in self.loci:
            alelos = sorted(set(l_p.alelos) | set(l_m.alelos))
            p_p = [l_p.frecuencias.get(a, 0.0) for a in alelos]
            p_m = [l_m.frecuencias.get(a, 0.0) for a in alelos]

            gam_p = self.rng.choice(alelos, size=n, p=p_p)
            gam_m = self.rng.choice(alelos, size=n, p=p_m)

            conteos = {}
            for g1, g2 in zip(gam_p, gam_m):
                gen = _etiqueta_genotipo(g1, g2)
                conteos[gen] = conteos.get(gen, 0) + 1

            hijos[nombre_l] = {
                "conteos": conteos,
                "frec_obs": {g: c / n for g, c in conteos.items()},
                "frec_esp": self.prob_f1[nombre_l],
                "frec_alelica_f1": self.frec_f1[nombre_l],
            }
        return hijos

    def mostrar(self, n: int = 10_000) -> None:
        res = self.simular(n)
        print("=" * 64)
        print(f"  CRUCE {self.padre.nombre.upper()} x {self.madre.nombre.upper()}"
              f"  (N = {n} hijos)")
        print("=" * 64)
        for nombre_l, _lp, _lm in self.loci:
            r = res[nombre_l]
            fa = ", ".join(f"{a}: {f:.3f}" for a, f in r["frec_alelica_f1"].items())
            print(f"  Locus {nombre_l} | frec. alélica F1: {{{fa}}}")
            gens = [g for g in r["frec_esp"] if r["frec_esp"][g] > 0]
            for gen in sorted(gens, key=lambda g: -r["frec_esp"][g]):
                obs = r["frec_obs"].get(gen, 0.0)
                print(f"    {gen:6s}: obs {obs:7.2%} | esp {r['frec_esp'][gen]:7.2%}")
        print()

def cruzar(padre: Casa, madre: Casa, n: int = 10_000, seed: int = None) -> CruceEntreCasas:
    """Atajo: crea el cruce, lo simula y lo muestra por pantalla."""
    c = CruceEntreCasas(padre, madre, seed=seed)
    c.mostrar(n)
    return c

def _gser(a: float, x: float) -> float:
    """Gamma incompleta inferior regularizada (serie)."""
    gln = math.lgamma(a)
    ap, s, delta = a, 1.0 / a, 1.0 / a
    for _ in range(200):
        ap += 1
        delta *= x / ap
        s += delta
        if abs(delta) < abs(s) * 1e-12:
            break
    return s * np.exp(-x + a * np.log(x) - gln)


def _gcf(a: float, x: float) -> float:
    """Gamma incompleta superior regularizada (fracción continua)."""
    gln = math.lgamma(a)
    b, c = x + 1 - a, 1e300
    d = 1.0 / b
    h = d
    for i in range(1, 200):
        an = -i * (i - a)
        b += 2
        d = an * d + b
        if abs(d) < 1e-300:
            d = 1e-300
        c = b + an / c
        if abs(c) < 1e-300:
            c = 1e-300
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-12:
            break
    return np.exp(-x + a * np.log(x) - gln) * h


def _gammaincc(a: float, x: float) -> float:
    """Q(a, x) = 1 - P(a, x): cola superior de la gamma (fallback puro)."""
    if x < a + 1.0:
        return 1.0 - _gser(a, x)
    return _gcf(a, x)


def chi2_sf(x: float, df: int) -> float:
    """Función de supervivencia de la chi-cuadrado: P(X > x), X ~ chi2(df)."""
    if df <= 0:
        return 1.0
    try:
        from scipy.stats import chi2
        return float(chi2.sf(x, df))
    except ImportError:
        return float(_gammaincc(df / 2.0, x / 2.0))

def test_chi2_locus(locus: LocusHardyWeinberg, conteos: dict) -> dict:
    """
    Test chi-cuadrado de bondad de ajuste al equilibrio H-W.

        H0: los conteos genotípicos observados provienen de las
            probabilidades teóricas del locus (p², 2pq, q²).

    Grados de libertad = #genotipos - #alelos (las frecuencias
    alélicas se reestiman de los datos: k(k+1)/2 - k = k(k-1)/2).
    Si el locus es fijo (un solo genotipo), gl = 0 y no hay test.
    """
    n = sum(conteos.values())
    genotipos = list(locus.prob_genotipicas)
    obs = np.array([conteos.get(g, 0) for g in genotipos], dtype=float)
    esp = np.array([n * locus.prob_genotipicas[g] for g in genotipos])

    chi2_stat = float(((obs - esp) ** 2 / esp).sum())
    k = len(locus.alelos)
    gl = len(genotipos) - k

    if gl <= 0:
        return {"chi2": 0.0, "gl": 0, "p": 1.0,
                "veredicto": "locus fijo (sin test)"}

    p_valor = chi2_sf(chi2_stat, gl)
    veredicto = "No se rechaza H0 (ajusta a H-W)" if p_valor > 0.05 \
        else "Se rechaza H0 (NO ajusta, p < 0.05)"
    return {"chi2": chi2_stat, "gl": gl, "p": p_valor, "veredicto": veredicto}


def test_chi2_casa(casa: Casa, n: int = 10_000) -> list:
    """Aplica el test H-W a ambos loci de una casa simulada."""
    res = casa.simular(n)
    salida = []
    for locus, marg, key in ((casa.locus_pelo, res["marginal_locus_1"], "Pelo"),
                             (casa.locus_ojos, res["marginal_locus_2"], "Ojos")):
        t = test_chi2_locus(locus, marg)
        t["locus"] = f"{casa.nombre}/{key}"
        salida.append(t)
    return salida


def test_chi2_cruce(cruce: CruceEntreCasas, n: int = 10_000) -> list:
    """Compara los conteos de hijos simulados con la distribución F1 exacta."""
    res = cruce.simular(n)
    salida = []
    for nombre_l in res:
        conteos = res[nombre_l]["conteos"]
        # solo genotipos con esperanza positiva (los de prob. 0 no aportan)
        gens = [g for g in res[nombre_l]["frec_esp"] if res[nombre_l]["frec_esp"][g] > 0]
        obs = np.array([conteos.get(g, 0) for g in gens], dtype=float)
        esp = np.array([n * res[nombre_l]["frec_esp"][g] for g in gens])
        chi2_stat = float(((obs - esp) ** 2 / esp).sum())
        gl = len(gens) - 1          # las frec. parentales son parámetros fijos
        if gl <= 0:
            salida.append({"locus": f"{cruce.padre.nombre}x{cruce.madre.nombre}/{nombre_l}",
                           "chi2": 0.0, "gl": 0, "p": 1.0,
                           "veredicto": "F1 fija en este locus (sin test)"})
            continue
        p_valor = chi2_sf(chi2_stat, gl)
        veredicto = "No se rechaza H0 (ajusta a F1)" if p_valor > 0.05 \
            else "Se rechaza H0 (p < 0.05)"
        salida.append({"locus": f"{cruce.padre.nombre}x{cruce.madre.nombre}/{nombre_l}",
                       "chi2": chi2_stat, "gl": gl, "p": p_valor,
                       "veredicto": veredicto})
    return salida

def test_chi2_independencia(res_dos_locus: dict, locus_1: LocusHardyWeinberg,
                            locus_2: LocusHardyWeinberg) -> dict:
    """
    Test chi-cuadrado de independencia (tabla de contingencia) entre
    los dos loci a partir de los conteos conjuntos simulados.

        H0: los loci son independientes (no hay ligamiento).

    gl = (k1 - 1)(k2 - 1). Si algún locus es fijo, gl = 0 y no hay test.
    """
    g1 = list(locus_1.prob_genotipicas)
    g2 = list(locus_2.prob_genotipicas)
    obs = np.array([[res_dos_locus["conteos_conjunta"].get(f"{a}-{b}", 0)
                     for b in g2] for a in g1], dtype=float)
    n = obs.sum()
    esp = np.outer(obs.sum(axis=1), obs.sum(axis=0)) / n

    chi2_stat = float(((obs - esp) ** 2 / esp).sum())
    gl = (len(g1) - 1) * (len(g2) - 1)

    if gl <= 0:
        return {"chi2": 0.0, "gl": 0, "p": 1.0,
                "veredicto": "algún locus fijo (sin test)"}

    p_valor = chi2_sf(chi2_stat, gl)
    veredicto = "No se rechaza H0 (loci independientes)" if p_valor > 0.05 \
        else "Se rechaza H0 (posible ligamiento, p < 0.05)"
    return {"chi2": chi2_stat, "gl": gl, "p": p_valor, "veredicto": veredicto,
            "tabla_obs": obs, "genotipos_1": g1, "genotipos_2": g2}


# ======================================================================
#  EXTENSIONES PARA EL JUEGO  (no están en el notebook)
# ======================================================================
from collections import Counter

LOCI = ("Pelo", "Ojos")
EPS = 1e-12


def individuos_f1(cruce: CruceEntreCasas, n: int) -> list:
    """
    Simula n INDIVIDUOS de la F1 (no solo conteos por locus) con el mismo
    modelo de CruceEntreCasas.simular(): un gameto de cada progenitor por
    locus, loci independientes. Devuelve [(geno_pelo, geno_ojos), ...].

    Casa única (endogamia) == CruceEntreCasas(casa, casa): con p_padre =
    p_madre la F1 es exactamente p², 2pq, q² (Hardy-Weinberg).
    """
    columnas = {}
    for nombre_l, l_p, l_m in cruce.loci:
        alelos = sorted(set(l_p.alelos) | set(l_m.alelos))
        p_p = [l_p.frecuencias.get(a, 0.0) for a in alelos]
        p_m = [l_m.frecuencias.get(a, 0.0) for a in alelos]
        gam_p = cruce.rng.choice(alelos, size=n, p=p_p)
        gam_m = cruce.rng.choice(alelos, size=n, p=p_m)
        columnas[nombre_l] = [_etiqueta_genotipo(str(a), str(b))
                              for a, b in zip(gam_p, gam_m)]
    return list(zip(columnas["Pelo"], columnas["Ojos"]))


def firma_hipotesis(prob_f1: dict) -> tuple:
    """Huella de una distribución F1: dos hipótesis con la misma huella son
    estadísticamente indistinguibles (p. ej. Stark x Stark == Stark)."""
    return tuple((loc, tuple(sorted((g, round(p, 9))
                                    for g, p in prob_f1[loc].items() if p > EPS)))
                 for loc in LOCI)


def test_chi2_hipotesis(prob_f1: dict, individuos: list) -> dict:
    """
    Test chi-cuadrado de bondad de ajuste de una muestra de individuos
    contra la distribución F1 de una hipótesis (casa única o cruce).

        H0: los individuos provienen de esa distribución F1.

    Igual que test_chi2_cruce del notebook (gl = #genotipos posibles - 1 por
    locus), pero:
      * ambos loci se combinan (loci independientes => chi2 y gl se suman);
      * si aparece un genotipo con probabilidad 0 bajo H0, H0 es IMPOSIBLE
        (se rechaza sin calcular chi2: la contribución sería infinita).
    """
    n = len(individuos)
    loci_out, chi2_tot, gl_tot = {}, 0.0, 0
    imposible, esperados_bajos = False, False

    for idx, loc in enumerate(LOCI):
        probs = prob_f1[loc]
        conteos = Counter(ind[idx] for ind in individuos)
        posibles = [g for g, p in probs.items() if p > EPS]
        filas, chi2_l = [], 0.0
        for g in sorted(posibles, key=lambda g: -probs[g]):
            esp = n * probs[g]
            obs = conteos.get(g, 0)
            contrib = (obs - esp) ** 2 / esp
            chi2_l += contrib
            esperados_bajos |= esp < 5
            filas.append({"g": g, "obs": obs, "esp": esp,
                          "prob": probs[g], "contrib": contrib})
        imp = [g for g in conteos if probs.get(g, 0.0) <= EPS]
        for g in imp:
            filas.append({"g": g, "obs": conteos[g], "esp": 0.0,
                          "prob": 0.0, "contrib": None})
        imposible |= bool(imp)
        gl_l = max(len(posibles) - 1, 0)
        loci_out[loc] = {"filas": filas, "chi2": chi2_l, "gl": gl_l,
                         "imposibles": imp}
        chi2_tot += chi2_l
        gl_tot += gl_l

    if imposible:
        p_valor, rechaza = 0.0, True
        veredicto = "Se rechaza H0: aparece un genotipo imposible bajo esta hipótesis"
    elif gl_tot == 0:
        p_valor, rechaza = 1.0, False
        veredicto = "Sin test (loci fijos): los datos son compatibles con H0"
    else:
        p_valor = chi2_sf(chi2_tot, gl_tot)
        rechaza = p_valor < 0.05
        veredicto = ("Se rechaza H0 (p < 0.05)" if rechaza
                     else "No se rechaza H0 (p ≥ 0.05)")

    return {"n": n, "chi2": None if imposible else chi2_tot, "gl": gl_tot,
            "p": p_valor, "rechaza": rechaza, "imposible": imposible,
            "esperados_bajos": esperados_bajos, "veredicto": veredicto,
            "loci": loci_out}


def test_chi2_independencia_robusto(res_dos_locus: dict, locus_1, locus_2) -> dict:
    """
    Igual que test_chi2_independencia del notebook, pero descarta filas y
    columnas vacías de la tabla (con muestras chicas puede faltar un genotipo
    y el esperado sería 0 -> división por cero).
    """
    g1, g2 = list(locus_1.prob_genotipicas), list(locus_2.prob_genotipicas)
    obs = np.array([[res_dos_locus["conteos_conjunta"].get(f"{a}-{b}", 0)
                     for b in g2] for a in g1], dtype=float)
    filas = obs.sum(axis=1) > 0
    cols = obs.sum(axis=0) > 0
    obs = obs[filas][:, cols]
    gl = (obs.shape[0] - 1) * (obs.shape[1] - 1)
    if gl <= 0:
        return {"chi2": 0.0, "gl": 0, "p": 1.0,
                "veredicto": "algún locus fijo (sin test)"}
    esp = np.outer(obs.sum(axis=1), obs.sum(axis=0)) / obs.sum()
    chi2_stat = float(((obs - esp) ** 2 / esp).sum())
    p_valor = chi2_sf(chi2_stat, gl)
    veredicto = ("No se rechaza H0 (loci independientes)" if p_valor > 0.05
                 else "Se rechaza H0 (posible ligamiento, p < 0.05)")
    return {"chi2": chi2_stat, "gl": gl, "p": p_valor, "veredicto": veredicto}


def test_chi2_locus_conocido(locus: LocusHardyWeinberg, conteos: dict) -> dict:
    """
    Bondad de ajuste a H-W con frecuencias alélicas CONOCIDAS (las teóricas
    de la casa). Corrige test_chi2_locus del notebook: allí los esperados se
    calculan con las frecuencias teóricas pero gl = #genotipos - #alelos, que
    solo corresponde si las frecuencias se reestimaran desde la muestra (y
    entonces habría que usar los esperados con p-hat). Mezclar ambas cosas
    rechaza ~15% de las veces una población que SÍ está en H-W (nominal: 5%).

    Con parámetros fijos: gl = #genotipos - 1, igual que test_chi2_cruce.
    """
    n = sum(conteos.values())
    genotipos = list(locus.prob_genotipicas)
    obs = np.array([conteos.get(g, 0) for g in genotipos], dtype=float)
    esp = np.array([n * locus.prob_genotipicas[g] for g in genotipos])
    gl = len(genotipos) - 1
    if gl <= 0:
        return {"chi2": 0.0, "gl": 0, "p": 1.0, "veredicto": "locus fijo (sin test)"}
    chi2_stat = float(((obs - esp) ** 2 / esp).sum())
    p_valor = chi2_sf(chi2_stat, gl)
    veredicto = ("No se rechaza H0 (ajusta a H-W)" if p_valor > 0.05
                 else "Se rechaza H0 (NO ajusta, p < 0.05)")
    return {"chi2": chi2_stat, "gl": gl, "p": p_valor, "veredicto": veredicto}
