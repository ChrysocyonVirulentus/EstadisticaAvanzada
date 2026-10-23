"""
app.py -- Servidor local del minijuego "Sangre y Linaje".

    pip install -r requirements.txt
    python app.py          ->  http://localhost:5000

Toda la estadística vive en hw_core.py (las clases de tu notebook).
Acá solo hay: metadatos visuales de las casas, la API JSON y el estado
de las partidas (en memoria; se pierde al reiniciar el servidor).
"""
import os
import threading
import uuid
import webbrowser

import numpy as np
from flask import Flask, jsonify, request, send_from_directory

import hw_core as hw

app = Flask(__name__, static_folder="static", static_url_path="/static")

# ----------------------------------------------------------------------
#  Metadatos visuales  (pos = coordenadas dentro del viewBox 600x1000 del mapa)
# ----------------------------------------------------------------------
META = {
    "Stark":     dict(emblema="🐺", color="#8a9aa6", color2="#e8edf0", lugar="Invernalia",            pos=(275, 235)),
    "Lannister": dict(emblema="🦁", color="#8b1a1a", color2="#d4af37", lugar="Roca Casterly",        pos=(172, 556)),
    "Targaryen": dict(emblema="🐉", color="#b01c1c", color2="#1a1a1a", lugar="Rocadragón",           pos=(478, 640)),
    "Martell":   dict(emblema="☀️", color="#d9622b", color2="#f2c200", lugar="Lanza del Sol",        pos=(392, 922)),
    "Blackfyre": dict(emblema="🖤", color="#3a0d12", color2="#b01c1c", lugar="Campo de Hierba Roja", pos=(322, 765)),
    "Baratheon": dict(emblema="🦌", color="#f2c200", color2="#111111", lugar="Bastión de Tormentas", pos=(432, 735)),
    "Tyrell":    dict(emblema="🌹", color="#3f7d3a", color2="#d9b84a", lugar="Jardín de Aguas",      pos=(228, 742)),
    "Hightower": dict(emblema="🗼", color="#6f7680", color2="#e8e8e8", lugar="Antigua",              pos=(150, 832)),
}

# Leyenda de alelos (del notebook). 'tipo' es solo informativo.
ALELOS = {
    "N": dict(nombre="pelo negro",    color="#1c1c1c"),
    "n": dict(nombre="pelo rubio",    color="#e6cf6e"),
    "r": dict(nombre="pelo rojo",     color="#b5371e"),
    "b": dict(nombre="pelo blanco",   color="#f1efe6"),
    "M": dict(nombre="pelo marrón",   color="#6b4423"),
    "C": dict(nombre="ojos marrones", color="#7a4a1e"),
    "v": dict(nombre="ojos verdes",   color="#3b8f4d"),
    "a": dict(nombre="ojos azules",   color="#3f7bbf"),
    "t": dict(nombre="ojos púrpura",  color="#7a3fa0"),
}

# Dificultad: individuos iniciales, y si Blackfyre ("versión difícil" del
# notebook) puede ser parte del cruce secreto.
NIVELES = {
    "escudero":  dict(titulo="Escudero",  n_inicial=5, blackfyre=False),
    "caballero": dict(titulo="Caballero", n_inicial=4, blackfyre=True),
    "señor":     dict(titulo="Señor de la guerra", n_inicial=3, blackfyre=True),
}
POOL_F1 = 200          # tamaño de la F1 simulada de cada partida
MAX_EXTRA = 10         # evidencias extra permitidas
COSTO_EVIDENCIA = 15
COSTO_TEST = 5
PUNTOS_BASE = 100
ALPHA = 0.05

RNG = np.random.default_rng()
JUEGOS = {}            # id -> estado


def num(x):
    """numpy -> tipos nativos (para JSON)."""
    return None if x is None else float(x)


def todas_las_hipotesis(casas):
    """8 casas únicas + 28 cruces = 36 hipótesis (A,B) con A<=B por orden."""
    nombres = list(casas)
    return [(a, nombres[j]) for i, a in enumerate(nombres) for j in range(i, len(nombres))]


def prob_hipotesis(casas, a, b):
    # Casa única == cruce de la casa consigo misma (p², 2pq, q²).
    return hw.CruceEntreCasas(casas[a], casas[b]).prob_f1


def etiqueta(a, b):
    return f"Casa {a}" if a == b else f"{a} × {b}"


# ----------------------------------------------------------------------
#  Atlas
# ----------------------------------------------------------------------
@app.get("/")
def index():
    return send_from_directory(app.static_folder, "index.html")


@app.get("/api/casas")
def api_casas():
    casas = hw.crear_casas(seed=None)
    out = []
    for nombre, casa in casas.items():
        m = META[nombre]
        loci = {}
        for key, locus in (("Pelo", casa.locus_pelo), ("Ojos", casa.locus_ojos)):
            loci[key] = {
                "alelos": [[a, f] for a, f in locus.frecuencias.items()],
                "genotipos": [[g, p] for g, p in locus.prob_genotipicas.items()],
            }
        out.append(dict(nombre=nombre, titulo=f"Casa {nombre}", emblema=m["emblema"],
                        color=m["color"], color2=m["color2"], lugar=m["lugar"],
                        pos=dict(x=m["pos"][0], y=m["pos"][1]), loci=loci))
    return jsonify(casas=out, alelos=ALELOS, niveles=NIVELES,
                   costos=dict(evidencia=COSTO_EVIDENCIA, test=COSTO_TEST,
                               base=PUNTOS_BASE, max_extra=MAX_EXTRA))


@app.get("/api/casa/<nombre>/simular")
def api_simular(nombre):
    casas = hw.crear_casas(seed=None)            # semilla None: muestra nueva cada vez
    if nombre not in casas:
        return jsonify(error="Casa desconocida"), 404
    n = max(10, min(int(request.args.get("n", 1000)), 100_000))
    casa = casas[nombre]
    res = casa.simular(n)
    loci = {}
    for key, locus, marg in (("Pelo", casa.locus_pelo, res["marginal_locus_1"]),
                             ("Ojos", casa.locus_ojos, res["marginal_locus_2"])):
        t = hw.test_chi2_locus_conocido(locus, marg)   # (ver nota en hw_core: el del notebook infla los falsos rechazos)
        obs_al = locus.frecuencias_alelicas_observadas(marg)
        loci[key] = {
            "genotipos": [dict(g=g, p=num(p), obs=int(marg[g]), esp=num(n * p))
                          for g, p in locus.prob_genotipicas.items()],
            "alelos_obs": [[a, num(f)] for a, f in obs_al.items()],
            "chi2": dict(chi2=num(t["chi2"]), gl=int(t["gl"]), p=num(t["p"]),
                         veredicto=t["veredicto"]),
        }
    ind = hw.test_chi2_independencia_robusto(res, casa.locus_pelo, casa.locus_ojos)
    return jsonify(n=n, loci=loci,
                   independencia=dict(chi2=num(ind["chi2"]), gl=int(ind["gl"]),
                                      p=num(ind["p"]), veredicto=ind["veredicto"]))


# ----------------------------------------------------------------------
#  Juego
# ----------------------------------------------------------------------
def ind_json(ind):
    return dict(pelo=ind[0], ojos=ind[1])


def estado_publico(j):
    return dict(id=j["id"], nivel=j["nivel"],
                individuos=[ind_json(i) for i in j["revelados"]],
                extras=j["extras"], tests=j["tests"],
                puntos=puntos_actuales(j))


def puntos_actuales(j):
    return max(10, PUNTOS_BASE - COSTO_EVIDENCIA * j["extras"] - COSTO_TEST * j["tests"])


@app.post("/api/juego/nuevo")
def juego_nuevo():
    nivel = (request.get_json(silent=True) or {}).get("nivel", "escudero")
    if nivel not in NIVELES:
        nivel = "escudero"
    cfg = NIVELES[nivel]
    casas = hw.crear_casas(seed=None)
    candidatas = [n for n in casas if cfg["blackfyre"] or n != "Blackfyre"]

    if RNG.random() < 0.5:                                   # endogamia de una casa
        a = b = str(RNG.choice(candidatas))
    else:                                                    # cruce de dos casas distintas
        a, b = (str(x) for x in RNG.choice(candidatas, size=2, replace=False))
        a, b = sorted((a, b), key=list(casas).index)

    cruce = hw.CruceEntreCasas(casas[a], casas[b], seed=int(RNG.integers(1 << 31)))
    pool = hw.individuos_f1(cruce, POOL_F1)                  # la F1 final simulada
    RNG.shuffle(pool)                                        # permutación de la F1

    j = dict(id=uuid.uuid4().hex[:12], nivel=nivel, secreto=(a, b), pool=pool,
             cursor=cfg["n_inicial"], revelados=list(pool[:cfg["n_inicial"]]),
             extras=0, tests=0, consultadas=set(), terminado=False)
    JUEGOS[j["id"]] = j
    return jsonify(estado_publico(j))


def _juego(jid):
    j = JUEGOS.get(jid)
    if j is None:
        return None, (jsonify(error="Partida inexistente (¿se reinició el servidor?)"), 404)
    if j["terminado"]:
        return None, (jsonify(error="La partida ya terminó"), 409)
    return j, None


@app.post("/api/juego/<jid>/evidencia")
def juego_evidencia(jid):
    j, err = _juego(jid)
    if err:
        return err
    if j["extras"] >= MAX_EXTRA:
        return jsonify(error="No hay más evidencia disponible"), 409
    if j["cursor"] < len(j["pool"]):                         # permutación (sin reposición)
        nuevo = j["pool"][j["cursor"]]
        j["cursor"] += 1
    else:                                                    # bootstrap (con reposición)
        nuevo = j["pool"][int(RNG.integers(len(j["pool"])))]
    j["revelados"].append(nuevo)
    j["extras"] += 1
    return jsonify(estado_publico(j))


def _leer_hipotesis(casas):
    d = request.get_json(silent=True) or {}
    a, b = d.get("a"), d.get("b")
    if a not in casas or b not in casas:
        return None
    return sorted((a, b), key=list(casas).index)


@app.post("/api/juego/<jid>/chi2")
def juego_chi2(jid):
    j, err = _juego(jid)
    if err:
        return err
    casas = hw.crear_casas(seed=None)
    h = _leer_hipotesis(casas)
    if h is None:
        return jsonify(error="Hipótesis inválida"), 400
    a, b = h
    clave = (a, b, len(j["revelados"]))
    if clave not in j["consultadas"]:                        # repetir la misma consulta no cuesta
        j["consultadas"].add(clave)
        j["tests"] += 1
    t = hw.test_chi2_hipotesis(prob_hipotesis(casas, a, b), j["revelados"])
    t.update(a=a, b=b, etiqueta=etiqueta(a, b), puntos=puntos_actuales(j),
             tests=j["tests"])
    return jsonify(t)


@app.post("/api/juego/<jid>/responder")
def juego_responder(jid):
    j, err = _juego(jid)
    if err:
        return err
    casas = hw.crear_casas(seed=None)
    h = _leer_hipotesis(casas)
    if h is None:
        return jsonify(error="Hipótesis inválida"), 400
    a, b = h
    sa, sb = j["secreto"]
    # acierta si es la misma hipótesis O una indistinguible (misma F1)
    acierto = (hw.firma_hipotesis(prob_hipotesis(casas, a, b))
               == hw.firma_hipotesis(prob_hipotesis(casas, sa, sb)))
    puntos = puntos_actuales(j) if acierto else 0
    j["terminado"] = True

    # Ranking de las 36 hipótesis contra TODA la evidencia vista
    ranking = []
    for ha, hb in todas_las_hipotesis(casas):
        t = hw.test_chi2_hipotesis(prob_hipotesis(casas, ha, hb), j["revelados"])
        ranking.append(dict(etiqueta=etiqueta(ha, hb), a=ha, b=hb, p=t["p"],
                            imposible=t["imposible"], rechaza=t["rechaza"],
                            es_verdad=hw.firma_hipotesis(prob_hipotesis(casas, ha, hb))
                            == hw.firma_hipotesis(prob_hipotesis(casas, sa, sb))))
    ranking.sort(key=lambda r: (r["imposible"], -r["p"]))

    # La F1 completa vs la esperada (para el gráfico de la revelación)
    esp = prob_hipotesis(casas, sa, sb)
    f1 = {}
    for idx, loc in enumerate(hw.LOCI):
        cont = hw.Counter(i[idx] for i in j["pool"])
        f1[loc] = [dict(g=g, esp=num(p), obs=cont.get(g, 0) / len(j["pool"]))
                   for g, p in esp[loc].items() if p > hw.EPS]
    return jsonify(acierto=acierto, puntos=puntos, secreto=dict(a=sa, b=sb, etiqueta=etiqueta(sa, sb),
                   tipo="casa" if sa == sb else "cruce"), respuesta=dict(a=a, b=b, etiqueta=etiqueta(a, b)),
                   n=len(j["revelados"]), extras=j["extras"], tests=j["tests"],
                   ranking=ranking[:6], f1=f1, pool=len(j["pool"]))


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    if os.environ.get("NO_BROWSER") != "1":
        threading.Timer(1.0, lambda: webbrowser.open(f"http://localhost:{port}")).start()
    print(f"\n  Sangre y Linaje -> http://localhost:{port}   (Ctrl+C para salir)\n")
    app.run(host="127.0.0.1", port=port, debug=False)
