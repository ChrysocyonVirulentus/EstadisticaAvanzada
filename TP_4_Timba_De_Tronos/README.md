# Timba de Tronos — Hardy-Weinberg en Poniente

Prototipo local: backend Python (tu notebook) + frontend web sin dependencias.

## Correr
```bash
pip install -r requirements.txt     # flask, numpy 
python app.py                       # abre http://localhost:5000
```
Variables opcionales: `PORT=8080`, `NO_BROWSER=1`.
Las fuentes medievales vienen de Google Fonts; sin internet cae a serif del sistema.

## Estructura
| Archivo | Qué es |
|---|---|
| `hw_core.py` | clases del notebook (verbatim) + sección final "EXTENSIONES PARA EL JUEGO" |
| `app.py` | API Flask, metadatos de casas (emblema, color, posición en el mapa), lógica de partida |
| `static/index.html` | Mapa SVG de Poniente |
| `static/app.js` | Atlas + juego + gráficos (SVG a mano) |
| `static/style.css` | Estilo |

## Cómo se juega
1. **Atlas**: clic en una casa → alelos de Pelo/Ojos, ecuación H-W, barras observado vs esperado y curva p²/2pq/q².
2. **Juego**: el servidor elige un linaje secreto (una casa o un cruce), simula su F1 (200 individuos) y te muestra 3–5.
   Calculás frecuencias alélicas, probás hipótesis con χ² y proclamás veredicto.
   Más evidencia = otro individuo de la F1 simulada (permutación sin reposición; bootstrap si se agotara).
3. Puntos: 100 − 15 por evidencia extra − 5 por consulta χ² nueva (mínimo 10 si acertás; 0 si fallás).

## Para tocar
- **Frecuencias / casas**: `crear_casas()` en `hw_core.py`. Si agregás una casa, agregala también en `META` (`app.py`) con su posición (x,y) en el viewBox 600×1000 del mapa.
- **Costos y niveles**: constantes arriba en `app.py` (`NIVELES`, `COSTO_*`, `POOL_F1`).
- **Mapa**: es un SVG estilizado dibujado a mano (`index.html`). Para usar una imagen propia, reemplazá el contenido del `<svg>` por un `<image>` y reajustá `pos` en `META`.
