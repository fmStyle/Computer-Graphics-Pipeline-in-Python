Short demo renderer written in Python (pygame + numpy). This repository contains a simple software pipeline that:

- Generates a Bézier surface via de Casteljau.
- Transforms and projects triangles with a basic projection matrix.
- Rasterizes triangles into a coarse grid with depth testing and a simple Blinn–Phong-like shading.

This project is intended as a portfolio piece. It is slow and experimental — not suitable for production.

**Requirements**

- Python 3.8+
- numpy
- pygame

Install with:

```
pip install -r requirements.txt
```

**Run**

```
python pipeline_complete.py
```

**Controls**

- Left / Right arrows: rotate the surface.
- Up / Down arrows: adjust vertical offset (height).

**Quick notes & parameters**

- The grid resolution is created in `pipeline_complete.py` with `grid = Grid(100, 100, 7, (5, 18), (19, 11))`. Lower the first two numbers to improve speed (for example `Grid(60, 60, ...)`).
- The Bézier sampling density is controlled by the `divisions` variable inside `Grid.casteljau()` (default is `5`). Lower that to reduce triangle count.
- Toggle `wireframe` in the script to skip expensive fragment shading and instead draw edges.

**Performance & limitations**

- This is a software rasterizer that paints into a pixel grid of coarse square fragments. It is intentionally unoptimized and therefore slow for larger grids and finer Bézier tessellations.
- Main bottlenecks: Python-level loops over fragments, per-triangle per-pixel math, and frequent object/array updates. Vectorizing more work in `numpy` or moving hot loops to C/C++ (or using numba) will give large speedups.

