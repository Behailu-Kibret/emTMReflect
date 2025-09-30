# emTMReflect — TM plane-wave reflection/transmission (MATLAB)

`emTMReflect` visualizes and animates **TM (parallel)** plane-wave incidence on a planar interface at _z = 0_.  
Medium 1 occupies _z ≤ 0_ (incident + reflected), Medium 2 occupies _z ≥ 0_ (transmitted or evanescent).  
Instantaneous fields use the convention **Re{· e^{+jωt}}**.

---

## Features
- Correct Fresnel **TM** coefficients, Snell’s law, Brewster & critical angle.
- Proper **evanescent** handling (decay into _z>0_ above θ_c).
- Optional overlays: **rays**, **angle arcs**, **angle value boxes**, **medium labels**.
- Rotated view for intuition (**+z → right**, **+x ↑**).
- **GIF/AVI** recording (Linux-friendly); **fixed-pixel** sizing for crisp output.
- Clean **name–value** API with sensible defaults.
- Optional **compute-only** mode for scripts/CI.

---

## Requirements
- MATLAB R2018b+ (R2020a+ recommended for `exportgraphics`)
- For video on Linux: use profile **`'Motion JPEG AVI'`** (default provided)

---

## Quickstart

~~~matlab
% Default demo (air → polystyrene, 1 GHz, 30°)
emTMReflect;

% Change media & angle
emTMReflect('Eps1',4,'Eps2',1,'ThetaInc',60);

% Bigger movie & save outputs (Linux-friendly)
emTMReflect('AxesSizePx',1000, ...
            'SaveGIF',true,'GIFFile','tm_demo.gif', ...
            'SaveVideo',true,'VideoProfile','Motion JPEG AVI');
~~~

---

## Function signature (all inputs optional)

~~~matlab
out = emTMReflect( ...
  'Eps1',1,'Mu1',1,'Eps2',2.56,'Mu2',1, ...
  'Freq',1e9,'ThetaInc',30,'E0',1, ...
  'GridSize',[500 500],'DomainLambda',[4 4], ...
  'RotateZRight',true,'AnimateMagnitude',true, ...
  'SlowMo',0.2,'NPeriods',12, ...
  'ShowRays',true,'ShowAngleArcs',true,'ShowAngleValues',true, ...
  'ShowMediumLabels',true,'Colormap','turbo','FreezeCLim',true, ...
  'FigSizePx',[1200 900],'AxesSizePx',900,'Renderer','opengl', ...
  'SaveGIF',false,'GIFFile','tm_wave.gif', ...
  'SaveVideo',false,'VideoFile','tm_wave.avi', ...
  'VideoProfile','Motion JPEG AVI','VideoFPS',30,'MaxFrames',300, ...
  'UpdateEvery',4,'ComputeOnly',false,'ReturnFields',false);
~~~

### Input parameters (name–value)
**Physics / source**
- `Eps1` (default `1`), `Mu1` (`1`), `Eps2` (`2.56`), `Mu2` (`1`)
- `Freq` Hz (`1e9`)
- `ThetaInc` deg from +z toward +x (`30`)
- `E0` incident amplitude (`1`)

**Spatial grid (in λ₁ units)**
- `GridSize` `[Nx Nz]` (`[500 500]`)
- `DomainLambda` `[Lx Lz]` in wavelengths (`[4 4]`)

**Visualization / overlays**
- `RotateZRight` (`true`)
- `AnimateMagnitude` (`true`)  % vs. Re{Ex}
- `SlowMo` (`0.2`), `NPeriods` (`12`)
- `ShowRays`, `ShowAngleArcs`, `ShowAngleValues`, `ShowMediumLabels` (all `true`)
- `Colormap` (`'turbo'`), `FreezeCLim` (`true`)

**Figure / axes sizing**
- `FigSizePx` `[W H]` (`[1200 900]`)
- `AxesSizePx` scalar (`900`)  ← **output pixels** for GIF/AVI
- `Renderer` (`'opengl'`)

**Recording**
- `SaveGIF` (`false`), `GIFFile` (`'tm_wave.gif'`)
- `SaveVideo` (`false`), `VideoFile` (`'tm_wave.avi'`)
- `VideoProfile` (`'Motion JPEG AVI'`)
- `VideoFPS` (`30`), `MaxFrames` (`300`)

**Performance**
- `UpdateEvery` (`4`)  % update UI & record every Nth step

**Expert**
- `ComputeOnly` (`false`)  % skip plotting/recording
- `ReturnFields` (`false`) % include Ex0/Ez0 in output

---

## Output (`struct out`)
- Fresnel: `Gamma_b`, `T_b`, `R`, `Tpow`, `Regime` = `'prop'|'crit'|'evan'`
- Angles: `ThetaBrew`, `ThetaCrit`, `ThetaInc`, `ThetaTrans` (if real)
- Wave numbers: `beta1`, `beta2`, `kx`, `kz1`, `kx2`, `kz2`
- Grid: `x`, `z`, `X`, `Z`
- Handles: `fig`, `ax` (if plotted)
- Fields (optional): `Ex0`, `Ez0` (when `'ReturnFields',true`)

---

## Examples

### 1) Critical angle exploration
~~~matlab
eps1 = 4; eps2 = 1;
theta_c = asind(sqrt(eps2/eps1));
emTMReflect('Eps1',eps1,'Eps2',eps2,'ThetaInc',theta_c, ...
            'ShowAngleValues',true,'AxesSizePx',900);
~~~

### 2) Compute-only (no graphics)
~~~matlab
out = emTMReflect('ComputeOnly',true,'Eps1',80,'Eps2',1,'ThetaInc',45);
fprintf('R=%.4f, T=%.4f, regime=%s\n', out.R, out.Tpow, out.Regime);
~~~

---

## Tips
- On Linux, prefer **`'Motion JPEG AVI'`** for `VideoProfile`.
- To increase movie size: set `AxesSizePx` (e.g., 1000 or 1200).
- If the UI feels slow, reduce `GridSize` (e.g., `[300 300]`) and/or increase `UpdateEvery`.

---

## Repository layout 
~~~text
em-tm-reflect/
├─ emTMReflect.m
├─ media/
│  ├─ screenshot.png
│  └─ demo.gif
├─ README.md
├─ LICENSE
└─ CHANGELOG.md
~~~

---

## Minimal tests

~~~matlab
% tests/test_emTMReflect.m
function tests = test_emTMReflect, tests = functiontests(localfunctions); end

function testEnergyConservationAtNormalIncidence(~)
  out = emTMReflect('ComputeOnly',true,'ThetaInc',0,'Eps1',2,'Eps2',1);
  assert(abs(out.R + out.Tpow - 1) < 1e-10);
end

function testCriticalAngle(~)
  out = emTMReflect('ComputeOnly',true,'Eps1',4,'Eps2',1, ...
                    'ThetaInc',asind(sqrt(1/4)));
  assert(strcmp(out.Regime,'crit'));
  assert(abs(out.Tpow) < 1e-12);
end
~~~

---

## License
Released under the **MIT License** (see `LICENSE`).

~~~text
MIT License
Copyright (c) 2025 Behailu Kibret

