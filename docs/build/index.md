
<a id='Documentation-1'></a>

# Documentation

- [Documentation](index.md#Documentation-1)
    - [Types](index.md#Types-1)
    - [Utility Functions](index.md#Utility-Functions-1)
        - [Differentiation](index.md#Differentiation-1)
        - [Galerkin](index.md#Galerkin-1)
        - [GUI](index.md#GUI-1)
        - [Session Control](index.md#Session-Control-1)
        - [ODE](index.md#ODE-1)
    - [Examples](index.md#Examples-1)
        - [Roessler](index.md#Roessler-1)
        - [Lorenz](index.md#Lorenz-1)
    - [Index](index.md#Index-1)


<a id='Types-1'></a>

## Types

<a id='Main.Branch' href='#Main.Branch'>#</a>
**`Main.Branch`** &mdash; *Type*.



a series of `Solution`s

<a id='Main.Solution' href='#Main.Solution'>#</a>
**`Main.Solution`** &mdash; *Type*.



todo

<a id='Main.ContinuationMethod' href='#Main.ContinuationMethod'>#</a>
**`Main.ContinuationMethod`** &mdash; *Type*.



```
ContinuationMethod..
```

continuation method implementations extend  `ContinuationMethod` and the corresponding `show` and `step`. responsible for changing the project itself

<a id='Main.PC' href='#Main.PC'>#</a>
**`Main.PC`** &mdash; *Type*.



a predictor-corrector-method with step-size adaption and GUI Controls

<a id='Main.SystemCore' href='#Main.SystemCore'>#</a>
**`Main.SystemCore`** &mdash; *Type*.



todo

<a id='Main.Galerkin' href='#Main.Galerkin'>#</a>
**`Main.Galerkin`** &mdash; *Type*.



todo

<a id='Main.Session' href='#Main.Session'>#</a>
**`Main.Session`** &mdash; *Type*.



todo


<a id='Utility-Functions-1'></a>

## Utility Functions


<a id='Differentiation-1'></a>

### Differentiation

<a id='mbNewton.centralDifference' href='#mbNewton.centralDifference'>#</a>
**`mbNewton.centralDifference`** &mdash; *Function*.



```
centralDifference(homotopy, v, epsilon)
```

numerical differentiation method

<a id='mbNewton.forwardDifference' href='#mbNewton.forwardDifference'>#</a>
**`mbNewton.forwardDifference`** &mdash; *Function*.



```
forwardDifference(homotopy, v, epsilon)
```

numerical differentiation method

<a id='mbNewton.broyden' href='#mbNewton.broyden'>#</a>
**`mbNewton.broyden`** &mdash; *Function*.



```
broyden(homotopy, jacobian)
```

<a id='mbNewton.newton' href='#mbNewton.newton'>#</a>
**`mbNewton.newton`** &mdash; *Function*.



```
newton(homotopy, jacobian, v₀, ...)
```


<a id='Galerkin-1'></a>

### Galerkin

<a id='Main.findCycle' href='#Main.findCycle'>#</a>
**`Main.findCycle`** &mdash; *Function*.



```
findCycle(H, t0, y0, transientIterations, transientStepSize,
    steadyStateIterations, steadyStateStepSize)
```

<a id='Main.findCyclePoincare' href='#Main.findCyclePoincare'>#</a>
**`Main.findCyclePoincare`** &mdash; *Function*.



```
findCyclePoincare(F, y[, plane, clusterRating, nIntersections,
    maxCycles, sampleSize, transientIterations, transientStepSize,
	steadyStateStepSize])
```

extracts a single cycle of the steady state of ode `F` using poincare cuts through the `plane`.

<a id='Main.prepareCycle' href='#Main.prepareCycle'>#</a>
**`Main.prepareCycle`** &mdash; *Function*.



```
prepareCycle(data, h, P[, fac])
```

cut single cycle of length P*fac from data, resample, shift s.t. X(0)≈0, Fourier transform.

<a id='mbInterpolate.interpolateLanczos' href='#mbInterpolate.interpolateLanczos'>#</a>
**`mbInterpolate.interpolateLanczos`** &mdash; *Function*.



```
interpolateLanczos(V, a::Integer)
```

simple periodic (!) Lanczos interpolation

<a id='mbInterpolate.interpolateTrigonometric' href='#mbInterpolate.interpolateTrigonometric'>#</a>
**`mbInterpolate.interpolateTrigonometric`** &mdash; *Function*.



```
interpolateTrigonometric(a₀, a, b)
```

returns trigonometric polynomial. use with 2a,-2b and divide by 2m+1 to use with rfft coefficients.


<a id='GUI-1'></a>

### GUI

<a id='Main.ctrl' href='#Main.ctrl'>#</a>
**`Main.ctrl`** &mdash; *Function*.



```
ctrl(D, x)
```

Tuple (name::String, ::Type, init, v...)

<a id='Main.mkControlGrid' href='#Main.mkControlGrid'>#</a>
**`Main.mkControlGrid`** &mdash; *Function*.



```
mkControlGrid(D, C)
```

creates a grid of controls with labels, handlers and encapsulated storage `c` in `C` is Tuple (name::String, ::Type, init, v...)


<a id='Session-Control-1'></a>

### Session Control

<a id='Main.create' href='#Main.create'>#</a>
**`Main.create`** &mdash; *Function*.



```
create(homotopy, jacobian, projection)
```

<a id='Main.save' href='#Main.save'>#</a>
**`Main.save`** &mdash; *Function*.



```
save(filename, session[, overwrite])
```

<a id='Main.load' href='#Main.load'>#</a>
**`Main.load`** &mdash; *Function*.



```
load(filename, homotopy, jacobian, projection)
```


<a id='ODE-1'></a>

### ODE

<a id='mbRK.rk' href='#mbRK.rk'>#</a>
**`mbRK.rk`** &mdash; *Function*.



```
rk(butcherTableau)
```

returns a runge-kutta method using the respective tableau:

```
function(f, t0, y0, h, pred[, init, callback])
```

e.g. `rk1`, or `rk4`. Examines the ode `f` starting from `t0`, `y0` with fixed stepsize `h` until `pred` evalutes to `false`.


<a id='Examples-1'></a>

## Examples


<a id='Roessler-1'></a>

### Roessler


<a id='Lorenz-1'></a>

### Lorenz


<a id='Index-1'></a>

## Index

- [`Main.Branch`](index.md#Main.Branch)
- [`Main.ContinuationMethod`](index.md#Main.ContinuationMethod)
- [`Main.Galerkin`](index.md#Main.Galerkin)
- [`Main.PC`](index.md#Main.PC)
- [`Main.Session`](index.md#Main.Session)
- [`Main.Solution`](index.md#Main.Solution)
- [`Main.SystemCore`](index.md#Main.SystemCore)
- [`Main.create`](index.md#Main.create)
- [`Main.ctrl`](index.md#Main.ctrl)
- [`Main.findCycle`](index.md#Main.findCycle)
- [`Main.findCyclePoincare`](index.md#Main.findCyclePoincare)
- [`Main.load`](index.md#Main.load)
- [`Main.mkControlGrid`](index.md#Main.mkControlGrid)
- [`Main.prepareCycle`](index.md#Main.prepareCycle)
- [`Main.save`](index.md#Main.save)
- [`mbInterpolate.interpolateLanczos`](index.md#mbInterpolate.interpolateLanczos)
- [`mbInterpolate.interpolateTrigonometric`](index.md#mbInterpolate.interpolateTrigonometric)
- [`mbNewton.broyden`](index.md#mbNewton.broyden)
- [`mbNewton.centralDifference`](index.md#mbNewton.centralDifference)
- [`mbNewton.forwardDifference`](index.md#mbNewton.forwardDifference)
- [`mbNewton.newton`](index.md#mbNewton.newton)
- [`mbRK.rk`](index.md#mbRK.rk)

