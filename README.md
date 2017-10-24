# orbel
Purpose: Convertion between cartesian coordinates/velocities and orbital elements

## Usage

  1. $ xv2el ! prompt ialpha, gm, a, e, inc, capom, omega, capm
  2. $ xv2el < input_file ! screen output
  3. $ xv2el < input_file > output_file ! output to file

### Input files must have at the begin

- ialpha: conic section type

    + hyperbole (+1)
    + parabola (0)
    + ellipse (-1)

- gm: factor solar mass. ex: 39.476926421373015 [au^3 a^-2]

### Cartesian coordinates and velocities:

- x, y, v     -- cartesian coordinats [ex. au]
- vx, vy, vz  -- velocities [ex. au a^-1]

### Orbital elements:

- gm      -- factor solar mass [ex.: 39.476926421373015 au^3 a^-2]
- a       -- semi-major axis [ex.: au]
- e       -- eccentricity
- inc     -- inclination [deg]
- capom   -- longitude of the ascending node [deg]
- omega   -- argument of periapsis [deg]
- capm    -- mean anomaly [deg]

## Examples

### To convert orbital elements to cartesian coordinates

input file: el.in

```
-1                   ---> ialpha
39.476926421373015   ---> gm
0.38709900           ---> x
0.00292400           ---> y
0.00011194           ---> z
4.91338388           ---> vx
0.84639517           ---> vy
2.18152877           ---> vz
```

```Fortran
$ el2xv < el.in > xv.out
```
output file: xv.out (a, e, inc, capom, omega, capm)

```
0.38225464  0.05341021  0.00000004 -1.40041096 10.03085350  0.00001976
```

