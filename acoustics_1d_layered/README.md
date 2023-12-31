
# Acoustic waves propagating in a layered medium 

As shown in Figures 9.8 and 9.9.

The transmission-based limiter of Fogarty and LeVeque is used when `mylim=1` in `setprob.data`.  This works better than the standard limiter, and is implemented in `trlimit.f` as called from `step1.f90` (a modified version of this library routine is in this directory).

See differences in `step1`

```
diff ./step1.f90 $(CLAW)/classic/src/1d/step1.f90
```

## Standard wave limiter

Set `mylim=0` in `setrun.py` and run

```shell
make .plots
open _plots/_PlotIndex.html
```

## Transmissive wave limiter

Set `mylim=1` in `setrun.py` and run

```shell
make .plots
open _plots/_PlotIndex.html
```

## Exercise

Run with 480 and 2400 cells, and make figures as in the book. Try both `dt = dx` and `dt = 0.8*dx`.
