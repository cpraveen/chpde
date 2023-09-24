# Linear advection equation

Run and start visclaw

```shell
python advection_1d.py iplot=1
```

Run and generate html plots

```shell
python advection_1d.py htmlplot=1
```

Only run

```shell
python advection_1d.py
```

Now you can plot intertactively like this; start `ipython`

```python
from clawpack.pyclaw import plot
plot.interactive_plot()
```

Run in parallel

```shell
mpirun -np 2 python ./advection_1d.py use_petsc=True
```
