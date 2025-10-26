# Isentropic vortex with AMR

This problem is similar to `../isentropic_vortex.py`.

```shell
make .output
```

The solution files should now be in `_output` directory. Start `ipython` and run the following commands

```shell
import clawpack.visclaw.Iplotclaw as Iplotclaw
ip = Iplotclaw()
ip.plotloop()
```

Press enter to get the first plot and for later frames.

Or, make html plots

```shell
make .plots
open _plots/_PlotIndex.html
```
