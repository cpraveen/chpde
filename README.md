# Computational Hyperbolic PDE

These codes are based on Clawpack and the book 

```text
R. J. LeVeque, Finite volume methods for hyperbolic problems, CUP.
```

Most of the codes have been taken from various Clawpack repositories and modified slightly in some cases.

To install Clawpack, see the instructions here

https://cpraveen.github.io/comp/clawpack.html

Directory [book](https://github.com/cpraveen/chpde/tree/master/book) contains original codes used for the book; they will not work with latest clawpack, you need to use v4.3 but I do not recommend this. Many of the codes have been converted here

https://github.com/clawpack/apps/tree/master/fvmbook

## Examples from book on "Riemann Problems and Jupyter Solutions"

Read it here 

https://www.clawpack.org/riemann_book/html/Index.html

or read the code on github

https://github.com/clawpack/riemann_book

## More examples

1. Shallow water
   1. [Tsunami from ocean onto shelf and beach](https://github.com/clawpack/geoclaw/tree/master/examples/1d_classic/ocean_shelf_beach)
   1. [Sloshing liquid in a bowl](https://github.com/clawpack/geoclaw/tree/master/examples/tsunami/bowl-slosh)

## Visualization

After generating the solution files, e.g., for classic solver using

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

---

* `Origin`: https://codeberg.org/cpraveen/chpde
* `Mirror`: https://git.sr.ht/~cpraveen/chpde
* `Mirror`: https://github.com/cpraveen/chpde
