# Replication of 'A Macroeconomic Model with a Financial Sector'

This repository contains a full computational replication of the foundational macroeconomic model presented by Markus K. Brunnermeier and Yuliy Sannikov (AER, 2014).

**Authors:**
* Mirco Tucci (mirco.tucci@carloalberto.org) - Polytechnic University of Turin, Collegio Carlo Alberto
* Davide Raciti (davide.raciti@carloalberto.org) - Polytechnic University of Turin, Collegio Carlo Alberto

**Course:**
This project was developed for the [Computational Economics Course](https://floswald.github.io/CompEcon/) in the PhD program at Collegio Carlo Alberto, taught by Prof. Florian Oswald.

---

## Overview of the Project
This package uses Julia to simulate the highly non-linear boundary value problem at the core of the Brunnermeier and Sannikov model. We utilize a precise Bisection and Shooting Method coupled with the Runge-Kutta solver Tsit5 to capture the unstable saddle-path dynamics of the economy's financial sector.

The package successfully replicates:
* **Figure 1:** Macroeconomic and financial variables with respect to experts' wealth share ($\eta$).
* **Figure 2:** The utility frontier between households and experts.

---

## How to Run the Code

To run this replication on your local machine, you need **Julia (v1.10+)**.

From the Julia REPL, position yourself in the folder Template-Replication-Package-Computational-Economics with 'cd("Template-Replication-Package-Computational-Economics")'. Then you should type 'using Pkg', then 'Pkg.activate("Our_Replication_package")', 'using Our_Replication_package', 'results= Our_Replication_package.run_all()'. If you want to visualize the plots you can write 'display(results.figure1)' and 'display(results.figure2)'


