# power_system_simulation

Power system small signal analysis, multi synchronous machine model linearisation

This code follow exactly from the formulation in the textbook Power System dynamics and Stability by Peter W.Sauer and M.A.Pai,
published in 1997 (updated in 2006)

Steps to obtain a small signal model are:
1. Get correct input data for bus, trasmission line, generators, exciters and power system stabilizers.
2. Run loadflow.m to get the operating terminal voltages and angles.
3. Get the initial condition for all variables from the large signal model.
4. Plug in all the values and construct different matrices.
5. Perform matrix operation such as Kron reduction to get matrix Asys and B.

-------------------------------------------

## Note:
This include a WSCC 9 Bus 3 Machine case and New England IEEE 39 Bus case. The code is proven to be consistence with the neumerical results in Sauer and Pai's book.
