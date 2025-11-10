This file contains names and definitions to explain the naming convention used within the code.



**Process Components:**

a = CO2, b = H2O, c = CO, d = H2, e = O2



**Co-Electrolysis\_FO**

* **user\_CO2:** The amount of input CO2 to the electrolyser as estimated within the "fuel demand\_FO" module
* **user\_temp:** Desired system temperature
* **"out","prod", "recy":** denote components in the product, outlet and recycle streams respectively, while **"feed"** and "**syn"** represent the feed and syngas streams respectively
* **"elec\_sp\_E", "heat\_sp\_E" :** Specific electricity and heat demand respectively



**Comidy\_FT\_FO**

* **"C2to8", "C9to16" and "c17plus":** denote Fischer-Tropsch product pseudocomponents i.e. light ends, middle distillates and heavy cuts



**Hydrocracking\_FO**

* **mass\_gases, mass\_mid\_distillates, mass\_heavy\_cuts:** Resulting masses of light ends, middle distillates and heavy cuts from the hydrocracker respectively



**Fuel\_demand\_FO**

* **slope\_prodvsCO2:** ratio of produced fuel to input CO2
* **slope\_ATF\_arr, y\_int\_arr:** representing the slope and y-intercept respectively



**Main\_FO\_20240316**

* Array naming convention describes the total masses of components in and out of a process unit.

E.g.

**"tl\_CO2\_input\_CoElec"** stores the CO2 requirement for the electrolyser as calculated within the fuel demand module

**"tl\_H2O\_input\_CoElec"** stores the inlet H2O to the electrolyser while

**"tl\_H2O\_required\_CoElec"** stores the stoichiometric H2O requirement for coelectrolysis



* **"sum\_useful":** represents the ATF produced
* **"sum\_heavy":** Total amount of heavy cuts produced in process
* "**sum\_heavy\_sell":** Amount of heavy cuts available to be sold at the end of the process
* **"sum\_light\_sell":** Amount of light ends available to be sold at the end of the process



* **H2\_stock:** A stock is created to hold the recycling that remains each month. It is assumed that components can readily be separated for reuse





