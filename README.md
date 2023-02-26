# SPHGravityMaster

This project serves as the foundation of a planet simulator. It serves to simulate three key areas.
1) Solve planetary fluid simulation and convective forcing of planet crust.
2) Simulate core heating from different sources.
3) Simulate radiative cooling of planet surface (black body radiation).
4) Simulate tectonic plate evolution without using user-function approximations.

An example calculation as of 2/25/23:

![7ch1ov](https://user-images.githubusercontent.com/62128346/221385819-70d791ed-3460-452d-ac82-bc60d87ec1c1.gif)

An example calculation of an group of particles indicating gravitation attraction, fluid repulsion, and temperature conduction.
The brightly colored red particle is set to a constant high temp in order to add heat to the system.

# Future Goals
1) Simplify code into moduled format for ease of model changes.
2) Reduce core utilities into a server format for longer processing times on a dedicated machine.
3) Convert project to OpenMPI format for a distributed computation.
4) Adjust pressure model to more accurately reflect Earth-like minerals and metals.
5) Build web UI for easier, and a more entertaining, view of the latest simulation state.
