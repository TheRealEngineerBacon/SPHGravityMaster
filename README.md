# SPHGravityMaster

This project serves as the foundation of a planet simulator. It serves to simulate three key areas.
1) Solve planetary fluid simulation and convective forcing of planet crust.
2) Simulate core heating from different sources.
3) Simulate radiative cooling of planet surface (black body radiation).
4) Simulate tectonic plate evolution without using user-function approximations.

An example calculation as of 2/22/23:
![simexample](https://user-images.githubusercontent.com/62128346/220812239-50be1e16-e3ba-4b7a-bd5c-a332efeb4760.gif)

# Future Goals
1) Simplify code into moduled format for ease of model changes.
2) Reduce core utilities into a server format for longer processing times on a dedicated machine.
3) Convert project to OpenMPI format for a distributed computation.
4) Adjust pressure model to more accurately reflect Earth-like minerals and metals.
5) Build web UI for easier, and a more entertaining, view of the latest simulation state.
