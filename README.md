# Non-Convective Stellar Model of the Sun
The is a RK4-shooting solver for the non-convective stellar structure ODE's, it works to produce a decent approximation for the sun's structure. It utilizes RK4 and shoots for the core values of pressure and temperature which are required for the initial boundary problem. Once it finds the correct values it uses RK4 a final time to give you the results. We cannot just use known solar values for the core pressure and temperature because convection in the real sun shifts them; therefore, they must be solved for first. 

The convection-less core temperature is ~92% accurate with the real core of the sun, and the rest of the model's values are in that range of accuracy:
P_c = 1.748e+16 Pa, T_c = 1.454e+07 K, rho_c = 8.87e+04 kg/m^3
R = 0.89 R_sun, M = 1.000 M_sun
L = 1.17 L_sun, T_surf = 4.98e+05 K.
Making this a useful approximation of the stellar structure. The solved stellar values are given by:
<img width="3000" height="1600" alt="solar-structure" src="https://github.com/user-attachments/assets/1876cfc1-889b-473d-ab3f-440c3094ac17" />
The resulting mappings for the temperature and energy released via fusion is:
<img width="1280" height="960" alt="temperturegraph" src="https://github.com/user-attachments/assets/ff00ab7e-eb5a-4212-afbd-131f5056e60e" />
<img width="1280" height="960" alt="fusiongraph" src="https://github.com/user-attachments/assets/3bea8e58-c1fc-43bb-b17d-8ea436e595ef" />

