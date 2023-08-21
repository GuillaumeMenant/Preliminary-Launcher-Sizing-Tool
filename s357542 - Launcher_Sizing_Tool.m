close all
clear all
clc

%% Constant parameters

g0 = 9.81; % m/s^2, standard gravity acceleration of orbiting object
mu = 3.986*10^14; % m3/s-2, gravitational parameter of orbiting object
R_E = 6378*10^3; % m, Earth radius

%% Inputs

%%%%%%% Launch site inputs %%%%%%%
launch_latitude = 5.16; %degrees

%%%%%%% Orbit inputs %%%%%%%

h_perigee = 200*10^3; % Altitude of the perigee
h_apogee =  200*10^3; % Altitude of the apogee
h = h_perigee; % altitude of insertion of payload, usually at perigee
a_orbit = (h_perigee+R_E+h_apogee+R_E)/2; % semi-major axis
inclination = 5.16; %degrees, final orbit inclination
FPA = 0; %degrees, flight path angle at orbit entry
azimuth = asind(cosd(inclination)/cosd(launch_latitude)); %launch azimuth angle

%%%%%% General launcher inputs %%%%%%%
m_payload = 19139; % desired payload mass, in kg.
n_stages = 2; % number of stages of the launcher, maximum available = 3


%%%%%%% First stage inputs %%%%%%%

D_1 = 4.5; % diameter of the first stage, in m.
Isp1 = 295; % Engine effective ISP, in seconds.
sigma1 = 0.05; % Structural factor, = m_structure/(m_structure+m_propellant).
Thrust_1 = 850*10^3; %Thrust of one main engine, in N.
Pressure_1 = 10*10^(6); %Pressure in the combustion chamber, in Pa.
T_W_1 = 1.5; %Thrust to weight ratio of the stage. 
n_engines_1 = 1; % number of engines, put "1" as input, the program will calculate the necessary amount automatically

t_start_1 = 1; % time to start the engine, in seconds. 

f_residual_1 = 0.02; % fraction of residual propellant remaining

f_ox_ratio_1 = 1/3.5; %fuel/oxidizer mass mixture ratio

AR_ox_1 = sqrt(2); % aspect ratio of the oxidizer tank, depending on the desired shape
AR_fuel_1 = sqrt(2); % aspect ratio of the fuel tank, depending on the desired shape

density_ox_1 = 1140; %density of oxidizer propellant, kg/m3
density_fuel_1 = 423; %density of fuel propellant, kg/m3

ullage_ox_1 = 0.05; %factor accounting for the added volume needed for the pressured gases 
ullage_fuel_1 = 0.05; %factor accounting for the added volume needed for the pressured gases

        %%% Only cryogenic cases

cryo_ox_1 = 1; %If the oxidizer is cryogenic, put "1" as an input, otherwise put "0"
cryo_fuel_1 = 1; %If the fuel is cryogenic, put "1" as an input, otherwise put "0"

alpha_ox_1 = 23.4*10^(-6); % Thermal expansion coefficient of the oxidizer tank material
alpha_fuel_1 = 23.4*10^(-6); % Thermal expansion coefficient of the fuel tank material
delta_temp_ox_1 = 293 - 54.4; % Difference between the boiling temperature of the cryogenic oxidizer
delta_temp_fuel_1 = 293 - 20.3; % Difference between the boiling temperature of the cryogenic fuel

        %%% Only hybrid cases
hybrid_1 = 0; % If the stage has a hybrid engine, put "1" as an input, otherwise put "0"

        %%% Only reusable cases
reusability_1 = 1; % if stage is reusable, select "1", else, select "0";

density_cold_gas_1 = 250; % density of the chosen cold gas, in kg/m3
ISP_cold_gas_1 = 80; % ISP of the cold gas engine, in seconds
delta_V_landing_per_thruster_1 = 10; % Available deltaV for a cold gas thruster, in m/s
T_cold_gas_1 = 500; % Thrust of cold gas thruster, in N 
n_thruster_1 = 8; % Number of cold gas thrusters

            %%%%%%%%% If VTVL (Vertical Take-Off Vertical Landing)

VTVL_1 = 1; % if the stage is a VTVL, select "1" otherwise select "0"
VTVL_DRL_1 = 1; % If the stage is Down-Range landing, select "1", otherwise select "0" if return to launch site landing

n_landing_legs_1 = 4; % number of landing legs on the stage
density_landing_leg_1 = 1750; % density of the material of the landing legs, in kg/m3
D_landing_leg_1 = 1; % Diameter of one landing leg, in m
thickness_landing_leg_1 = 0.05; % Thickness of one landing leg, defined as a hollow cylinder, in m

n_fins_1 = 4; % Number of aerodynamic fins on the stage
density_aerodynamic_fin_1 = 150; % density of the material of the aerodynamic fins, in kg/m3
length_aerodynamic_fin_1 = 1; % Length of the side of one aerodynamic fin, in m
thickness_aerodynamic_fin_1 = 0.3; % Thickness of one aerodynamic fin, in m

density_TPS_VTVL_1 = 200; % density of the material used for TPS, in kg/m3
thickness_TPS_VTVL_1 = 0.02; % thickness of the TPS layer, in m 

            %%%%%%%%% If VTHL (Vertical Take-Off Horizontal Landing)

VTHL_1 = 0;  % if the stage is a VTHL, select "1" otherwise select "0"
VTHL_1_FB = 0; % If the stage is Flying Back autonomously to runway, select "1" otherwise select "0" if air-capture

density_wings_1 = 1800; % density of the material used for the wings, in kg/m3
thickness_wings_1 = 1.5; % thickness of one wing's root, in m
thickness_body_flap_1 = 0.5; % thickness of the body flap, in m 
length_body_flap_1 = 1; % length of the body flap, in m

density_TPS_VTHL_1 = 96; % density of the material used for TPS, in kg/m3
thickness_TPS_VTHL_1_wings = 0.01; % thickness of the TPS layer on the wings, in m 
thickness_TPS_VTHL_1_body = 0.01; % thickness of the TPS layer on the body, in m 

m_engine_FB_1 = 988; % mass of the chosen airbreathing engine for the FB configuration, in kg
Thrust_per_engine_FB_1 = 54000; % thrust of the chosen airbreathing engine for the FB configuration, in N
t_landing_FB_1 = 600; % duration of the autonomous fly-back manoeuvre for the FB configuration, in s
specific_fuel_cons_FB_1 = 8.1*10^(-6); % the specific fuel consumption of the airbreathing engine for the FB configuration, in kg/N.s


%%%%%%% Second stage inputs %%%%%%%

D_2 = 4.5; % diameter of the second stage, in m.
Isp2 = 360; % Engine effective ISP, in seconds.
sigma2 = 0.064; % Structural factor, = m_structure/(m_structure+m_propellant).
Thrust_2 = 980*10^3; % Thrust of one main engine, in N.
Pressure_2 = 10*10^(6); % Pressure in the combustion chamber, in Pa.
T_W_2 = 0.7; % Thrust to weight ratio of the stage.
n_engines_2 = 1; % number of engines, put "1" as input, the program will calculate the necessary amount automatically

t_start_2 = 1; % time to start the engine, in seconds.

f_residual_2 = 0.02; % fraction of residual propellant remaining

f_ox_ratio_2 = 1/3.5; % fuel/oxidizer mass mixture ratio

AR_ox_2 = sqrt(2); % aspect ratio of the tank, depending on the desired shape
AR_fuel_2 = sqrt(2); % aspect ratio of the tank, depending on the desired shape

density_ox_2 = 1140; % density of oxidizer propellant, kg/m3
density_fuel_2 = 423; % density of fuel propellant, kg/m3

ullage_ox_2 = 0.05; % factor accounting for the added volume needed for the pressured gases
ullage_fuel_2 = 0.05; % factor accounting for the added volume needed for the pressured gases

        %%% Only cryogenic cases
cryo_ox_2 = 1; %If the oxidizer is cryogenic, put "1" as an input, otherwise put "0"
cryo_fuel_2 = 1; %If the fuel is cryogenic, put "1" as an input, otherwise put "0"

alpha_ox_2 = 23.4*10^(-6); % Thermal expansion coefficient of the oxidizer tank material
alpha_fuel_2 = 23.4*10^(-6); % Thermal expansion coefficient of the fuel tank material
delta_temp_ox_2 = 293 - 54.4; % Difference between the boiling temperature of the cryogenic oxidizer
delta_temp_fuel_2 = 293 - 20.3; % Difference between the boiling temperature of the cryogenic fuel

        %%% Only hybrid cases
hybrid_2 = 0; % If the stage has a hybrid engine, put "1" as an input, otherwise put "0"

        %%% Only reusable cases

reusability_2 = 1; % if stage is reusable, select "1", else, select "0";

density_cold_gas_2 = 250; % density of the chosen cold gas, in kg/m3
ISP_cold_gas_2 = 80; % ISP of the cold gas engine, in seconds
delta_V_landing_per_thruster_2 = 10; % Available deltaV for a cold gas thruster, in m/s
T_cold_gas_2 = 500; % Thrust of cold gas thruster, in N
n_thruster_2 = 8; % Number of cold gas thrusters

            %%%%%%%%% If VTVL (Vertical Take-Off Vertical Landing)

VTVL_2 = 0; % if the stage is a VTVL, select "1" otherwise select "0"
VTVL_DRL_2 = 0; % If the stage is Down-Range landing, select "1", otherwise select "0" if return to launch site landing

n_landing_legs_2 = 4; % number of landing legs on the stage
density_landing_leg_2 = 1750; % density of the material of the landing legs, in kg/m3
D_landing_leg_2 = 2; % Diameter of one landing leg, in m
thickness_landing_leg_2 = 0.1; % Thickness of one landing leg, defined as a hollow cylinder, in m

n_fins_2 = 4; % Number of aerodynamic fins on the stage
density_aerodynamic_fin_2 = 150; % density of the material of the aerodynamic fins, in kg/m3
length_aerodynamic_fin_2 = 1; % Length of the side of one aerodynamic fin, in m
thickness_aerodynamic_fin_2 = 0.3; % Thickness of one aerodynamic fin, in m

density_TPS_VTVL_2 = 200; % density of the material used for TPS, in kg/m3
thickness_TPS_VTVL_2 = 0.02; % thickness of the TPS layer, in m

            %%%%%%%%% If VTHL (Vertical Take-Off Horizontal Landing)

VTHL_2 = 0;  % if VTHL, select "1" otherwise select "0"
VTHL_2_FB = 0; % If the stage is Flying Back autonomously to runway, select "1" otherwise select "0" if air-capture

density_wings_2 = 1800; % density of the material used for the wings, in kg/m3
thickness_wings_2 = 1.5; % thickness of one wing's root, in m
thickness_body_flap_2 = 0.5; % thickness of the body flap, in m
length_body_flap_2 = 1; % length of the body flap, in m

density_TPS_VTHL_2 = 96; % density of the material used for TPS, in kg/m3
thickness_TPS_VTHL_2_wings = 0.01; % thickness of the TPS layer on the wings, in m 
thickness_TPS_VTHL_2_body = 0.01; % thickness of the TPS layer on the body, in m 

m_engine_FB_2 = 988; % mass of the chosen airbreathing engine for the FB configuration, in kg
Thrust_per_engine_FB_2 = 54000; % thrust of the chosen airbreathing engine for the FB configuration, in N
t_landing_FB_2 = 3600; % duration of the autonomous fly-back manoeuvre for the FB configuration, in s
specific_fuel_cons_FB_2 = 8.1*10^(-3); % the specific fuel consumption of the airbreathing engine for the FB configuration, in kg/N.s

%%%%%%% Thirdstage inputs %%%%%%%

D_3 = 2.2; % diameter of the first stage, in m.
Isp3 = 314; % Engine effective ISP, in seconds.
sigma3 = 0.18; % Structural factor, = m_structure/(m_structure+m_propellant).
Thrust_3 = 25*10^3; % Thrust of one main engine, in N.
Pressure_3 = 30*10^(6); % Pressure in the combustion chamber, in Pa.
T_W_3 = 0.6; % Thrust to weight ratio of the stage.
n_engines_3 = 1; % number of engines, put "1" as input, the program will calculate the necessary amount automatically

t_start_3 = 0.5; % time to start the engine, in seconds.

f_residual_3 = 0.02; % fraction of residual propellant remaining

f_ox_ratio_3 = 3.1; % fuel/oxidizer mass mixture ratio

AR_ox_3 = sqrt(2); % aspect ratio of the oxidizer tank, depending on the desired shape
AR_fuel_3 = sqrt(2); % aspect ratio of the fuel tank, depending on the desired shape

density_ox_3 = 1140; % density of oxidizer propellant, kg/m3
density_fuel_3 = 900; % density of fuel propellant, kg/m3

ullage_ox_3 = 0.05; % factor accounting for the added volume needed for the pressured gases
ullage_fuel_3 = 0.05; % factor accounting for the added volume needed for the pressured gases

        %%% Only cryogenic cases
cryo_ox_3 = 1; % If the oxidizer is cryogenic, put "1" as an input, otherwise put "0"
cryo_fuel_3 = 0; % If the fuel is cryogenic, put "1" as an input, otherwise put "0"

alpha_ox_3 = 23.4*10^(-6); % Thermal expansion coefficient of the oxidizer tank material
alpha_fuel_3 = 23.4*10^(-6); % Thermal expansion coefficient of the fuel tank material
delta_temp_ox_3 = 293 - 20.7; % Difference between the boiling temperature of the cryogenic oxidizer
delta_temp_fuel_3 = 293 - 20.7; % Difference between the boiling temperature of the cryogenic fuel

        %%% Only hybrid cases
hybrid_3 = 0; % If the stage has a hybrid engine, put "1" as an input, otherwise put "0"

        %%% Only reusable cases

reusability_3 = 0; % if stage is reusable, select "1", else, select "0";

density_cold_gas_3 = 250; % density of the chosen cold gas, in kg/m3
ISP_cold_gas_3 = 80; % ISP of the cold gas engine, in seconds
delta_V_landing_per_thruster_3 = 10; % Available deltaV for a cold gas thruster, in m/s
T_cold_gas_3 = 500; % Thrust of cold gas thruster, in N
n_thruster_3 = 8; % Number of cold gas thrusters

           %%%%%%%%% If VTVL (Vertical Take-Off Vertical Landing)

VTVL_3 = 0; % if the stage is a VTVL, select "1" otherwise select "0"
VTVL_DRL_3 = 0; % If the stage is Down-Range landing, select "1", otherwise select "0" if return to launch site landing

n_landing_legs_3 = 4; % number of landing legs on the stage
density_landing_leg_3 = 1750; % density of the material of the landing legs, in kg/m3
D_landing_leg_3 = 1; % Diameter of one landing leg, in m
thickness_landing_leg_3 = 0.1; % Thickness of one landing leg, defined as a hollow cylinder, in m

n_fins_3 = 4; % Number of aerodynamic fins on the stage
density_aerodynamic_fin_3 = 150;% density of the material of the aerodynamic fins, in kg/m3
length_aerodynamic_fin_3 = 1; % Length of the side of one aerodynamic fin, in m
thickness_aerodynamic_fin_3 = 0.3; % Thickness of one aerodynamic fin, in m

density_TPS_VTVL_3 = 200; % density of the material used for TPS, in kg/m3
thickness_TPS_VTVL_3 = 0.02; % thickness of the TPS layer, in m

            %%%%%%%%% If VTHL (Vertical Take-Off Horizontal Landing)

VTHL_3 = 0;  % if the stage is a VTHL, select "1" otherwise select "0"
VTHL_3_FB = 0; % If the stage is Flying Back autonomously to runway, select "1" otherwise select "0" if air-capture

density_wings_3 = 1800; % density of the material used for the wings, in kg/m3
thickness_wings_3 = 1.5; % thickness of one wing's root, in m
thickness_body_flap_3 = 0.5; % thickness of the body flap, in m
length_body_flap_3 = 1; % length of the body flap, in m

density_TPS_VTHL_3 = 96; % density of the material used for TPS, in kg/m3
thickness_TPS_VTHL_3_wings = 0.01; % thickness of the TPS layer on the wings, in m 
thickness_TPS_VTHL_3_body = 0.01; % thickness of the TPS layer on the body, in m 


m_engine_FB_3 = 988; % mass of the chosen airbreathing engine for the FB configuration, in kg
Thrust_per_engine_FB_3 = 54000; % thrust of the chosen airbreathing engine for the FB configuration, in N
t_landing_FB_3 = 3600; % duration of the autonomous fly-back manoeuvre for the FB configuration, in s
specific_fuel_cons_FB_3 = 8.1*10^(-3); % the specific fuel consumption of the airbreathing engine for the FB configuration, in kg/N.s
  
%%%%%%% Booster stage inputs %%%%%%%
n_boosters = 0; % number of boosters for the whole launcher
Isp_booster = 250; % ISP of one booster
sigma_booster = 0.145; % Structural factor of one booster, = m_structure/(m_structure+m_propellant)
T_W_booster = 1.9; %Thrust to weight ratio of the boosters
density_prop_0 = 1000; %density of the solid propellant

%% Delta-V calculations

%%%%% Orbital Velocities %%%%%%
V_orbit_burnout = sqrt(mu*(2/(h_perigee+R_E) - 1/a_orbit)); % Orbit burnout velocity, in m/s
V_orbit_burnout_south = -V_orbit_burnout*cosd(FPA)*cosd(azimuth); % Orbit burnout velocity, south component, in m/s
V_orbit_burnout_east = V_orbit_burnout*cosd(FPA)*sind(azimuth); % Orbit burnout velocity, east component, in m/s
V_orbit_burnout_zenith = V_orbit_burnout*sind(FPA); % Orbit burnout velocity, zenith component, in m/s

V_launch_site = 465.1*cosd(launch_latitude); %Launch site speed due to earth rotation, in m/s

V_needed_south = V_orbit_burnout_south; % Needed velocity with added launch site velocity, south component, in m/s
V_needed_east = V_orbit_burnout_east - V_launch_site; % Needed velocity with added launch site velocity, east component, in m/s
V_needed_zenith = V_orbit_burnout_zenith; % Needed velocity with added launch site velocity, zenith component, in m/s

V_needed = sqrt(V_needed_south^2 + V_needed_east^2 + V_needed_zenith^2); % Needed velocity with added launch site velocity, in m/s

%%%%% Losses Delta-Vs %%%%%%

V_grav_losses = 0.8*sqrt(2*g0*h/(1+h/R_E)); % Gravity losses deltaV, multiplied by 0.8 for more accurate estimation, in m/s
V_drag_losses =  0.05*V_orbit_burnout; % Drag losses, in m/s
V_steering_losses = 0.1*V_orbit_burnout; %Steering losses, in m/s

V_losses = V_grav_losses + V_drag_losses + V_steering_losses; % Total deltaV losses, in m/s

%%%%% Final DeltaV %%%%%%
delta_V_needed = V_needed + V_losses; % The total and final deltaV needed for the launcher to perform the mission, in m/s

%%%%% DeltaV distribution between stages %%%%%%
boosters_delta_V = 0; % percentage of total DeltaV ensured by the boosters
stage_1_deltaV = 0.45; % percentage of total DeltaV ensured by the first stage
stage_2_deltaV = 0.55; % percentage of total DeltaV ensured by the second stage
stage_3_deltaV = 0; % percentage of total DeltaV ensured by the third stage

delta_V_1_per_booster = boosters_delta_V*delta_V_needed/n_boosters; % DeltaV ensured per booster
delta_V_1_core = stage_1_deltaV*delta_V_needed; % DeltaV ensured by the first stage

delta_V_2 = stage_2_deltaV*delta_V_needed; % DeltaV ensured by the second stage

delta_V_3 = stage_3_deltaV*delta_V_needed; % DeltaV ensured by the third stage

% Definition of the DeltaV values depending on the number of stages 

if n_stages == 2

    delta_V_3 = 0;

elseif n_stages == 1
    
    delta_V_2 = 0;
    delta_V_3 = 0;

end

%% Initial sizing calculations - masses of propellant and overall structure
index = 3; % variable initiated to ensure the running of the while loop and the update of the overall payload mass to add the reusability to the system
count = 1; % count of the iterations needed to converge once the reusability systems have been added

while index == 3

    %%%% Stage 3 masses
    
    if n_stages == 3

        m_payload_step_3 = m_payload; % definition of the step 3 payload mass

        m_ratio_3 = exp(delta_V_3/(g0*Isp3)); % definition of the mass ratio step 3
        m_prop_3 = m_payload*((m_ratio_3 - 1)*(1 - sigma3)/(1 - m_ratio_3*sigma3)); % definition of the propellant mass step 3
        m_struct_3 = m_prop_3*sigma3/(1-sigma3); % definition of the structure mass step 3
        
        m_prop_3_startup = T_W_3*(m_prop_3+m_struct_3+m_payload_step_3)*t_start_3/Isp3; %definition of the startup propellant needed step 3
        m_prop_3 = m_prop_3 + m_prop_3_startup; %adding up the startup propellant needed to the overall propellant mass of step 3
        
        m_prop_residual_3 = m_prop_3*f_residual_3; % definition of the residual propellant step 3
        m_prop_3 = m_prop_3 + m_prop_residual_3; % adding up the residual propellant needed to the overall propellant mass of step 3
       
        m_payload_step_2 = m_struct_3 + m_prop_3 + m_payload; % definition of the 'payload' for step 2
    
    else
    
        m_payload_step_2 = m_payload; % if there is no third stage, then the payload mass of step 2 is the overall payload mass
    
        m_prop_3 = 0; % since there is no stage 3, no propellant mass defined
        m_struct_3 = 0; % since there is no stage 3, no structure mass defined
    
    end
    
    %%%%% Stage 2 masses
    
    if n_stages >= 2
    
        m_ratio_2 = exp(delta_V_2/(g0*Isp2)); % definition of the mass ratio step 2
        m_prop_2 = m_payload_step_2*((m_ratio_2 - 1)*(1 - sigma2)/(1 - m_ratio_2*sigma2)); % definition of the propellant mass step 2
        m_struct_2 = m_prop_2*sigma2/(1-sigma2); % definition of the structure mass step 2
        
        m_prop_2_startup = T_W_2*(m_prop_2+m_struct_2+m_payload_step_2)*t_start_2/Isp2; %definition of the startup propellant needed step 2
        m_prop_2 = m_prop_2 + m_prop_2_startup; %adding up the startup propellant needed to the overall propellant mass of step 2
        
        m_prop_residual_2 = m_prop_2*f_residual_2; % definition of the residual propellant step 2
        m_prop_2 = m_prop_2 + m_prop_residual_2; % adding up the residual propellant needed to the overall propellant mass of step 2
     
        m_payload_step_1 = m_struct_2 + m_prop_2 + m_payload_step_2; % definition of the 'payload' for step 1
            
    else
                
        m_payload_step_1 = m_payload; % if there is no second stage, then the payload mass of step 1 is the overall payload mass
        
        m_prop_2 = 0; % since there is no stage 2, no propellant mass defined
        m_struct_2 = 0;% since there is no stage 2, no structure mass defined
    
    end
    
    %%%%% Stage 1 masses
    
    m_ratio_1 = exp(delta_V_1_core/(g0*Isp1)); % definition of the mass ratio step 1
    m_prop_1 = m_payload_step_1*((m_ratio_1 - 1)*(1 - sigma1)/(1 - m_ratio_1*sigma1)); % definition of the propellant mass step 1
    m_struct_1 = m_prop_1*sigma1/(1-sigma1);% definition of the structure mass step 1

    m_prop_1_startup = T_W_1*(m_prop_1+m_struct_1+m_payload_step_1)*t_start_1/Isp1; %definition of the startup propellant needed step 1
    m_prop_1 = m_prop_1 + m_prop_1_startup; %adding up the startup propellant needed to the overall propellant mass of step 1
    
    m_prop_residual_1 = m_prop_1*f_residual_1; % definition of the residual propellant step 1
    m_prop_1 = m_prop_1 + m_prop_residual_1; % adding up the residual propellant needed to the overall propellant mass of step 1
    
   
    %%%%% Stage 0 masses, booster stage
    
    if n_boosters > 0
    
        m_payload_step_0 = m_payload_step_1 + m_struct_1 + m_prop_1; % definition of the 'payload' for step 0
            
        m_ratio_0 = exp(delta_V_1_per_booster/(g0*Isp_booster)); % definition of the mass ratio step 0
        m_prop_0 = m_payload_step_0*((m_ratio_0 - 1)*(1 - sigma_booster)/(1 - m_ratio_0*sigma_booster)); % Definition of the propellant mass step 0
        m_struct_0 = m_prop_0*sigma_booster/(1-sigma_booster); % Definition of the structure mass step 0
    
    end

          
    %%%%%% Final mass with added reusability if needed

    if reusability_3 == 1 && count >= 2 % If the third stage is reusable and the count is above 2 to already have an expendable configuration to work with
            
        if VTVL_3 == 1 % If the third stage is a VTVL vehicle
        
            %%%% Propellant penalty for landing thrust
        
            if VTVL_DRL_3 == 1 % Different propellant penalty applied depending ont eh type of VTVL vehicle
                m_prop_penalty_reus_3 = 0.1*m_prop_3; % If DRL, then only 10% penalty applied, in kg
            else
                m_prop_penalty_reus_3 = 0.2*m_prop_3; % If RTLS, then 20% penalty applied, in kg
            end
        
            %%%% Attitude control thruster
                     
            m_prop_cold_gas_3 = n_thruster_3*(m_struct_3 + m_prop_penalty_reus_3)*(exp(delta_V_landing_per_thruster_3/(g0*ISP_cold_gas_3)) - 1); % propellant mass of the cold gas thrusters, in kg
            vol_tank_cold_gas_3 = m_prop_cold_gas_3/density_cold_gas_3; % volume of the cold gas tanks, in m3
            m_tank_cold_gas_3 = 12.16*vol_tank_cold_gas_3; % mass of the cold gas tanks, in kg
            m_engine_cold_gas_3 = T_cold_gas_3*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20))+59; % mass of the engine of the cold gas thrusters, in kg
        
            %%%%% Landing surfaces
            
            length_per_landing_leg_3 = 18/47 * height_stage3; % length of one landing leg, in m
            vol_per_landing_leg_3 = pi*(D_landing_leg_3/2)^2*length_per_landing_leg_3 - pi*((D_landing_leg_3-2*thickness_landing_leg_3)/2)^2*length_per_landing_leg_3;% volume of one landing leg, in m3
            m_per_landing_leg_3 = vol_per_landing_leg_3*density_landing_leg_3; % mass of one landing leg, in kg
            m_total_landing_legs_3 = n_landing_legs_3 * m_per_landing_leg_3; % mass of all the landing legs, in kg
        
            %%%% Aerodynamic surfaces
        
            vol_aerodynamic_fin_3 = length_aerodynamic_fin_3^2*thickness_aerodynamic_fin_3; % volume of one aerodynamic fin, in m3
            m_per_aerodynamic_fin_3 = vol_aerodynamic_fin_3*density_aerodynamic_fin_3; % mass of one aerodynamic fin, in kg
            m_total_aerodynamic_fin_3 = m_per_aerodynamic_fin_3*n_fins_3; % mass of all the aerodynamic fins, in kg
    
            %%%% TPS 
    
            vol_TPS_3 = thickness_TPS_VTVL_3*PI*(D_3/2)^2; % volume of the TPS layer, in m3
            m_TPS_3 = vol_TPS_3*density_TPS_VTVL_3; % mass of the TPS layer, in kg
        
        
            %%%% Added reusability masses to the overall system
        
            m_prop_3_added_reusability = m_prop_penalty_reus_3 + m_prop_cold_gas_3; % the added reusability mass for propellant purposes, in kg
            m_struct_3_added_reusability = m_TPS_3 + m_total_aerodynamic_fin_3 + m_total_landing_legs_3 + m_tank_cold_gas_3; % the added reusability mass for structure purposes, in kg

            %%%% Total mass reusability
            m_added_reusability_3 =  m_total_aerodynamic_fin_3 + m_TPS_3 +  m_total_landing_legs_3 + m_tank_cold_gas_3 + m_engine_cold_gas_3 + m_prop_penalty_reus_3 + m_prop_cold_gas_3; % the total reusability mass to be added to the system, in kg


            %%%% Propellant and structure mass with reusability

            m_prop_3 = m_prop_3 + m_prop_3_added_reusability; % the final propellant mass with the added reusability propellant, in kg
            m_struct_3 = m_struct_3 + m_struct_3_added_reusability;  % the final structure mass with the added reusability structure, in kg   
            

         
        elseif VTHL_3 == 1

            %%%% Air breathing engines

            if VTHL_3_FB == 1 % if the vehicle is a FB configuration
                

                Thrust_FB_3 = g0*1/3*(m_engine_FB_3); % Thrust of the airbreathing engine, in N
                n_engines_FB_3 = Thrust_FB_3/Thrust_per_engine_FB_3; % number of needed airbreathing engines to ensure the desire thrust is produced
                m_struct_engines_FB_3 = m_engine_FB_3*n_engines_FB_3; % the total dry mass of the airbreathing engines, in kg
                m_prop_FB_3 = specific_fuel_cons_FB_3*Thrust_FB_3*t_landing_FB_3; % the mass of propellant needed for the airbreathing engines, in kg

                vol_tank_FB_3 = m_prop_FB_3/71; % the volume of the propellant tank for airbreathing engine, in m3; with hydrogen density used
                m_tank_FB_3 = 12.16*vol_tank_FB_3; % the mass of the propellant tank for airbreathing engines, in kg

                h_dome_FB_3 = (D_3/2)/AR_ox_3; % the height of the dome for the propellant tank of airbreathing engine, in m
                h_cyl_FB_3 = 4*vol_tank_FB_3/(pi*D_3^2) - 2*D_3/(3*AR_ox_3); % the height of the cylinder for the propellant tank of airbreathing engine, in m
                
                A_tank_FB_3 = pi*D_3*h_cyl_FB_3 + 2*1.6234*pi*(D_3/2)^2; % the area of the tank of air breathing engine propellant, in m2
                
                if h_cyl_FB_3 < 0 % if the height of the cylinder is negative, that means the height of two domes is enough to contain all the propellant needed
                    h_cyl_FB_3 = 0;
                end

                                
            else
                m_prop_FB_3 = 0;
                m_tank_FB_3 = 0;
                m_struct_engines_FB_3 = 0;
                h_dome_FB_3 = 0;
                h_cyl_FB_3 = 0;
            end             
            
            %%%% Attitude control thruster
                     
            m_prop_cold_gas_3 = n_thruster_3*(m_struct_3 + m_prop_FB_3)*(exp(delta_V_landing_per_thruster_3/(g0*ISP_cold_gas_3)) - 1);% propellant mass of the cold gas thrusters, in kg
            vol_tank_cold_gas_3 = m_prop_cold_gas_3/density_cold_gas_3; % volume of the tanks for the cold gas thrusters, in m3
            m_tank_cold_gas_3 = 12.16*vol_tank_cold_gas_3; % mass of the tank for the cold gas thrusters propellant, in kg
            m_engine_cold_gas_3 = T_cold_gas_3*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20))+59; % mass of the engine for the cold gas thrusters, in kg

            %%%% Tail wing

            length_tail_wing_3 = (6.8/32)*height_stage3; % length of one tail wing, normalized based on the space shuttle dimensions, in m
            width_tail_wing_3 = length_tail_wing_3*tand(45);% width of one tail wing, normalized based on the space shuttle dimensions, in m
            m_tail_wing_3 = 1/2*(110*(m_struct_3_stored(1)+m_prop_FB_3)*0.454*(1.5*4.4)*(2*width_tail_wing_3+D_3)*3.28*(width_tail_wing_3*length_tail_wing_3*3.28^2 + length_tail_wing_3*D_3*3.28^2)^(0.77)/(thickness_wings_3*3.28))*10^(-6); %mass of the tail wings, in kg
            
            %%%% Wings

            length_wing_3 = (10/32)*height_stage3; % length of one wing, normalized based on the space shuttle dimensions, in m
            width_wing_3 = length_wing_3*tand(45);% width of one wing, normalized based on the space shuttle dimensions, in m
            m_wings_3 = (110*(m_struct_3_stored(1)+m_prop_FB_3)*0.454*(1.5*4.4)*(2*width_wing_3+D_3)*3.28*(width_wing_3*length_wing_3*3.28^2 + length_wing_3*D_3*3.28^2)^(0.77)/(thickness_wings_3*3.28))*10^(-6); %mass of both wings, in kg

            %%%% Body flap
    
            vol_body_flap_3 = length_body_flap_3*D_3*thickness_body_flap_3;% volume of one body flap, in m3
            m_body_flap_3 = vol_body_flap_3*density_wings_3; % mass of one body flap, in kg

            %%%% TPS
            
            length_tail_wing_TPS_3 = length_tail_wing_3/10; % length of the TPS along the tail wings, in m
            length_wing_TPS_3 = length_wing_3/10; % length of the TPS along the wings, in m
            vol_TPS_3 = length_tail_wing_TPS_3*thickness_TPS_VTHL_3_wings*width_tail_wing_3 + 2*length_wing_TPS_3*thickness_TPS_VTHL_3_wings*width_wing_3 + thickness_TPS_VTHL_3_body*PI*(D_3/2)^2; % volume of the TPS layers applied including the body TPS layer, in m3
            m_TPS_3 = vol_TPS_3*density_TPS_VTHL_3; % mass of the TPS layers, in kg

            %%%% Landing gear

            m_landing_gear_3 = 0.05*(m_struct_3_stored(1)+m_prop_3_stored(1)); % mass of the landing gear, in kg

            %%%% Added reusability masses to the overall system
        
            m_prop_3_added_reusability = m_prop_cold_gas_3 + m_prop_FB_3; % mass of the propellant added reusability mass, in kg
            m_struct_3_added_reusability = m_struct_engines_FB_3 + m_landing_gear_3 + m_TPS_3 + m_body_flap_3 + m_wings_3 + m_tail_wing_3 + m_engine_cold_gas_3 + m_tank_cold_gas_3 + m_tank_FB_3; % mass of the added structural reusability mass, in kg

            %%%% Total mass reusability
           
            m_added_reusability_3 =  m_struct_engines_FB_3 + m_landing_gear_3 + m_tail_wing_3 + m_wings_3 + m_body_flap_3 + m_TPS_3 + m_tank_cold_gas_3 + m_engine_cold_gas_3 + m_prop_cold_gas_3 + m_tank_FB_3 + m_prop_FB_3; % the total reusability mass to add to the overall system, in kg

            %%%% Propellant and structure mass with reusability

            m_prop_3 = m_prop_3 + m_prop_3_added_reusability; % the final propellant mass with the added reusability propellant, in kg
            m_struct_3 = m_struct_3 + m_struct_3_added_reusability;   % the final structural mass with the added reusability structure, in kg  
        end

    else

        m_added_reusability_3 = 0;

        h_cyl_FB_3 = 0;
        h_dome_FB_3 = 0; 
    end


    
     if reusability_2 == 1 && count >= 2 % If the second stage is reusable and the count is above 2 to already have an expendable configuration to work with
    
        if VTVL_2 == 1 % If the second stage is a VTVL vehicle
        
            %%%% Propellant penalty for landing thrust
        
            if VTVL_DRL_2 == 1 % Different propellant penalty applied depending on the type of VTVL vehicle
                m_prop_penalty_reus_2 = 0.1*m_prop_2; % If DRL, then only 10% penalty applied, in kg
            else
                m_prop_penalty_reus_2 = 0.2*m_prop_2; % If RTLS, then 20% penalty applied, in kg
            end
        
            %%%% Attitude control thruster
                     
            m_prop_cold_gas_2 = n_thruster_2*(m_struct_2 + m_prop_penalty_reus_2)*(exp(delta_V_landing_per_thruster_2/(g0*ISP_cold_gas_2)) - 1); % propellant mass of the cold gas thrusters, in kg
            vol_tank_cold_gas_2 = m_prop_cold_gas_2/density_cold_gas_2; % volume of the cold gas tanks, in m3
            m_tank_cold_gas_2 = 12.16*vol_tank_cold_gas_2; % mass of the cold gas tanks, in kg
            m_engine_cold_gas_2 = T_cold_gas_2*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20))+59; % mass of the engine of the cold gas thrusters, in kg
        
            %%%% Landing surfaces
            
            length_per_landing_leg_2 = 18/47 * height_stage2; % length of one landing leg, in m
            vol_per_landing_leg_2 = pi*(D_landing_leg_2/2)^2*length_per_landing_leg_2 - pi*((D_landing_leg_2-2*thickness_landing_leg_2)/2)^2*length_per_landing_leg_2; % volume of one landing leg, in m3
            m_per_landing_leg_2 = vol_per_landing_leg_2*density_landing_leg_2; % mass of one landing leg, in kg
            m_total_landing_legs_2 = n_landing_legs_2 * m_per_landing_leg_2; % mass of all the landing legs, in kg
        
            %%%% Aerodynamic surfaces
        
            vol_aerodynamic_fin_2 = length_aerodynamic_fin_2^2*thickness_aerodynamic_fin_2; % volume of one aerodynamic fin, in m3
            m_per_aerodynamic_fin_2 = vol_aerodynamic_fin_2*density_aerodynamic_fin_2; % mass of one aerodynamic fin, in kg
            m_total_aerodynamic_fin_2 = m_per_aerodynamic_fin_2*n_fins_2; % mass of all the aerodynamic fins, in kg
    
            %%%% TPS 
    
            vol_TPS_2 = thickness_TPS_VTVL_2*pi*(D_2/2)^2; % volume of the TPS layer, in m3
            m_TPS_2 = vol_TPS_2*density_TPS_VTVL_2; % mass of the TPS layer, in kg
        
        
            %%%% Added reusability masses to the overall system
        
            m_prop_2_added_reusability = m_prop_penalty_reus_2 + m_prop_cold_gas_2; % the added reusability mass for propellant purposes, in kg
            m_struct_2_added_reusability = m_TPS_2 + m_total_aerodynamic_fin_2 + m_total_landing_legs_2 + m_tank_cold_gas_2; % the added reusability mass for structure purposes, in kg
    
            %%%% Total mass reusability
            m_added_reusability_2 =  m_total_aerodynamic_fin_2 + m_TPS_2 +  m_total_landing_legs_2 + m_tank_cold_gas_2 + m_engine_cold_gas_2 + m_prop_penalty_reus_2 + m_prop_cold_gas_2; % the total reusability mass to be added to the system, in kg

            %%%% Propellant and structure mass with reusability

            m_prop_2 = m_prop_2 + m_prop_2_added_reusability; % the final propellant mass with the added reusability propellant, in kg
            m_struct_2 = m_struct_2 + m_struct_2_added_reusability; % the final structure mass with the added reusability structure, in kg      
    
         
        elseif VTHL_2 == 1
    
            %%%% Air breathing engines
    
                       
            if VTHL_2_FB == 1 % if the vehicle is a FB configuration
                
                Thrust_FB_2 = g0*1/3*(m_engine_FB_2); % Thrust of the airbreathing engine, in N
                n_engines_FB_2 = Thrust_FB_2/Thrust_per_engine_FB_2; % number of needed airbreathing engines to ensure the desired thrust is produced
                m_struct_engines_FB_2 = n_engines_FB_2*m_engine_FB_2; % the total dry mass of the airbreathing engines, in kg
                m_prop_FB_2 = specific_fuel_cons_FB_2*Thrust_FB_2*t_landing_FB_2; % the mass of propellant needed for the airbreathing engines, in kg
    
                vol_tank_FB_2 = m_prop_FB_2/71; % the volume of the propellant tank for airbreathing engine, in m3; with hydrogen density used
                m_tank_FB_2 = 12.16*vol_tank_FB_2; % the mass of the propellant tank for airbreathing engines, in kg
    
                h_dome_FB_2 = (D_2/2)/AR_ox_2; % the height of the dome for the propellant tank of airbreathing engine, in m
                h_cyl_FB_2 = 4*vol_tank_FB_2/(pi*D_2^2) - 2*D_2/(3*AR_ox_2); % the height of the cylinder for the propellant tank of airbreathing engine, in m
                
                A_tank_FB_2 = pi*D_2*h_cyl_FB_2 + 2*1.6234*pi*(D_2/2)^2; % the area of the tank of air breathing engine propellant, in m2
                
                if h_cyl_FB_2 < 0 % if the height of the cylinder is negative, that means the height of two domes is enough to contain all the propellant needed
                    h_cyl_FB_2 = 0;
                end
    
            else
                m_prop_FB_2 = 0;
                m_tank_FB_2 = 0;
                m_struct_engines_FB_2 = 0;
                h_cyl_FB_2 = 0;
                h_dome_FB_2 = 0; 
    
            end 
            
            %%%% Attitude control thruster
                     
            m_prop_cold_gas_2 = n_thruster_2*(m_struct_2 + m_prop_FB_2)*(exp(delta_V_landing_per_thruster_2/(g0*ISP_cold_gas_2)) - 1); % propellant mass of the cold gas thrusters, in kg
            vol_tank_cold_gas_2 = m_prop_cold_gas_2/density_cold_gas_2; % volume of the cold gas tanks, in m3
            m_tank_cold_gas_2 = 12.16*vol_tank_cold_gas_2; % mass of the cold gas tanks, in kg
            m_engine_cold_gas_2 = T_cold_gas_2*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20))+59; % mass of the engine of the cold gas thrusters, in kg
    
            %%%% Tail wing
    
            length_tail_wing_2 = (6.8/32)*height_stage2; % length of one tail wing, normalized based on the Space Shuttle data, in m 
            width_tail_wing_2 = length_tail_wing_2*tand(40); % width of one tail wing, normalized based on the Space Shuttle data, in m 
            m_tail_wing_2 = 1/2*(110*(m_struct_2_stored(1)+m_prop_FB_2)*0.454*(1.5*4.4)*(2*width_tail_wing_2+D_2)*3.28*(width_tail_wing_2*length_tail_wing_2*3.28^2 + length_tail_wing_2*D_2*3.28^2)^(0.77)/(thickness_wings_2*3.28))*10^(-6); %mass of the tail wings, in kg
            
    
            %%%% Wings
    
            length_wing_2 = (10/32)*height_stage2; % length of one wing, normalized based on the space shuttle dimensions, in m
            width_wing_2 = length_wing_2*tand(40); % width of one wing, normalized based on the space shuttle dimensions, in m
            m_wings_2 = (110*(m_struct_2_stored(1)+m_prop_FB_2)*0.454*(1.5*4.4)*(2*width_wing_2+D_2)*3.28*(width_wing_2*length_wing_2*3.28^2 + length_wing_2*D_2*3.28^2)^(0.77)/(thickness_wings_2*3.28))*10^(-6); %mass of both wings
                
            %%%% Body flap
    
            vol_body_flap_2 = length_body_flap_2*D_2*thickness_body_flap_2; % volume of a body flap, in m3
            m_body_flap_2 = vol_body_flap_2*density_wings_2; % mass of a body flap, in kg
    
            %%%% TPS
            
            length_tail_wing_TPS_2 = length_tail_wing_2/10; % length of the TPS layer applied on the tail wings, in m
            length_wing_TPS_2 = length_wing_2/10; % length of the TPS layer applied on the wings, in m
            vol_TPS_2 = length_tail_wing_TPS_2*thickness_TPS_VTHL_2_wings*width_tail_wing_2 + 2*length_wing_TPS_2*thickness_TPS_VTHL_2_wings*width_wing_2 + thickness_TPS_VTHL_2_body*PI*(D_2/2)^2; % volume of the TPS layers applied, including the layer on the body, in m3
            m_TPS_2 = vol_TPS_2*density_TPS_VTHL_2; % mass of the TPS, in kg

            %%%% Landing gear

            m_landing_gear_2 = 0.05*(m_struct_2_stored(1)+m_prop_2_stored(1)); % mass of the landing gear, in kg
    
            %%%% Added reusability masses to the overall system
        
            m_prop_2_added_reusability = m_prop_cold_gas_2 + m_prop_FB_2; % mass of the propellant added reusability mass
            m_struct_2_added_reusability = m_struct_engines_FB_2 + m_landing_gear_2 + m_TPS_2 + m_body_flap_2 + m_wings_2 + m_tail_wing_2 + m_engine_cold_gas_2 + m_tank_cold_gas_2 + m_tank_FB_2; % mass of the added structural reusability mass
    
            %%%% Total mass reusability
           
            m_added_reusability_2 = m_struct_engines_FB_2 + m_landing_gear_2 + m_tail_wing_2 + m_wings_2 + m_body_flap_2 + m_TPS_2 + m_tank_cold_gas_2 + m_engine_cold_gas_2 + m_prop_cold_gas_2 + m_tank_FB_2 + m_prop_FB_2; % total added reusability mass to the expendable configuration, in kg
    
            %%%% Propellant and structure mass with reusability

            m_prop_2 = m_prop_2 + m_prop_2_added_reusability; % the final propellant mass with the added reusability propellant
            m_struct_2 = m_struct_2 + m_struct_2_added_reusability; % the final structural mass with the added reusability structure         
    
        end
    
        else
            m_added_reusability_2 = 0;
            h_dome_FB_2 = 0;
            h_cyl_FB_2 = 0;  

                  
     end
        
     if reusability_1 == 1 && count >= 2 % If the first stage is reusable
                
            if VTVL_1 == 1 % If the first stage is a VTVL vehicle
            
                %%%% Propellant penalty for landing thrust
    
                if VTVL_DRL_1 == 1 % Different propellant penalty applied depending on the type of VTVL vehicle
                    m_prop_penalty_reus_1 = 0.1*m_prop_1_stored(1); % If DRL, then only 10% penalty applied, in kg
                else
                    m_prop_penalty_reus_1 = 0.2*m_prop_1_stored(1); % If RTLS, then 20% penalty applied, in kg
                end
                
            
                %%%%% Attitude control thruster
                         
                m_prop_cold_gas_1 = n_thruster_1*(m_struct_1 + m_prop_penalty_reus_1)*(exp(delta_V_landing_per_thruster_1/(g0*ISP_cold_gas_1)) - 1); % propellant mass of the cold gas thrusters, in kg
                vol_tank_cold_gas_1 = m_prop_cold_gas_1/density_cold_gas_1; % volume of the cold gas tanks, in m3
                m_tank_cold_gas_1 = 12.16*vol_tank_cold_gas_1; % mass of the cold gas tanks, in kg
                m_engine_cold_gas_1 = T_cold_gas_1*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20))+59; % mass of the engine of the cold gas thrusters, in kg
            
                %%%%% Landing surfaces
                
                length_per_landing_leg_1 = 18/47 * height_stage1; % length of one landing leg, in m
                vol_per_landing_leg_1 = pi*(D_landing_leg_1/2)^2*length_per_landing_leg_1 - pi*((D_landing_leg_1-2*thickness_landing_leg_1)/2)^2*length_per_landing_leg_1; % volume of one landing leg, in m3
                m_per_landing_leg_1 = vol_per_landing_leg_1*density_landing_leg_1; % mass of one landing leg, in kg
                m_total_landing_legs_1 = n_landing_legs_1 * m_per_landing_leg_1; % mass of all the landing legs, in kg
            
                %%%% Aerodynamic surfaces
            
                vol_aerodynamic_fin_1 = length_aerodynamic_fin_1^2*thickness_aerodynamic_fin_1; % volume of one aerodynamic fin, in m3
                m_per_aerodynamic_fin_1 = vol_aerodynamic_fin_1*density_aerodynamic_fin_1;% mass of one aerodynamic fin, in kg
                m_total_aerodynamic_fin_1 = m_per_aerodynamic_fin_1*n_fins_1; % mass of all the aerodynamic fins, in kg
        
                %%%%% TPS 
        
                vol_TPS_1 = thickness_TPS_VTVL_1*pi*(D_1/2)^2; % volume of the TPS layer, in m3
                m_TPS_1 = vol_TPS_1*density_TPS_VTVL_1; % mass of the TPS layer, in kg
                    
                %%%%% Added reusability masses to the overall system
    
                m_prop_1_added_reusability = m_prop_penalty_reus_1 + m_prop_cold_gas_1; % the added reusability mass for propellant purposes, in kg
                m_struct_1_added_reusability = m_TPS_1 + m_total_aerodynamic_fin_1 + m_total_landing_legs_1 + m_tank_cold_gas_1; % the added reusability mass for structure purposes, in kg
    
                %%%%% Total mass reusability
                m_added_reusability_1 =  m_total_aerodynamic_fin_1 + m_TPS_1 +  m_total_landing_legs_1 + m_tank_cold_gas_1 + m_engine_cold_gas_1 + m_prop_penalty_reus_1 + m_prop_cold_gas_1; % the total reusability mass to be added to the system, in kg

                %%%%% Propellant and structure mass with reusability

                m_prop_1 = m_prop_1 + m_prop_1_added_reusability; % the final propellant mass with the added reusability propellant, in kg
                m_struct_1 = m_struct_1 + m_struct_1_added_reusability; % the final structure mass with the added reusability structure, in kg
                
    
             
            elseif VTHL_1 == 1
    
                %%%%% Air breathing engines
    
                if VTHL_1_FB == 1 % if the vehicle is a FB configuration
                    
                    Thrust_FB_1 = g0*1/3*(m_struct_1_stored(1)+m_prop_1_stored(1)); % Thrust of the airbreathing engine, in N
                    n_engines_FB_1 = Thrust_FB_1/Thrust_per_engine_FB_1; % number of needed airbreathing engines to ensure the desired thrust is produced
                    m_struct_engines_FB_1 = n_engines_FB_1*m_engine_FB_1; % the total dry mass of the airbreathing engines, in kg
                    m_prop_FB_1 = specific_fuel_cons_FB_1*Thrust_FB_1*t_landing_FB_1; % the mass of propellant needed for the airbreathing engines, in kg
    
                    vol_tank_FB_1 = m_prop_FB_1/70; % the volume of the propellant tank for airbreathing engine, in m3; with hydrogen density used
                    m_tank_FB_1 = 12.16*vol_tank_FB_1; % the mass of the propellant tank for airbreathing engines, in kg
    
                    h_dome_FB_1 = (D_1/2)/AR_ox_1; % the height of the dome for the propellant tank of airbreathing engine, in m
                    h_cyl_FB_1 = 4*vol_tank_FB_1/(pi*D_1^2) - 2*D_1/(3*AR_ox_1); % the height of the cylinder for the propellant tank of airbreathing engine, in m
                    
                    A_tank_FB_1 = pi*D_1*h_cyl_FB_1 + 2*1.6234*pi*(D_1/2)^2; % the area of the tank of air breathing engine propellant, in m2

                    if h_cyl_FB_1 < 0 % if the height of the cylinder is negative, that means the height of two domes is enough to contain all the propellant needed
                        h_cyl_FB_1 = 0;
                    end
    
                else


                    m_prop_FB_1 = 0;
                    m_tank_FB_1 = 0;
                    n_engines_FB_1 = 0;
                    m_struct_engines_FB_1 = 0;
    
                end
      
                
                %%%%% Attitude control thruster
                         
                m_prop_cold_gas_1 = n_thruster_1*(m_struct_1 + m_prop_FB_1)*(exp(delta_V_landing_per_thruster_1/(g0*ISP_cold_gas_1)) - 1);% propellant mass of the cold gas thrusters, in kg
                vol_tank_cold_gas_1 = m_prop_cold_gas_1/density_cold_gas_1;% volume of the cold gas tanks, in m3
                m_tank_cold_gas_1 = 12.16*vol_tank_cold_gas_1;% mass of the cold gas tanks, in kg
                m_engine_cold_gas_1 = T_cold_gas_1*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20))+59;% mass of the engine of the cold gas thrusters, in kg
    
                %%%%% Tail wing
    
                length_tail_wing_1 = (6.8/32)*height_stage1;% length of one tail wing, normalized based on the Space Shuttle data, in m 
                width_tail_wing_1 = length_tail_wing_1*tand(40);% width of one tail wing, normalized based on the Space Shuttle data, in m 
                m_tail_wing_1 = 2*(110*((m_struct_1_stored(1)+m_prop_FB_1)*2.20*(1.5*4.4)*(2*width_tail_wing_1+D_1)*3.28*(width_tail_wing_1*length_tail_wing_1*3.28^2 + length_tail_wing_1*D_1*3.28^2))^(0.77)/(thickness_wings_1*3.28))*10^(-6)*0.454; %mass of all tail wings, in kg
    
                %%%%% Wings
    
                length_wing_1 = (10/32)*height_stage1;% length of one wing, normalized based on the Space Shuttle data, in m 
                width_wing_1 = length_wing_1*tand(40);% width of one tail wing, normalized based on the Space Shuttle data, in m 
                m_wings_1 = 2*(110*((m_struct_1_stored(1)+m_prop_FB_1)*2.20*(1.5*4.4)*(2*width_wing_1+D_1)*3.28*(width_wing_1*length_wing_1*3.28^2 + length_wing_1*D_1*3.28^2))^(0.77)/(thickness_wings_1*3.28))*10^(-6)*0.454; %mass of all wings, in kg
                                   
                %%%%% Body flap
        
                vol_body_flap_1 = length_body_flap_1*D_1*thickness_body_flap_1;% volume of a body flap, in m3
                m_body_flap_1 = vol_body_flap_1*density_wings_1;% mass of a body flap, in kg
    
                %%%%% TPS
                
                length_tail_wing_TPS_1 = length_tail_wing_1/10;% length of the TPS layer applied on the tail wings, in m
                length_wing_TPS_1 = length_wing_1/10;% length of the TPS layer applied on the wings, in m
                vol_TPS_1 = length_tail_wing_TPS_1*thickness_TPS_VTHL_1_wings*width_tail_wing_1 + 2*length_wing_TPS_1*thickness_TPS_VTHL_1_wings*width_wing_1 + thickness_TPS_VTHL_1_body*pi*(D_1/2)^2; % volume of the TPS layers applied, including the layer on the body, in m3
                m_TPS_1 = vol_TPS_1*density_TPS_VTHL_1;% mass of the TPS layers, in kg

                %%%%% Landing gear

                m_landing_gear_1 = 0.05*(m_struct_1_stored(1)+m_prop_1_stored(1));% mass of the landing gear, in kg
    
                %%%%% Added reusability masses to the overall system
            
                m_prop_1_added_reusability = m_prop_cold_gas_1 + m_prop_FB_1;% mass of the propellant added reusability mass
                m_struct_1_added_reusability = m_struct_engines_FB_1 + m_landing_gear_1 + m_TPS_1 + m_body_flap_1 + m_wings_1 + m_tail_wing_1 + m_engine_cold_gas_1 + m_tank_cold_gas_1 + m_tank_FB_1;% mass of the added structural reusability mass
    
                %%%%%% Total mass reusability
               
                m_added_reusability_1 =  m_struct_engines_FB_1 + m_landing_gear_1 + m_tail_wing_1 + m_wings_1 + m_body_flap_1 + m_TPS_1 + m_tank_cold_gas_1 + m_engine_cold_gas_1 + m_prop_cold_gas_1 + m_tank_FB_1 + m_prop_FB_1;% total added reusability mass to the expendable configuration, in kg


                %%%%%% Propellant and structure mass with reusability

                m_prop_1 = m_prop_1 + m_prop_1_added_reusability;% the final propellant mass with the added reusability propellant
                m_struct_1 = m_struct_1 + m_struct_1_added_reusability;% the final structural mass with the added reusability structure

            end

     else 
           
            m_added_reusability_1 = 0;
            h_dome_FB_1 = 0;
            h_cyl_FB_1 = 0;
             
        
     end

    m_struct_1_stored(count) = m_struct_1; %storing the structural masses through the different iterations to have access to them 
    m_prop_1_stored(count) = m_prop_1; %storing the propellant masses through the different iterations to have access to them 

    m_struct_2_stored(count) = m_struct_2;%storing the structural masses through the different iterations to have access to them 
    m_prop_2_stored(count) = m_prop_2;%storing the propellant masses through the different iterations to have access to them 

    m_struct_3_stored(count) = m_struct_3;%storing the structural masses through the different iterations to have access to them 
    m_prop_3_stored(count) = m_prop_3;%storing the propellant masses through the different iterations to have access to them 

    m_fuel_1 = m_prop_1*(f_ox_ratio_1/(1+f_ox_ratio_1));% calculation of the fuel mass of the first stage based on the previous calculation, in kg
    m_ox_1 = m_prop_1*(1/(1+f_ox_ratio_1));% calculation of the oxidizer mass of the first stage based on the previous calculation, in kg

    if n_stages >= 2
        m_fuel_2 = m_prop_2*(f_ox_ratio_2/(1+f_ox_ratio_2)); % calculation of the fuel mass of the second stage based on the previous calculation, in kg
        m_ox_2 = m_prop_2*(1/(1+f_ox_ratio_2));% calculation of the oxidizer mass of the second stage based on the previous calculation, in kg
    else
        m_fuel_2 = 0;
        m_ox_2 = 0;
    end

    if n_stages == 3

        m_fuel_3 = m_prop_3*(f_ox_ratio_3/(1+f_ox_ratio_3));% calculation of the fuel mass of the third stage based on the previous calculation, in kg
        m_ox_3 = m_prop_3*(1/(1+f_ox_ratio_3));% calculation of the oxidizer mass of the third stage based on the previous calculation, in kg
    else

        m_fuel_3 = 0;
        m_ox_3 = 0;
    end



    %%%%% Lift-off mass calculation 

    if n_boosters > 0
        m_lift_off(count) = n_boosters*m_prop_0 + n_boosters*m_struct_0 + m_struct_1 + m_prop_1 + m_struct_2 + m_prop_2 + m_struct_3 + m_prop_3 + m_payload; %lift off mass if there is boosters, in kg
        
    elseif n_boosters == 0
        m_lift_off(count) = m_struct_1 + m_prop_1 + m_struct_2 + m_prop_2 + m_struct_3 + m_prop_3 + m_payload; %lift off mass if there isn't boosters, in kg
        
    end    
    
    %% Initial sizing calculations - volumes of propellant tanks
    
    %%%%% Stage 3 volume 
    
        % Ideal volume calculation
    if hybrid_3 == 1 % if the vehicle is hybrid, different volume calculation 
        vol_tank_fuel_3 = 1.12*m_fuel_3/(density_fuel_3); %solid propellant tank volume calculation, in m3
        vol_tank_ox_3 = m_ox_3/density_ox_3; %liquid propellant oxidizer tank volume calculation, in m3
         
            % Ullage effect
        
        vol_ullage_fuel_3 = 0;% since the fuel is solid, no ullage effect
        vol_ullage_ox_3 = vol_tank_ox_3*ullage_ox_3; % ullage effect for the liquid propellant, in m3
        
            % Shrinkage effect
        
        if cryo_ox_3 == 1 % shrinkage effect only happen sfor cryogenic propellant, oxidizer only considered since the fuel is solid 
            vol_shrinkage_ox_3 = vol_tank_ox_3*3*alpha_ox_3*delta_temp_ox_3; % volume of shrinkage effect, in m3
        else 
            vol_shrinkage_ox_3 = 0;
        end
        
        vol_shrinkage_fuel_3 = 0; % fuel is solid, so no shrinkage effect to consider
               
          
    elseif hybrid_3 == 0
        vol_tank_fuel_3 = m_fuel_3/density_fuel_3; %stage 3 oxidizer tank ideal volume, in m3
        vol_tank_ox_3 = m_ox_3/density_ox_3; %stage 3 oxidizer tank ideal volume, in m3

             % Ullage effect
        
        vol_ullage_ox_3 = vol_tank_ox_3*ullage_ox_3; % ullage effect volume for the oxidizer, in m3
        vol_ullage_fuel_3 = vol_tank_fuel_3*ullage_fuel_3;% ullage effect volume for the fuel, in m3
        
            % Shrinkage effect
        
        % Shrinkage effect only impactful for cryogenic propellants
        if cryo_ox_3 == 1
            vol_shrinkage_ox_3 = vol_tank_ox_3*3*alpha_ox_3*delta_temp_ox_3; % shrinkage volume for the oxidizer, in m3
        else 
            vol_shrinkage_ox_3 = 0;
        end
        
        if cryo_fuel_3 == 1
            vol_shrinkage_fuel_3 = vol_tank_fuel_3*3*alpha_fuel_3*delta_temp_fuel_3; % shrinkage volume for the fuel, in m3
        else 
            vol_shrinkage_fuel_3 = 0;
        end
        
        
    end
    
    vol_tank_ox_3_total = vol_tank_ox_3 + vol_shrinkage_ox_3 + vol_shrinkage_ox_3; % total oxidizer tank volume with the real-life added effects, in m3
    vol_tank_fuel_3_total = vol_tank_fuel_3 + vol_shrinkage_fuel_3 + vol_shrinkage_fuel_3; % total fuel tank volume with the real-life added effects, in m3

    % Stage 2 volume 
    
            % Ideal volume calculation
    if hybrid_2 == 1 % if the vehicle is hybrid, different volume calculation
        vol_tank_fuel_2 = 1.12*m_fuel_2/(density_fuel_2); % solid propellant tank volume calculation, in m3
        vol_tank_ox_2 = m_ox_2/density_ox_2; % liquid propellant oxidizer tank volume calculation, in m3
         
            % Ullage effect
        
        vol_ullage_fuel_2 = 0; % since the fuel is solid, no ullage effect
        vol_ullage_ox_2 = vol_tank_fuel_2*ullage_fuel_2; % ullage effect for the liquid propellant, in m3
        
            % Shrinkage effect
        
        if cryo_ox_2 == 1 % shrinkage effect only happens for cryogenic propellant, oxidizer only considered since the fuel is solid
            vol_shrinkage_ox_2 = vol_tank_ox_2*3*alpha_ox_2*delta_temp_ox_2; % volume of shrinkage effect, in m3
        else 
            vol_shrinkage_ox_2 = 0;
        end
        
        vol_shrinkage_fuel_2 = 0; % fuel is solid, so no shrinkage effect to consider
               
          
    elseif hybrid_2 == 0
        vol_tank_fuel_2 = m_fuel_2/density_fuel_2; % stage 2 fuel tank ideal volume, in m3
        vol_tank_ox_2 = m_ox_2/density_ox_2; %stage 2 oxidizer tank ideal volume, in m3

             % Ullage effect
        
        vol_ullage_ox_2 = vol_tank_ox_2*ullage_ox_2; % ullage effect volume for the oxidizer, in m3
        vol_ullage_fuel_2 = vol_tank_fuel_2*ullage_fuel_2; % ullage effect volume for the fuel, in m3
        
            % Shrinkage effect
        % Shrinkage effect only impactful for cryogenic propellants
        if cryo_ox_2 == 1 
            vol_shrinkage_ox_2 = vol_tank_ox_2*3*alpha_ox_2*delta_temp_ox_2; % shrinkage volume for the oxidizer, in m3
        else 
            vol_shrinkage_ox_2 = 0;
        end
        
        if cryo_fuel_2 == 1
            vol_shrinkage_fuel_2 = vol_tank_fuel_2*3*alpha_fuel_2*delta_temp_fuel_2; % shrinkage volume for the fuel, in m3
        else 
            vol_shrinkage_fuel_2 = 0;
        end
        
        
    end
    
    vol_tank_ox_2_total = vol_tank_ox_2 + vol_shrinkage_ox_2 + vol_shrinkage_ox_2; % total oxidizer tank volume with the real-life added effects, in m3
    vol_tank_fuel_2_total = vol_tank_fuel_2 + vol_shrinkage_fuel_2 + vol_shrinkage_fuel_2;% total fuel tank volume with the real-life added effects, in m3
      
    % Stage 1 volume 
    
            % Ideal volume calculation
    if hybrid_1 == 1 % if the vehicle is hybrid, different volume calculation
        vol_tank_fuel_1 = 1.12*m_fuel_1/(density_fuel_1); % solid propellant tank volume calculation, in m3
        vol_tank_ox_1 = m_ox_1/density_ox_1; % liquid propellant oxidizer tank volume calculation, in m3
         
            % Ullage effect
        
        vol_ullage_fuel_1 = 0;% since the fuel is solid, no ullage effect
        vol_ullage_ox_1 = vol_tank_fuel_1*ullage_fuel_1;% ullage effect for the liquid propellant, in m3
        
            % Shrinkage effect
        % shrinkage effect only happens for cryogenic propellant, oxidizer only considered since the fuel is solid
        if cryo_ox_3 == 1
            vol_shrinkage_ox_1 = vol_tank_ox_1*3*alpha_ox_1*delta_temp_ox_1;% volume of shrinkage effect, in m3
        else 
            vol_shrinkage_ox_1 = 0;
        end
        
        vol_shrinkage_fuel_1 = 0;% fuel is solid, so no shrinkage effect to consider
               
          
    elseif hybrid_1 == 0
        vol_tank_fuel_1 = m_fuel_1/density_fuel_1; %stage 1 fuel tank ideal volume, in m3
        vol_tank_ox_1 = m_ox_1/density_ox_1; %stage 1 oxidizer tank ideal volume, in m3

             % Ullage effect
        
        vol_ullage_ox_1 = vol_tank_ox_1*ullage_ox_1; % ullage effect volume for the oxidizer, in m3
        vol_ullage_fuel_1 = vol_tank_fuel_1*ullage_fuel_1; % ullage effect volume for the fuel, in m3
        
            % Shrinkage effect
        % Shrinkage effect only impactful for cryogenic propellants
        if cryo_ox_1 == 1
            vol_shrinkage_ox_1 = vol_tank_ox_1*3*alpha_ox_1*delta_temp_ox_1;% shrinkage volume for the oxidizer, in m3
        else 
            vol_shrinkage_ox_1 = 0;
        end
        
        if cryo_fuel_1 == 1
            vol_shrinkage_fuel_1 = vol_tank_fuel_1*3*alpha_fuel_1*delta_temp_fuel_1;% shrinkage volume for the fuel, in m3
        else 
            vol_shrinkage_fuel_1 = 0;
        end
        
        
    end
    
    vol_tank_ox_1_total = vol_tank_ox_1 + vol_shrinkage_ox_1 + vol_shrinkage_ox_1;% total oxidizer tank volume with the real-life added effects, in m3
    vol_tank_fuel_1_total = vol_tank_fuel_1 + vol_shrinkage_fuel_1 + vol_shrinkage_fuel_1;% total fuel tank volume with the real-life added effects, in m3
    
    % Stage 0 - Solid boosters tank volumes 
    
    if n_boosters > 0
    
        vol_tank_0_per_booster = 1.12*m_prop_0/(density_prop_0);% solid boosters casing volume, in m3
    else
    
        vol_tank_0_per_booster = 0;
    end
            
    %% Initial sizing calculations - dimensions of launcher
    
    %%%%% Tanks sizing 

    h_dome_ox_1 = (D_1/2)/AR_ox_1;% height of the oxidizer tank dome, in m
    h_cyl_ox_1 = 4*vol_tank_ox_1_total/(pi*D_1^2) - 2*D_1/(3*AR_ox_1); % height of the oxidizer tank cylinder, in m
    
    A_tank_ox_1 = pi*D_1*h_cyl_ox_1 + 2*1.6234*pi*(D_1/2)^2; % Area of the oxidizer tank, in m2
    if h_cyl_ox_1 < 0 % if the cylinder height is below zero, it means the two domes heights generate a high enough volume to store the desired propellant 
        h_cyl_ox_1 = 0;
    end
    
    h_cyl_fuel_1 = 4*vol_tank_fuel_1_total/(pi*D_1^2) - 2*D_1/(3*AR_fuel_1);% height of the fuel tank cylinder, in m
    h_dome_fuel_1 = (D_1/2)/AR_fuel_1;% height of the fuel tank dome, in m
    A_tank_fuel_1 = pi*D_1*h_cyl_fuel_1 + 2*1.6234*pi*(D_1/2)^2;% Area of the fuel tank, in m2
    if h_cyl_fuel_1 < 0 % if the cylinder height is below zero, it means the two domes heights generate a high enough volume to store the desired propellant 
        h_cyl_fuel_1 = 0;
    end
    
    h_cyl_ox_2 = 4*vol_tank_ox_2_total/(pi*D_2^2) - 2*D_2/(3*AR_ox_2);% height of the oxidizer tank cylinder, in m
    h_dome_ox_2 = (D_2/2)/AR_ox_2;% height of the oxidizer tank dome, in m
    A_tank_ox_2 = pi*D_2*h_cyl_ox_2 + 2*1.6234*pi*(D_2/2)^2;% Area of the oxidizer tank, in m2
    if h_cyl_ox_2 < 0 % if the cylinder height is below zero, it means the two domes heights generate a high enough volume to store the desired propellant 
        h_cyl_ox_2 = 0;
    end
    
    h_cyl_fuel_2 = 4*vol_tank_fuel_2_total/(pi*D_2^2) - 2*D_2/(3*AR_fuel_2);% height of the fuel tank cylinder, in m
    h_dome_fuel_2 = (D_2/2)/AR_fuel_2;% height of the fuel tank dome, in m
    A_tank_fuel_2 = pi*D_2*h_cyl_fuel_2 + 2*1.6234*pi*(D_2/2)^2;% Area of the fuel tank, in m2
    if h_cyl_fuel_2 < 0% if the cylinder height is below zero, it means the two domes heights generate a high enough volume to store the desired propellant 
        h_cyl_fuel_2 = 0;
    end
    
    h_cyl_ox_3 = 4*vol_tank_ox_3_total/(pi*D_3^2) - 2*D_3/(3*AR_ox_3);% height of the oxidizer tank cylinder, in m
    h_dome_ox_3 = (D_3/2)/AR_ox_3; % height of the oxidizer tank dome, in m
    A_tank_ox_3 = pi*D_3*h_cyl_ox_3 + 2*1.6234*pi*(D_3/2)^2; % Area of the oxidizer tank, in m2
    if h_cyl_ox_3 < 0 % if the cylinder height is below zero, it means the two domes heights generate a high enough volume to store the desired propellant 
        h_cyl_ox_3 = 0;
    end
    
    h_cyl_fuel_3 = 4*vol_tank_fuel_3_total/(pi*D_3^2) - 2*D_3/(3*AR_fuel_3); % height of the fuel tank cylinder, in m
    h_dome_fuel_3 = (D_3/2)/AR_fuel_3;% height of the fuel tank dome, in m
    A_tank_fuel_3 = pi*D_3*h_cyl_fuel_3 + 2*1.6234*pi*(D_3/2)^2;% Area of the fuel tank, in m2
    if h_cyl_fuel_2 < 0% if the cylinder height is below zero, it means the two domes heights generate a high enough volume to store the desired propellant 
        h_cyl_fuel_2 = 0;
    end
    
    
    %%%% Stage 1 dimensions
    
    L_1_aft_skirt = D_1; % after skirt length, in m
    L_1_intertank = 1/4*D_1 + h_dome_ox_1 + h_dome_fuel_1; % Intertank length, in m
    
    %%%% Stage 2 dimensions
    
    L_2_aft_skirt = 1/3*D_2 + h_dome_ox_2; % after skirt length, in m
    L_2_intertank = 1/4*D_2 + h_dome_ox_2 + h_dome_fuel_2; % Intertank length, in m
    
    %%%%% Stage 3 dimensions
    
    L_3_aft_skirt = 1/3*D_3 + h_dome_ox_3;% after skirt length, in m
    L_3_intertank = 1/4*D_3 + h_dome_ox_3 + h_dome_fuel_3; % Intertank length, in m
    
    if n_stages == 1
        
        L_1_interstage = 1/3*D_1 + h_dome_fuel_1; % first interstage length, in m
        L_2_interstage = 0;
        L_3_interstage = 0;
        L_payload_fairing = 2*D_1; % payload fairing assigned to first stage since there is only one stage, in m
    
        L_2_aft_skirt = 0;
        L_2_intertank = 0;
    
        L_3_aft_skirt = 0;
        L_3_intertank = 0;
    
        h_cyl_ox_2 = 0;
        h_dome_ox_2 = 0;
    
        h_cyl_fuel_2 = 0;
        h_dome_fuel_2 = 0;
        
        h_cyl_ox_3 = 0;
        h_dome_ox_3 = 0;
        
        h_cyl_fuel_3 = 0;
        h_dome_fuel_3 = 0;
    
    elseif n_stages == 2
    
        L_1_interstage = 5/4*D_1; % first interstage length, in m
        L_2_interstage = 1/3*D_2 + h_dome_fuel_2; % second interstage length, in m
        L_3_interstage = 0;
        L_payload_fairing = 2*D_2;% payload fairing assigned to second stage since there is two stages, in m
    
        L_3_aft_skirt = 0;
        L_3_intertank = 0; 
    
        h_cyl_ox_3 = 0;
        h_dome_ox_3 = 0;
        
        h_cyl_fuel_3 = 0;
        h_dome_fuel_3 = 0;
    
    elseif n_stages == 3
    
        L_1_interstage = 5/4*D_1; % first interstage length, in m
        L_2_interstage = 5/4*D_2; % second interstage length, in m
        L_3_interstage = 1/3*D_3 + h_dome_fuel_3;% third interstage length, in m
        L_payload_fairing = 2*D_3;% payload fairing assigned to second stage since there is three stages, in m
    end

    height_ox_1 = h_cyl_ox_1 + 2*h_dome_ox_1; % total height of the first stage oxidizer tank, in m
    height_ox_2 = h_cyl_ox_2 + 2*h_dome_ox_2;% total height of the second stage oxidizer tank, in m
    height_ox_3 = h_cyl_ox_3 + 2*h_dome_ox_3;% total height of the third stage oxidizer tank, in m

    height_fuel_1 = h_cyl_fuel_1 + 2*h_dome_fuel_1;% total height of the first stage fuel tank, in m
    height_fuel_2 = h_cyl_fuel_2 + 2*h_dome_fuel_2;% total height of the second stage fuel tank, in m
    height_fuel_3 = h_cyl_fuel_3 + 2*h_dome_fuel_3;% total height of the third stage fuel tank, in m

    height_stage1 = L_1_aft_skirt + h_cyl_ox_1 + h_cyl_fuel_1 + L_1_intertank + L_1_interstage; % total height of the first stage, in m
    
    height_stage2 = L_2_aft_skirt + h_cyl_ox_2 + h_cyl_fuel_2 + L_2_intertank + L_2_interstage + L_payload_fairing; % total height of the second stage, in m (necessity to assign the payload fairing to the last stage)
    
    height_stage3 = L_3_aft_skirt + h_cyl_ox_3 + h_cyl_fuel_3 + L_3_intertank + L_3_interstage;% total height of the third stage, in m (necessity to assign the payload fairing to the last stage)
    
    height_launcher = height_stage1 + height_stage2 + height_stage3 + 2*h_dome_FB_1 + h_cyl_FB_1 + 2*h_dome_FB_2 + h_cyl_FB_2; % total height of the launcher, in m
    
    
    %% Initial sizing calculations - mass breakdown calculation
    
    %%%%% Liquid tanks masses
    
    if cryo_ox_1 == 1 % different calculations for cryo and non cryo propellants
        m_tank_ox_1 = 9.09*vol_tank_ox_1_total; % mass of oxidizer tank, in m3
    else
        m_tank_ox_1 = 12.16*vol_tank_ox_1_total; % mass of oxidizer tank, in m3
    end
    
    if cryo_fuel_1 == 1 % different calculations for cryo and non cryo propellants
        m_tank_fuel_1 = 9.09*vol_tank_fuel_1_total; % mass of fuel tank, in m3
    else
        m_tank_fuel_1 = 12.16*vol_tank_fuel_1_total; % mass of fuel tank, in m3
    end
    
    
    if cryo_ox_2 == 1 % different calculations for cryo and non cryo propellants
        m_tank_ox_2 = 9.09*vol_tank_ox_2_total; % mass of oxidizer tank, in m3
    else
        m_tank_ox_2 = 12.16*vol_tank_ox_2_total; % mass of oxidizer tank, in m3
    end
    
    if cryo_fuel_2 == 1 % different calculations for cryo and non cryo propellants
        m_tank_fuel_2 = 9.09*vol_tank_fuel_2_total; % mass of fuel tank, in m3
    else
        m_tank_fuel_2 = 12.16*vol_tank_fuel_2_total; % mass of fuel tank, in m3
    end
    
    
    if cryo_ox_3 == 1 % different calculations for cryo and non cryo propellants
        m_tank_ox_3 = 9.09*vol_tank_ox_3_total; % mass of oxidizer tank, in m3
    else
        m_tank_ox_3 = 12.16*vol_tank_ox_3_total; % mass of oxidizer tank, in m3
    end
    
    if cryo_fuel_3 == 1 % different calculations for cryo and non cryo propellants
        m_tank_fuel_3 = 9.09*vol_tank_fuel_3_total; % mass of fuel tank, in m3
    else
        m_tank_fuel_3 = 12.16*vol_tank_fuel_3_total; % mass of fuel tank, in m3
    end
    
    %%%% Solid casing masses
    % Sorting out the cases were solid propellant are applicable

    if hybrid_3 == 1
       m_tank_fuel_3 = 0.135*m_fuel_3;
    end

    if hybrid_2 == 1
    m_tank_fuel_2 = 0.135*m_fuel_2;
    end

    if hybrid_1 == 1
    m_tank_fuel_1 = 0.135*m_fuel_1;
    end

   if n_boosters > 0
        
        m_casing_per_booster = 0.135 * m_prop_0;
    
    else
    
        m_casing_per_booster = 0;
    end
    
    
    %%%% Avionics masses
    
    %distribution of the avionics mass depending on the amount of stage, to put 80% of the mass on the upper stage 
    
    if n_stages == 1
        m_avionics_1 = 350;
        m_avionics_2 = 0;
        m_avionics_3 = 0;
        m_avionics_total = m_avionics_1;
    elseif n_stages == 2
        m_avionics_1 = 0.2*350;
        m_avionics_2 = 0.8*350;
        m_avionics_3 = 0;
        m_avionics_total = m_avionics_1 + m_avionics_2;
    elseif n_stages == 3
        m_avionics_1 = 0.1*350;
        m_avionics_2 = 0.1*350;
        m_avionics_3 = 0.8*350;
        m_avionics_total = m_avionics_1 + m_avionics_2 + m_avionics_3;
    end
    
    
    %%%% Cryogenic tank insulation mass
    %Tank insulation is only needed for cryo propellants and applied in that case only
    if cryo_ox_1 == 1
        m_tank_ox_insulation_1 = 2*A_tank_ox_1; % Usage of the areal density and the area of the tank, mass in kg
    else
        m_tank_ox_insulation_1 = 0;
    end
    
    if cryo_fuel_1 == 1
        m_tank_fuel_insulation_1 = 2*A_tank_fuel_1; % Usage of the areal density and the area of the tank, mass in kg
    else
        m_tank_fuel_insulation_1 = 0;
    end
    
    if cryo_ox_2 == 1
        m_tank_ox_insulation_2 = 2*A_tank_ox_2; % Usage of the areal density and the area of the tank, mass in kg
    else
        m_tank_ox_insulation_2 = 0;
    end
    
    if cryo_fuel_2 == 1
        m_tank_fuel_insulation_2 = 2*A_tank_fuel_2; % Usage of the areal density and the area of the tank, mass in kg
    else
        m_tank_fuel_insulation_2 = 0;
    end
    
    if cryo_ox_3 == 1
        m_tank_ox_insulation_3 = 2*A_tank_ox_3; % Usage of the areal density and the area of the tank, mass in kg
    else
        m_tank_ox_insulation_3 = 0;
    end
    
    if cryo_fuel_3 == 1
        m_tank_fuel_insulation_3 = 2*A_tank_fuel_3; % Usage of the areal density and the area of the tank, mass in kg
    else
        m_tank_fuel_insulation_3 = 0;
    end
    
    
    %%%%% Masses of skirts, intertanks, interfairing and fairing
    
    %%% Interstages and intertanks and after skirts
    
    if n_stages == 3
        A_aft_skirt_1 = 2*pi*D_1*L_1_aft_skirt; %Area of afterskirt between first and second stage, in m2
        A_aft_skirt_2 = 2*pi*D_2*L_2_aft_skirt; %Area of afterskirt between second and third stage, in m2
        A_aft_skirt_3 = 2*pi*D_3*L_3_aft_skirt; %Area of afterskirt between third and fairing stage, in m2
    
        %Usage of areal density to calculate the mass of the after skirts, depending on the metallic or composite material chosen
        m_aft_skirt_1_metal = 13.3*A_aft_skirt_1;
        m_aft_skirt_1_composite = 9.89*A_aft_skirt_1;
        
        m_aft_skirt_2_metal = 13.3*A_aft_skirt_2;
        m_aft_skirt_2_composite = 9.89*A_aft_skirt_2;
       
        m_aft_skirt_3_metal = 13.3*A_aft_skirt_3;
        m_aft_skirt_3_composite = 9.89*A_aft_skirt_3;
    
        A_interstage_1 = 2*pi*D_1*L_1_interstage; %Area of interstage between first and second stage, in m2
        A_interstage_2 = 2*pi*D_2*L_2_interstage; %Area of interstage between second and third stage, in m2
        A_interstage_3 = 2*pi*D_3*L_3_interstage; %Area of interstage between third and fairing stage, in m2
        
        %Usage of areal density to calculate the mass of the after skirts, depending on the metallic or composite material chosen
        m_interstage_1_metal = 13.3*A_interstage_1;
        m_interstage_1_composite = 9.89*A_interstage_1;
        
        m_interstage_2_metal = 13.3*A_interstage_2;
        m_interstage_2_composite = 9.89*A_interstage_2;
        
        m_interstage_3_metal = 13.3*A_interstage_3;
        m_interstage_3_composite = 9.89*A_interstage_3;
        
        A_intertank_1 = 2*pi*D_1*L_1_intertank;%Area of intertank between oxidizer and fuel tank in first stage, in m2
        A_intertank_2 = 2*pi*D_2*L_2_intertank;%Area of intertank between oxidizer and fuel tank in second stage, in m2
        A_intertank_3 = 2*pi*D_3*L_3_intertank;%Area of intertank between oxidizer and fuel tank in third stage, in m2
        
        %Usage of areal density to calculate the mass of the after skirts, depending on the metallic or composite material chosen
        m_intertank_1_metal = 13.3*A_intertank_1;
        m_intertank_1_composite = 9.89*A_intertank_1;
    
        m_intertank_2_metal = 13.3*A_intertank_2;
        m_intertank_2_composite = 9.89*A_intertank_2;
    
        m_intertank_3_metal = 13.3*A_intertank_3;
        m_intertank_3_composite = 9.89*A_intertank_3;
    
    elseif n_stages == 2
        A_aft_skirt_1 = 2*pi*D_1*L_1_aft_skirt; %Area of afterskirt between first and second stage, in m2
        A_aft_skirt_2 = 2*pi*D_2*L_2_aft_skirt; %Area of afterskirt between second and third stage, in m2
        A_aft_skirt_3 = 0; %Area of afterskirt between third and fairing stage, in m2
        
        %Usage of areal density to calculate the mass of the after skirts, depending on the metallic or composite material chosen
        m_aft_skirt_1_metal = 13.3*A_aft_skirt_1;
        m_aft_skirt_1_composite = 9.89*A_aft_skirt_1;
        
        m_aft_skirt_2_metal = 13.3*A_aft_skirt_2;
        m_aft_skirt_2_composite = 9.89*A_aft_skirt_2;
       
        m_aft_skirt_3_metal = 13.3*A_aft_skirt_3;
        m_aft_skirt_3_composite = 9.89*A_aft_skirt_3;
    
        A_interstage_1 = 2*pi*D_1*L_1_interstage;%Area of interstage between first and second stage, in m2
        A_interstage_2 = 2*pi*D_2*L_2_interstage;%Area of interstage between second and third stage, in m2
        A_interstage_3 = 0;%Area of intertank between oxidizer and fuel tank in third stage, in m2
        
        %Usage of areal density to calculate the mass of the after skirts, depending on the metallic or composite material chosen
        m_interstage_1_metal = 13.3*A_interstage_1;
        m_interstage_1_composite = 9.89*A_interstage_1;
        
        m_interstage_2_metal = 13.3*A_interstage_2;
        m_interstage_2_composite = 9.89*A_interstage_2;
        
        m_interstage_3_metal = 13.3*A_interstage_3;
        m_interstage_3_composite = 9.89*A_interstage_3;
            
        A_intertank_1 = 2*pi*D_1*L_1_intertank;%Area of intertank between oxidizer and fuel tank in first stage, in m2
        A_intertank_2 = 2*pi*D_2*L_2_intertank;%Area of intertank between oxidizer and fuel tank in second stage, in m2
        A_intertank_3 = 0;%Area of intertank between oxidizer and fuel tank in third stage, in m2
        
        %Usage of areal density to calculate the mass of the after skirts, depending on the metallic or composite material chosen
        m_intertank_1_metal = 13.3*A_intertank_1;
        m_intertank_1_composite = 9.89*A_intertank_1;
    
        m_intertank_2_metal = 13.3*A_intertank_2;
        m_intertank_2_composite = 9.89*A_intertank_2;
    
        m_intertank_3_metal = 13.3*A_intertank_3;
        m_intertank_3_composite = 9.89*A_intertank_3;
    
    elseif n_stages == 1
        A_aft_skirt_1 = 2*pi*D_1*L_1_aft_skirt; %Area of afterskirt between first and second stage, in m2
        A_aft_skirt_2 = 0; %Area of afterskirt between second and third stage, in m2
        A_aft_skirt_3 = 0; %Area of afterskirt between third and fairing stage, in m2
        
        %Usage of areal density to calculate the mass of the after skirts, depending on the metallic or composite material chosen
        m_aft_skirt_1_metal = 13.3*A_aft_skirt_1;
        m_aft_skirt_1_composite = 9.89*A_aft_skirt_1;
        
        m_aft_skirt_2_metal = 13.3*A_aft_skirt_2;
        m_aft_skirt_2_composite = 9.89*A_aft_skirt_2;
       
        m_aft_skirt_3_metal = 13.3*A_aft_skirt_3;
        m_aft_skirt_3_composite = 9.89*A_aft_skirt_3;
    
        A_interstage_1 = 2*pi*D_1*L_1_interstage;%Area of interstage between first and second stage, in m2
        A_interstage_2 = 0;%Area of interstage between second and third stage, in m2
        A_interstage_3 = 0;%Area of intertank between oxidizer and fuel tank in third stage, in m2
        
        %Usage of areal density to calculate the mass of the after skirts, depending on the metallic or composite material chosen
        m_interstage_1_metal = 13.3*A_interstage_1;
        m_interstage_1_composite = 9.89*A_interstage_1;
        
        m_interstage_2_metal = 13.3*A_interstage_2;
        m_interstage_2_composite = 9.89*A_interstage_2;
        
        m_interstage_3_metal = 13.3*A_interstage_3;
        m_interstage_3_composite = 9.89*A_interstage_3;
    
        A_intertank_1 = 2*pi*D_1*L_1_intertank;%Area of intertank between oxidizer and fuel tank in first stage, in m2
        A_intertank_2 = 0;%Area of intertank between oxidizer and fuel tank in second stage, in m2
        A_intertank_3 = 0;%Area of intertank between oxidizer and fuel tank in third stage, in m
        
        %Usage of areal density to calculate the mass of the after skirts, depending on the metallic or composite material chosen
        m_intertank_1_metal = 13.3*A_intertank_1;
        m_intertank_1_composite = 9.89*A_intertank_1;
    
        m_intertank_2_metal = 13.3*A_intertank_2;
        m_intertank_2_composite = 9.89*A_intertank_2;
    
        m_intertank_3_metal = 13.3*A_intertank_3;
        m_intertank_3_composite = 9.89*A_intertank_3;
    end
    
    
    %%%% Fairing mass
    
    % calculation of the fairing mass in kg based on the number of stages, to know where to attach the fairing (the most upper stage)    
    if n_stages == 3
    
        m_fairing = pi*(D_3/2)*sqrt((D_3/2)^2+L_payload_fairing^2);
    
    elseif n_stages == 2
        
        m_fairing = pi*(D_2/2)*sqrt((D_2/2)^2+L_payload_fairing^2);
    
    elseif n_stages == 1
    
        m_fairing = pi*(D_1/2)*sqrt((D_1/2)^2+L_payload_fairing^2);
    
    end
    
    
    %%%%% Masses of liquid rocket engine and thrust structure and definition of number of engines
    
    if n_stages == 3
        
        m_engine_1 = Thrust_1*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20)) + 59; % mass of the engine, in kg
        m_thrust_struct_1 = 2.55*10^(-4)*Thrust_1; % mass of the thrust structure, in kg
               
        m_engine_2 = Thrust_2*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20)) + 59; % mass of the engine, in kg
        m_thrust_struct_2 = 2.55*10^(-4)*Thrust_2;% mass of the thrust structure, in kg
    
        m_engine_3 = Thrust_3*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20)) + 59;% mass of the engine, in kg
        m_thrust_struct_3 = 2.55*10^(-4)*Thrust_3;% mass of the thrust structure, in kg
        
        % if the engine is hybrid, application of a mass reduction coefficient due to the difference in components
        if hybrid_3 == 1
            m_engine_3 = 0.68*m_engine_3;
        end

        if hybrid_2 == 1
            m_engine_2 = 0.68*m_engine_2;
        end

        if hybrid_1 == 1
            m_engine_1 = 0.68*m_engine_1;
        end

        % Calculation of the number of engines based on the needed thrust to weight ratio and therefore the amount of engines needed to achieve it 
        if n_engines_1*Thrust_1/(m_lift_off(1)*g0) < T_W_1
            while abs(n_engines_1*Thrust_1/(m_lift_off(1)*g0) - T_W_1) > 10^(-1)
            
            n_engines_1 = n_engines_1 + 0.1;
    
            end
        end

    
        if n_engines_2*Thrust_2/((m_lift_off(1) - (m_struct_1_stored(1) + m_prop_1_stored(1)))*g0) < T_W_2
            while abs(n_engines_2*Thrust_2/((m_lift_off(1) - (m_struct_1_stored(1) + m_prop_1_stored(1)))*g0) - T_W_2) > 10^(-1)
            
            abs(Thrust_2/((m_lift_off(1) - (m_struct_1_stored(1) + m_prop_1_stored(1)))*g0) - T_W_2)
            n_engines_2 = n_engines_2 + 0.1;
            
            end
        end

        
        if n_engines_3*Thrust_3/((m_lift_off(1) - (m_struct_1_stored(1) + m_prop_1_stored(1)) - (m_struct_2_stored(1) + m_prop_2_stored(1)))*g0) < T_W_3

            while abs(n_engines_3*Thrust_3/((m_lift_off(1) - (m_struct_1_stored(1) + m_prop_1_stored(1)) - (m_struct_2_stored(1) + m_prop_2_stored(1)))*g0) - T_W_3) > 10^(-1)
            
            n_engines_3 = n_engines_3 + 0.1;
            
            end
        end

    
    elseif n_stages == 2
    
        m_engine_1 = Thrust_1*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20)) + 59;% mass of the engine, in kg
        m_thrust_struct_1 = 2.55*10^(-4)*Thrust_1;% mass of the thrust structure, in kg
               
        m_engine_2 = Thrust_2*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20)) + 59;% mass of the engine, in kg
        m_thrust_struct_2 = 2.55*10^(-4)*Thrust_2;% mass of the thrust structure, in kg

        m_engine_3 = 0;
        n_engines_3 = 0;

        % if the engine is hybrid, application of a mass reduction coefficient due to the difference in components
        if hybrid_2 == 1
        m_engine_2 = 0.68*m_engine_2;
        end

        if hybrid_1 == 1
        m_engine_1 = 0.68*m_engine_1;
        end
        
        % Calculation of the number of engines based on the needed thrust to weight ratio and therefore the amount of engines needed to achieve it
        if n_engines_1*Thrust_1/(m_lift_off(1)*g0) < T_W_1

            while abs(n_engines_1*Thrust_1/(m_lift_off(1)*g0) - T_W_1) > 10^(-1)

            n_engines_1 = n_engines_1 + 0.1;
    
            end
        end

        if n_engines_2*Thrust_2/((m_lift_off(1) - (m_struct_1_stored(1) + m_prop_1_stored(1)))*g0) < T_W_2
            while abs(n_engines_2*Thrust_2/((m_lift_off(1) - (m_struct_1_stored(1) + m_prop_1_stored(1)))*g0) - T_W_2) > 10^(-1)
            
            n_engines_2 = n_engines_2 + 0.01;
    
            end
        end
                
    
    elseif n_stages == 1
        
        m_engine_1 = Thrust_1*(7.81*10^(-4) + 3.37*10^(-5)*sqrt(20)) + 59;% mass of the engine, in kg
        m_thrust_struct_1 = 2.55*10^(-4)*Thrust_1;% mass of the thrust structure, in kg
        
        % if the engine is hybrid, application of a mass reduction coefficient due to the difference in components
        if hybrid_1 == 1
        m_engine_1 = 0.68*m_engine_1;
        end
        
        % Calculation of the number of engines based on the needed thrust to weight ratio and therefore the amount of engines needed to achieve it
        if n_engines_1*Thrust_1/(m_lift_off(1)*g0) < T_W_1
            while abs(n_engines_1*Thrust_1/(m_lift_off(1)*g0) - T_W_1) > 10^(-1)
    
            n_engines_1 = n_engines_1 + 0.1;
    
            end
        end


        m_engine_2 = 0;

        m_engine_3 = 0;
    
    end

    %%%% Mass of gimbals 

    m_gimbals_1 = 237.8*(Thrust_1*n_engines_1/Pressure_1)^(0.9375);% mass of gimbals, in kg
    m_gimbals_2 = 237.8*(Thrust_2*n_engines_2/Pressure_2)^(0.9375);% mass of gimbals, in kg
    m_gimbals_3 = 237.8*(Thrust_3*n_engines_3/Pressure_3)^(0.9375);% mass of gimbals, in kg

    %%%%% Mass of payloach attach fitting
    
    m_PAF = 0.0755*m_payload + 50;% PAF mass in kg
    
    
    %%%%% Mass of electrical wiring
    
    m_wiring = 1.43*height_launcher; % mass of wiring for the whole launcher, in kg
    
    
    %%%%% Comparison between structure mass and mass breakdown calculations
    
    if n_stages == 3
        m_struct_1_breakdown_sum = m_prop_residual_1 + m_gimbals_1 + m_wiring/n_stages + n_engines_1*m_engine_1 + m_intertank_1_metal + m_aft_skirt_1_metal + m_interstage_1_metal + n_engines_1*m_thrust_struct_1 + m_tank_fuel_1 + m_tank_ox_1 + m_tank_ox_insulation_1 + m_tank_fuel_insulation_1 + m_avionics_1;
        m_struct_2_breakdown_sum = m_prop_residual_2 + m_gimbals_2 + m_wiring/n_stages + n_engines_2*m_engine_2 + m_intertank_2_metal + m_aft_skirt_2_metal + m_interstage_2_metal + n_engines_2*m_thrust_struct_2 + m_tank_fuel_2 + m_tank_ox_2 + m_tank_ox_insulation_2 + m_tank_fuel_insulation_2 + m_avionics_2;
        m_struct_3_breakdown_sum = m_prop_residual_3 + m_gimbals_3 + m_wiring/n_stages + n_engines_3*m_engine_3 + m_intertank_3_metal + + m_aft_skirt_3_metal + m_interstage_3_metal + n_engines_3*m_thrust_struct_3 + m_tank_fuel_3 + m_tank_ox_3 + m_tank_ox_insulation_3 + m_tank_fuel_insulation_3 + m_avionics_3 + m_PAF;
    
    elseif n_stages == 2
        m_struct_1_breakdown_sum = m_prop_residual_1 + m_gimbals_1 + m_wiring/n_stages + n_engines_1*m_engine_1 + m_intertank_1_metal +  m_aft_skirt_1_metal + m_interstage_1_metal + n_engines_1*m_thrust_struct_1 + m_tank_fuel_1 + m_tank_ox_1 + m_tank_ox_insulation_1 + m_tank_fuel_insulation_1 + m_avionics_1;
        m_struct_2_breakdown_sum = m_prop_residual_2 + m_gimbals_2 + m_wiring/n_stages + n_engines_2*m_engine_2 + m_intertank_2_metal +  m_aft_skirt_2_metal  + m_interstage_2_metal + n_engines_2*m_thrust_struct_2 + m_tank_fuel_2 + m_tank_ox_2 + m_tank_ox_insulation_2 + m_tank_fuel_insulation_2 + m_avionics_2 + m_PAF;
        m_struct_3_breakdown_sum = 0;
    
    elseif n_stages == 1
        m_struct_1_breakdown_sum = m_prop_residual_1 + m_gimbals_1 + m_wiring/n_stages + n_engines_1*m_engine_1 + m_intertank_1_metal +  m_aft_skirt_1_metal + m_interstage_1_metal + n_engines_1*m_thrust_struct_1 + m_tank_fuel_1 + m_tank_ox_1 + m_tank_ox_insulation_1 + m_tank_fuel_insulation_1 + m_avionics_1 + m_PAF;
        m_struct_2_breakdown_sum = 0;
        m_struct_3_breakdown_sum = 0;
    end
    
    
    
    % Iteration of the program to find the payload mass reduction due to the added reusability mass -> objective is to converge towards the same expendable GLOM with te added reusability mass by reducing the payload
    if count >= 2
        if abs(m_lift_off(count) - m_lift_off(1)) > 5000

            m_payload = m_payload - (count-0.1);
            count = count + 1;
                        
        else
            index = 2;
        end
    else
        count = count + 1;
        
    end
end

n_engines_1 = ceil(n_engines_1);
n_engines_2 = ceil(n_engines_2);
n_engines_3 = ceil(n_engines_3);

struct_factor_1 = m_struct_1/(m_struct_1 + m_prop_1);
struct_factor_2 = m_struct_2/(m_struct_2 + m_prop_2);

fprintf('height stage 1 = %d \n\n', height_stage1)
fprintf('height ox tank 1 = %d \n\n', height_ox_1)
fprintf('height fuel tank 1 = %d \n\n', height_fuel_1)
fprintf('number of main engines 1 = %d \n\n', n_engines_1)
fprintf('main engine dry mass 1 = %d \n\n', m_engine_1)
fprintf('airbreathing engine dry mass 1 = %d \n\n', m_engine_FB_1)

if VTHL_1_FB == 1
    fprintf('number of airbreathing engines 1 = %d \n\n', Thrust_FB_1/Thrust_per_engine_FB_1)
    fprintf('height added tank FB 1 = %d \n\n', h_cyl_FB_1 + 2*h_dome_FB_1)
end

fprintf('stage dry mass 1 = %d \n\n', m_struct_1)
fprintf('stage oxidizer mass 1 = %d \n\n', m_ox_1)
fprintf('stage fuel mass 1 = %d \n\n', m_fuel_1)
fprintf('stage propellant mass 1 = %d \n\n', m_prop_1)

fprintf('height stage 2 = %d \n\n', height_stage2)
fprintf('height ox tank 2 = %d \n\n', height_ox_2)
fprintf('height fuel tank 2 = %d \n\n', height_fuel_2)
fprintf('number of main engines 2 = %d \n\n', n_engines_2)
fprintf('main engine dry mass 2 = %d \n\n', m_engine_2)
fprintf('airbreathing engine dry mass 2 = %d \n\n', m_engine_FB_2)
if VTHL_2_FB == 1
    fprintf('number of airbreathing engines 1 = %d \n\n', Thrust_FB_2/Thrust_per_engine_FB_2)
    
end
fprintf('stage dry mass 2 = %d \n\n', m_struct_2)
fprintf('stage oxidizer mass 2 = %d \n\n', m_ox_2)
fprintf('stage fuel mass 2 = %d \n\n', m_fuel_2)
fprintf('stage propellant mass 2 = %d \n\n', m_prop_2)

fprintf('GLOM mass = %d \n\n', m_prop_1 + m_struct_1 + m_prop_2 + m_struct_2 + m_payload)
fprintf('launcher total height = %d \n\n', height_launcher)
fprintf('launcher dry mass = %d \n\n', m_struct_1 + m_struct_2)
fprintf('launcher propellant mass = %d \n\n', m_prop_1 + m_prop_2)

fprintf('struct factor 1 = %d \n\n', struct_factor_1)
fprintf('struct factor 2 = %d \n\n', struct_factor_2)
fprintf('payload mass = %d \n\n', m_payload)
