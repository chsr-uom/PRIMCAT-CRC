<img src="PRIMCAT-logo-small.jpg" align="right"/> <br>
This repository contains the R code associated with our manuscript https://doi.org/10.1016/j.jval.2024.06.006 and is part of the <a href="https://mdhs.unimelb.edu.au/centre-for-cancer-research/flagships/primcat-predicting-the-population-health-economic-impact-of-current-and-new-cancer-treatments" target="_blank">PRIMCAT research programme.</a>

### Predicting the Population Health Economic Impact of Current and New Cancer Treatments for Colorectal Cancer: A Data-Driven Whole Disease Simulation Model for Predicting the Number of Patients with Colorectal Cancer by Stage and Treatment Line in Australia.


------------------------------------------------------------------------

#### 📖 Abstract

**Objectives**: Effective healthcare planning, resource allocation, and budgeting require accurate predictions of the number of patients needing treatment at specific cancer stages and treatment lines. The Predicting the Population Health Economic Impact of Current and New Cancer Treatments (PRIMCAT) for Colorectal Cancer (CRC) simulation model (PRIMCAT-CRC) was developed to meet this requirement for all CRC stages and relevant molecular profiles in Australia.

**Methods**: Real-world data were used to estimate treatment utilization and time-to-event distributions. This populated a discrete-event simulation, projecting the number of patients receiving treatment across all disease stages and treatment lines for CRC and forecasting the number of patients likely to utilize future treatments. Illustrative analyses were undertaken, estimating treatments across disease stages and treatment lines over a 5-year period (2022-2026). We demonstrated the model s applicability through a case study introducing pembrolizumab as a first-line treatment for mismatch-repair-deficient stage IV.

**Results**: Clinical registry data from 7163 patients informed the model. The model forecasts 15 738 incident and 2821 prevalent cases requiring treatment in 2022, rising to 15 921 and 2871, respectively, by 2026. Projections show that over 2022 to 2026, there will be a total of 116 752 treatments initiated, with 43% intended for stage IV disease. The introduction of pembrolizumab is projected for 706 patients annually, totaling 3530 individuals starting treatment with pembrolizumab over the forecasted period, without significantly altering downstream utilization of subsequent treatments.

**Conclusion**: PRIMCAT-CRC is a versatile tool that can be used to estimate the eligible patient populations for novel cancer therapies, thereby reducing
uncertainty for policymakers in decisions to publicly reimburse new treatments.   

------------------------------------------------------------------------

#### 🔍 Repository content

``` bash
├─ Input                               # Input files required to run PRIMCAT-CRC
│  ├─ AIHW_CRC_Incidence.xlsx          # Incidence CRC in Australia 2010-2021
│  ├─ CVDL_CRC_StageDistribution.xlsx  # Stage distribution CRC from CVDL 2010-2019
│  ├─ pars_deterministic.RDS           # Estimated Parameters from ACCORD & TRACC analyses
│  └─ DES_function.R                   # R script with DES functions 
│
├─ Script
│  ├─ 1_PRIMCAT_CRC_function.R         # R script defining the PRIMCAT-CRC model
│  ├─ 2_model_validation.R             # R script used for model validation
│  ├─ 3_base_case.R                    # Base case scenario analysis
│  ├─ 4_scenario_pembro.R              # Scenario analysis for Pembrolizumab introduction
│  └─ 5_manuscript.R                   # Code used to produce manuscript tables and figures
│
├─ PRIMCAT-logo-small.jpg
├─ primcat_collab.png
├─ LICENCE
└─ README.md
```

------------------------------------------------------------------------

#### 🔧 What are the functions?
There are several functions used to run PRIMCAT-CRC. We provide an overview of these.

The following functions are defined:
- PRIMCAT_CRC_Parameters      function that defines a parameters object to run the model
- PRIMCAT_CRC_Population      function that defines the population to be simulation
- PRIMCAT_CRC_Simulation      function that contains the actual DES
- PRIMCAT_CRC                 function that calls the above functions multiple times in parallel to perform
                              multiple runs of the DES and average the output

1. The PRIMCAT_CRC_Parameters function loads the parameters that have been exported from the 1_parameters.R script
and returns them in a list that can be used to run the simulation. A list of parameters is defined for each time
point at which the parameters are to be defined/updated. The deterministic basecase values are used to define ALL
parameters at the start time, which is specified through the "starttime" argument. These values can be updated
through the list that can be passed to the "updates" argument. This way, parameter value changes over time can be
incorporated. The list that is returned will be used in the simulation to define events at the times defined by
the names of the pars list, at which the parameters will be loaded/updated. The "updates" argument can be used,
for example, to define the scenarios for the horizon scanning. The "updates" arguments requires a list of list,
each of which should be named according to the time at which the updates are to be performed, with the list
itself being a named list of parameters to be updated and their values (which can again be a list, for example to
define a distribution). Note that parameters are not changed for a period of time, but from the moment at which
they are updated for the remainder of the simulation or until the are updated again. an example of the "update"
list is as follows:
  list('1' = list(p_C1_SG = 0.5,p_C1_NEW = 0.5,p_C1_NEW_recur = 0.05,), '2' = list(p_C1_SG = 0, p_C1_NEW = 1))
This example updates the three parameters at time 1, and two of the three parameters again at time 2. The
p_C1_NEW_recur parameter that is updated only once, remains at it's updated value of 0.05.

2. The PRIMCAT_CRC_Population function defines the patient population that is to be simulated. For each year that is
to be simulated, it defines the number of individuals and their disease stage at diagnosis and at what time they
enter the simulation. The resulting data.frame can be used to run the simulation. All the information about the
population characteristics for the years that are to be simulated are provided through the 'm_years' argument.
This function input m_years should be a matrix or data.frame with a row for each year for which a population is
to be generated.

If function input by_total_incidence is TRUE, then the population is defined by applying a stage distribution
based on probability to the total incidence, for which the following columns are required:
- t       the simulation time at which the year should start
- n       the total incidence of that year
- p_C1    the probability of being C1 in that year
- P_C2    the probability of being C2 in that year
- p_C3    the probability of being C3 in that year
- p_R1    the probability of being R1 in that year
- P_R2    the probability of being R2 in that year
- p_R3    the probability of being R3 in that year
- p_CR4   the probability of being CR4 in that year

If function input by_total_incidence is FALSE, then the population is defined by stage-specific incidences using
the following columns:
- t       the simulation time at which the year should start
- n_C1    the incidence of C1 in that year
- n_C2    the incidence of C2 in that year
- n_C3    the incidence of C3 in that year
- n_R1    the incidence of R1 in that year
- n_R2    the incidence of R2 in that year
- n_R3    the incidence of R3 in that year
- n_CR4   the incidence of CR4 in that year

3. The PRIMCAT_CRC_Simulation function is the actual DES that simulated the individual patient trajectories. Time in
the simulation is defined in years It has the following input arguments:
- start     data.frame or matrix that specifies the number of individuals (rows) and must contains a startTime
            column that specifies the time at which the individual enters the simulation and a firstEvent column
            that specifies the first event that is to be simulated for the individual. Such a data.frame is
            returned by the PRIMCAT_CRC_Population function
- pars      list containing all the parameters used in the simulation, which is returned from the
            PRIMCAT_CRC_Parameters function
- until     numeric defining the time until which the simulation is to be run, with the default being Inf
            resulting in a simulation until all events have occurred
- seed      random number seed on for reproducibility

The function has three main parts:
  1. INITIALIZATION     loading the general DES functions and defining all the objects required in the simulation
  2. EVENT HANDLING     processing all the events that define the simulation model structure
  3. SIMULATION         running the simulation

It returns the output of the runAll() function, which is a list containing both the eventLog and resourceLog,
which contain a individual-level record of the events occurred and resources used, respectively.

1. INITIALIZATION ----
  Loading the general DES functions:
  - draw            randomly samples an event/index from a vector of probabilities
  - rcustom         convenience function to perform a random draw from any distribution
  - rmixture        custom function to sample from a mixture distribution
  - rgompertz       link to the flexsurv rgompertz function
  - rweibullPH      link to the flexsurv weibullPH function
  - priorityQueue   defines a event queue as backbone for a discrete event simulation
  - schedule        adds an event to the event queue
  - seize           seizes a resource by recording the start and end time
  - runOne          runs the discrete event simulation for one individual
  - runAll          runs the complete discrete event simulation for all individuals

2. EVENT HANDLING ----

  The goal of the init function is to ensure all global variable are reset or set to the value specified for that individual by the start input matrix/data.frame, as well as to schedule the first event and the 'Until'
  event based on the until argument. The handleEvent function defines what actions have to be performed when a particular event occurs. It is
  ordered according to the different disease stages and treatment lines. A switch statement is used to identify
  the event that needs to be processed, because this is much faster than a large if-elseif-elseif-elseif...
  statement. The names of the events correspond to the names of the treatment lines etc. in the conceptual model.
  Note that extensive comments are provided for the first events and for subsequent events when relevant. To
  reduce the length of the code, comments are not repeated too much.

4. The PRIMCAT_CRC function runs the simulation multiple times and summarizes and averages the number of events and
resource utilization in certain time periods. It has the following input arguments:
- n_runs                the number of runs to be performed, each using a different random number seed
- n_nodes               the number of CPU cores to be used to perform the runs of the simulation in parallel
- m_years               see the PRIMCAT_CRC_Population function for more information
- timegrid              numerical vector defining the breaks of the time periods for which the events and use of
                        resources is to be summarized. If NULL (default), year quarters are used corresponding to
                        four time periods per year
- gridnames             vector of names for the time periods specified through the timegrid argument. If NULL
                        (default), the quarters are named as Q1, Q2 etc for each year
- by_total_incidence    see the PRIMCAT_CRC_Population function for more information
- until                 see the PRIMCAT_CRC_Simulation function for more information

------------------------------------------------------------------------

#### 📫 Help needed?

Please let us know if you need help with PRIMCAT-CRC.

<img src="primcat_collab.png" align="center"/>

