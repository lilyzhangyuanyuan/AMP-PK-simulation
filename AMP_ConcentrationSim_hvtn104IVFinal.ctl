;################################################################################
;#The code of this paper includes three files:
;# (1)"dataShell_AMP_sim.R" : simulating the infusion and non-infusion visits 
;#                            according to the AMP protocol specification with 
;#                            different adherences and sample sizes.
;# (2)"AMP_ConcentrationSim_hvtn104IVFinal.ctl": simulating concentrations of AMP
;#                                               participants by using the master
;#                                               popPK model, which is the model 
;#                                               that best describes the HVTN104
;#                                               data.
;# (3)"AMP_estimate_PK_parameters": estimating popPK model parameters based on 
;#                                  the simulated concentration of AMP 
;#                                  participants.
;#
;################################################################################

$PROB 1
$INPUT ID TIME AMT RATE DV WT AGE sex f_vis dose visit_n boot
$DATA ../1.csv IGNORE=#
$SUBROUTINES ADVAN3 TRANS4
$PK

   TVCL=THETA(1)*EXP(THETA(5)*(WT-74.5))
   CL=TVCL*EXP(ETA(1))

   TVQ=THETA(2)
   Q=TVQ*EXP(ETA(2))

   TVV2=THETA(3)
   V2=TVV2*EXP(ETA(3))
   
   TVV1=THETA(4)*EXP(THETA(6)*(WT-74.5))
   V1=TVV1

   S1=V1

$ERROR
   Y = F + F*EPS(1)+EPS(2)
   IPRE=F
 
$THETA 0.412 0.871 5.32 1.92 0.00716 0.0106
$OMEGA BLOCK(3) 0.0753 0.0502 0.0909 0.0620 0.106 0.138
$SIGMA 0.0435 0.428

$SIMULATION (9215690) ONLYSIM SUBPROBLEM=1
$TABLE ID TIME AMT RATE WT AGE sex f_vis dose visit_n boot ONEHEADER FILE=./1.TAB NOPRINT


    
 


