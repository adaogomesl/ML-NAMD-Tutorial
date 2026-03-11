NequIP-NAC model are trained and predicted in eV and eV/Å units.


  -------------------------------------------------------------------
                              /\
   |\    /|                  /++\
   ||\  /||                 /++++\
   || \/ || ||             /++++++\
   ||    || ||            /PyRAI2MD\
   ||    || ||           /++++++++++\                    __
            ||          /++++++++++++\    |\ |  /\  |\/| | \
            ||__ __    *==============*   | \| /--\ |  | |_/

                          Python Rapid
                     Artificial Intelligence
                  Ab Initio Molecular Dynamics



                      Author @Jingbai Li
                      
                      
    2022 – present   Hoffmann Institute of Advanced Materials
                     Shenzhen Polytechnic, Shenzhen, China    
                                
    2019 – 2022      Department of Chemistry and Chemical Biology
                     Northeastern University, Boston, USA

                          version:   2.5.1

  With contributions from (in alphabetic order):
    Jingbai Li     - Fewest switches surface hopping
                     Zhu-Nakamura surface hopping
                     Velocity Verlet
                     OpenMolcas interface
                     OpenMolcas/Tinker interface
                     BAGEL interface
                     ORCA interface
                     OpenQP interface
                     GFN-xTB interface
                     Adaptive sampling
                     Grid search
                     Multilayer ONIOM (QM:QM' and QM:QM':MM)
                     Periodic boundary condition
                     Wall potential
                     QC/ML hybrid NAMD
                     Excited-State Neural Network Potentials (ESNNP)

    Patrick Reiser - Fully connected neural networks (pyNNsMD)
                     SchNet (pyNNsMD)

  Special acknowledgement to:
    Steven A. Lopez            - Project co-founder
    Pascal Friederich          - Project co-founder




  &control
-------------------------------------------------------
  Title:                      dbh-12    
  NPROCS for ML:              1         
  NPROCS for QC:              1         
  NPROCS for Multiscale:      1          
  Seed:                       1         
  Job: 	                      md        
  QM:          	       	      nequip    
  Ab initio:                  molcas    
-------------------------------------------------------

        
  &molecule
-------------------------------------------------------
  States:                     [2]       
  Spin:                       [0]       
  Interstates:                [[1, 2]]  
  QMMM keyfile:               None      
  QMMM xyzfile:               Input     
  High level region:          []         ...
  Middel level region:        []         ...
  Boundary:                   []         ...
  Embedding charges:          False     
  Reading charge:             False     
  Frozen atoms:               []        
  Constrained atoms:          []        
  Restrained bonds:           []        
  Restrained angles:          []        
  Restrained dihedrals:       []        
  Restrained bond value:      []        
  Restrained angle value:     []        
  Restrained dihedral value:  []        
  Bond potential scale:       10.0      
  Angle potential scale:      0.005     
  Dihedral potential scale:   1e-06     
  External potential shape:   ['ellipsoid']
  External potential radius:  []        
  External potential factor:  [40]      
  External potential scale:   [1.0]       
  Compress potential volume:  []        
  Potential center:           []        
  Potential center type:      xyz       
  Constrained groups          []        
  Track geometry type:        None      
  Track indices:              []        
  Track threshold:            []        
  Track stop task:            0         
  Lattice vectors:            []        
  Lattice constant:           []        
  Periodic boundary condition []        
-------------------------------------------------------

        
  &initial condition
-------------------------------------------------------
  Generate initial condition: 0         
  Number:                     5000      
  Method:                     wigner     
  Format:                     molden    
-------------------------------------------------------

        
  &md
-------------------------------------------------------
  Initial state:              2         
  Initialize random velocity  0         
  Temperature (K):            300       
  Step:                       2000      
  Dt (au):                    20.67     
  Only active state grad      0         
  Surface hopping:            fssh      
  NAC type:                   nac       
  Phase correction            1         
  Population reset            0         
  Substep:                    20        
  Integrate probability       0         
  Decoherence:                0.1       
  Adjust velocity:            2         
  Reflect velocity:           1         
  Maxhop:                     10        
  IC hopping gap threshold    0.5       
  ISC hopping gap threshold   0.5       
  Thermodynamic:              off       
  Thermodynamic delay:        200       
  Print level:                0         
  Direct output:              2000      
  Buffer output:              0         
  Record MD data:             1         
  Record MD steps:            0         
  Checkpoint steps:           1          
  Restart function:           0         
  Additional steps:           0         
-------------------------------------------------------

        
  &md velocity control
-------------------------------------------------------
  Excess kinetic energy       0         
  Scale kinetic energy        1         
  Target kinetic energy       0         
  Gradient descent path       0         
  Reset velocity:             0         
  Reset step:                 0         
-------------------------------------------------------

        
  &nequip
-------------------------------------------------------
  Train data:                 None      
  Prediction data:            None      
  Train mode:                 training  
  Silent mode:                1         
  Data splits:                10        
  NN EG type:                 1         
  NN NAC type:                0         
  NN SOC type:                0         
  Multiscale:                 []        
  Shuffle data:               False     
  EG unit:                    si        
  NAC unit:                   si        
  SOC unit:                   cm-1
  Select EG output:           []        
  Select NAC output:          []        
  Select SOC output:          []          
  Data permutation            No        
  Elements                              
-------------------------------------------------------

        

  NequIP-NAC (Energy + Gradient + NAC)

  &hyperparameters            Energy+Gradient      Nonadiabatic         Spin-orbit
----------------------------------------------------------------------------------------------
  Model type:                 NequIPeg            NequIPnac           n/a          
  Model path:                 dbh-12.nequip.pth    dbh-12.nequip.pth    n/a                 
  GPU:                        1                    1                    n/a                 
  Chemical symbols:           auto                 auto                 n/a                 
----------------------------------------------------------------------------------------------
        


 *---------------------------------------------------*
 |                                                   |
 |             NequIP-NAC model only                 |
 |                                                   |
 *---------------------------------------------------*

 Number of atoms:  15
 Number of state:  [2]
 Number of NAC:    1
 Number of SOC:    None (not supported yet)

 Device found: cpu
 Running device: cpu
 

