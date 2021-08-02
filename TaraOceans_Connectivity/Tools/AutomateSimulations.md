
<h2>Automate Simulation runs in remote server</h2>

<h5>This is a guide to set up automated simulation runs of a python code on a remote server. This will allow one to run qsub jobs in series.</h5>

<h3>Step 1: Create a new executable schell script</h3>

<h5>write an executable shell script (here called Run_loop.sh) over "experiment.sh" file to pass arguments. Eg below: run a simulation for each month from 2009-2012 passing the year, month and day as parameters to the other shell script.</h5>
   
    #!/bin/bash
    startday=15
    for yr in {2009..2012}
    do
        for mon in {1..12}
        do
            # Job name is set here- allows configuration of the name for each job

            qsub -V -N FullTara_Res5_TS_${startday}${mon}${yr}_dt600 experiment.sh $yr $mon $startday

            # sleep time should be set based on the average run time of your simulation  
            # avoid running too many simulations on the shared gemini server 
            sleep 70m     
        done    
    done
<h3> Step 2: Pass arguments received in experiment.sh file to the Python simulation</h3>

    echo 'Loading parcels ...'
    echo 'Running simulation....'
    echo 'year: '$1, 'month: '$2, 'day: '$3
    cd ${HOME}/<projectfolder>/
    python3 <program>.py $1 $2 $3
    echo 'Finished computation.'
    
> Note: remove the Jobname from the experiment.sh file as it has already been assigned in the qsub call

<h3> Step 3: setup to be completed before executing the run</h3>

<h5>Make the shell script as executable. Run the following in your terminal</h5>

    chmod +x Run_loop.sh
<h5>Run the script as a background process from the server terminal</h5>

    nohup ./Run_loop.sh > outputfile.out &
<h5>in case you encounter a problem, get the process ID and kill the process</h5>

    top -u <userid>
    kill <PID>