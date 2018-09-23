for tau_mix in 0.01 0.05 0.1 0.3
do
    for tau_ent_main in 0.1 0.2 1 2 3
    do
        for tau_ent_sec in 0.05 0.1 0.2 0.5 1.0 2.0 3.0 
        do
            filename="tauMix_$tau_mix-entMain_$tau_ent_main-entSec_$tau_ent_sec.pbs"
            out_filename="tauMix_$tau_mix-entMain_$tau_ent_main-entSec_$tau_ent_sec.out"
            output_dir="/nv/hp16/egoh6/data/09-17-2018_FiniteMixing_and_Entrainment"
            echo "#PBS -N $filename" >> $filename
            echo "#PBS -l nodes=1:ppn=4" >> $filename
            echo "#PBS -l pmem=2gb" >> $filename
            echo "#PBS -l walltime=01:30:00" >> $filename
            echo "#PBS -q force-6" >> $filename
            echo "#PBS -j oe" >> $filename
            echo "#PBS -o $out_filename" >> $filename
            echo "#PBS -m ae" >> $filename
            echo "#PBS -M edwin.goh@gatech.edu" >> $filename
            echo "" >> $filename
            echo "source activate canteraEnv" >> $filename
            echo "cd /nv/hp16/egoh6/data/BatchPaSR/Scripts" >> $filename
            echo "time python FiniteMixing_and_FiniteEntrainment.py $tau_mix $tau_ent_main $tau_ent_sec $output_dir" >> $filename
            echo "echo \"Finished running $filename\""
        done
    done    
done