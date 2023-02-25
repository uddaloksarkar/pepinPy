#!/bin/bash

#filespos="satcomp17-nolimits"
module unload gcc/4.9.3
module load gcc/5.1.0

filespos="non-proj-bnn"
#filespos="delphinium-stp-cnf"

ulimit -t unlimited
shopt -s nullglob
rm -f todo
touch todo

test="pbcount"

solver_arr=(
"cryptominisat5"
)

solver_opt_arr=(
"cryptominisat5 --branchstr=vsids1 "
"cryptominisat5 --branchstr=vsids2 "
"cryptominisat5 --branchstr=maple1 "
)

output="pbcount"
tlimit="100"
#5GB mem limit
memlimit="5000000"
numthreads=$((OMPI_COMM_WORLD_SIZE))

SERVER=$PBS_O_HOST
WORKDIR="scratch/${PBS_JOBID}_${OMPI_COMM_WORLD_RANK}"
output="${output}-${PBS_JOBID}"

mkdir -p "${WORKDIR}"
cd "${WORKDIR}" || exit

files=$(ls ${PBS_O_WORKDIR}/${filespos}/*.cnf.xz | shuf --random-source=${PBS_O_WORKDIR}/myrnd)
outputdir="${PBS_O_WORKDIR}/${test}-main"

for solver in "${solver_arr[@]}"
do
    echo "${solver}"
    cp ${PBS_O_WORKDIR}/${solver} .
done
# create todo
rm -f todo
mkdir -p ${output}
numlines=0
at_opt=0

for solver_opt in "${solver_opt_arr[@]}"
do
    echo "${solver_opt}"
    for file in $files
    do
        mkdir -p "${output}-${at_opt}" || exit
        filename=$(basename "$file")
        filenameunzipped=${filename%.xz}

        # copy CNF and decompress
        echo "cp ${PBS_O_WORKDIR}/${filespos}/${filename} ." >> todo
        echo "unxz ${filename}" >> todo

        # run
        baseout="${output}-${at_opt}/${filename}"
        echo "( ulimit -t $tlimit; /usr/bin/time --verbose -o ${baseout}.timeout ./${solver_opt} ${filenameunzipped} > ${baseout}.out 2>&1; )" >> todo
        echo "( ulimit -t $tlimit; /usr/bin/time --verbose -o ${baseout}.timeout ./${solver_opt} ${filenameunzipped} > ${baseout}.out 2>&1; )"

        echo "rm ${filenameunzipped}" >> todo

        #copy back result
        echo "mkdir -p  ${outputdir}/${output}-${at_opt}" >> todo

        echo "xz ${baseout}.out*" >> todo
        echo "xz ${baseout}.timeout*" >> todo
        echo "rm core.*" >> todo

        echo "mv ${baseout}.out*      ${outputdir}/${output}-${at_opt}/" >> todo
        echo "mv ${baseout}.timeout*  ${outputdir}/${output}-${at_opt}/" >> todo
        echo "mv core.* ${outputdir}/${output}/" >> todo

        # delete what's left
        echo "rm ${baseout}.timeout*" >> todo
        echo "rm ${baseout}.out*" >> todo

        # todos: 3+4+4+4+4+1+2=22

        numlines=$((numlines+1))
    done
    at_opt=$((at_opt+1))
done
todomylines=13

# create per-core todos
numper=$((numlines/numthreads))
remain=$((numlines-numper*numthreads))
if [[ $remain -ge 1 ]]; then
    numper=$((numper+1))
fi
remain=$((numlines-numper*(numthreads-1)))

mystart=0
for ((myi=0; myi < numthreads ; myi++))
do
    rm -f todo_$myi.sh
    touch todo_$myi.sh
    echo "#!/bin/bash" > todo_$myi.sh
    echo "ulimit -v $memlimit" >> todo_$myi.sh
    echo "ulimit -c 0" >> todo_$myi.sh
    echo "set -x" >> todo_$myi.sh
    typeset -i myi
    typeset -i numper
    typeset -i mystart
    mystart=$((mystart + numper))
    if [[ $myi -lt $((numthreads-1)) ]]; then
        if [[ $mystart -gt $((numlines+numper)) ]]; then
            # echo "No need, over the limit by more than numper"
            sleep 0
        else
            if [[ $mystart -lt $numlines ]]; then
                myp=$((numper*todomylines))
                mys=$((mystart*todomylines))
                head -n $mys todo | tail -n $myp >> todo_$myi.sh
            else
                #we are at boundary, e.g. numlines is 100, numper is 3, mystart is 102
                #we must only print the last numper-(mystart-numlines) = 3-2 = 1
                mys=$((mystart*todomylines))
                p=$(( numper-mystart+numlines ))
                if [[ $p -gt 0 ]]; then
                    myp=$((p*todomylines))
                    head -n $mys todo | tail -n $myp >> todo_$myi.sh
                fi
            fi
        fi
    else
        if [[ $remain -gt 0 ]]; then
            mys=$((mystart*todomylines))
            mr=$((remain*todomylines))
            head -n $mys todo | tail -n $mr >> todo_$myi.sh
        fi
    fi
    echo "exit 0" >> todo_$myi.sh
    chmod +x todo_$myi.sh
done
# echo "Done."

# Execute todos
echo "This is MPI exec number $OMPI_COMM_WORLD_RANK"
rm -f ${output}/out_${OMPI_COMM_WORLD_RANK}
./todo_${OMPI_COMM_WORLD_RANK}.sh > ${output}/out_${OMPI_COMM_WORLD_RANK}
echo "Finished waiting rank $OMPI_COMM_WORLD_RANK"

rm -f ${solver}*
