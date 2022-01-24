sample=$1      # sample name, for example CEN
thread=$2      # number of threads to run the pipeline
input_type=$3  # fastq or fasta
type_number=5  # analyze reads up to four rearrangements in a single read

# bash ../../scripts/run_stage_all.sh $sample $thread $input_type
# bash ../../scripts/run_stage_1.sh $sample $thread $input_type
# bash ../../scripts/run_stage_2.sh $sample $thread $input_type
# bash ../../scripts/run_stage_3.sh $sample $thread $input_type
# bash ../../scripts/run_stage_4.sh $sample $thread $input_type
# bash ../../scripts/run_stage_5.sh $sample $thread $input_type
bash ../../scripts/run_stage_test.sh $sample $thread $input_type