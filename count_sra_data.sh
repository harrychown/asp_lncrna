source activate lncRNA_env

GTF=$(echo /home/harry/Documents/lncRNA_August_21/ST_v9.1_merged_transcript_updated.gtf)

cd /media/harry/Elements/MSc/ox_iron/alignments
mkdir -f ../htseq_merged/
for file in *_sorted.bam
do
PREFIX=$( echo $file | sed 's/_sorted.bam//' )
htseq-count -s no -m intersection-strict -i transcript_id --nonunique all -f bam $file $GTF 1>../htseq_merged/"$PREFIX".counts
done

cd /media/harry/Elements/MSc/time_series/alignments
mkdir -f ../htseq_merged/
for file in *_sorted.bam
do
PREFIX=$( echo $file | sed 's/_sorted.bam//' )
PREFIX=$( echo $PREFIX | cut -d"_" -f2 )
htseq-count -s no -m intersection-strict -i transcript_id --nonunique all -f bam $file $GTF 1>../htseq_merged/"$PREFIX".counts
done

