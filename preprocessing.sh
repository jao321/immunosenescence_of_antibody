dir="path/to/folder/with_R1_and_R2-fastq_files/"


for entry in "$dir"/*
do

   pasta="$entry"


	foldername=${pasta##*/}
	#IFS='-'
	read -a sample <<< "$foldername"

	python3 /local_usr/AssemblePairs.py align -1 "$pasta"/R1.fastq -2 "$pasta"/R2.fastq --nproc 24 --rc tail --coord illumina --outdir "$pasta" --outname "${sample[0]}" --minlen 50 --log "${sample[0]}".log --failed

	python3 /local_usr/FilterSeq.py quality -s "$pasta"/"${sample[0]}"_assemble-pass.fastq -q 30 --outdir "$pasta" --outname "${sample[0]}" --log "${sample[0]}"_quality.log --failed

	python3 /local_usr/MaskPrimers.py score -s "$pasta"/"${sample[0]}"_quality-pass.fastq -p /primers_seq.fasta --outdir "$pasta" --outname "${sample[0]}" --mode tag --fasta --log "${sample[0]}"_primer.log --failed

	python3 /local_usr/AssignGenes.py igblast -s "$pasta"/"${sample[0]}"_primers-pass.fasta -b /local_usr/igblast --organism human --loci ig --format blast --exec /local_usr/igblastn --outdir "$pasta" --nproc 24

	python3 /local_usr/MakeDb.py igblast -s "$pasta"/"${sample[0]}"_primers-pass.fasta -i "$pasta"/"${sample[0]}"_primers-pass_igblast.fmt7 --format airr -r /imgt/human/vdj --outdir "$pasta" --outname "${sample[0]}" --extended 
	
	rm "$pasta"/"${sample[0]}"_primers-pass_igblast.fmt7




done