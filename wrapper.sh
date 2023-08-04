#!/usr/bin/env bash

inputDir="/path/to/cell/files"
outputDir="/path/for/.SAM/files/output/folder"
sivDB="/path/to/siv/database"
mmulDB="/path/to/macaca/mulatta/database"
annotation="/path/to/annotation/file/.gtf"
threads=6
macacaSplicing="/path/to/splice/file/.ss"
trim_galore="/path/to/trim_galore"
fileoutputDir="/path/for/chimeric/sequences/Excel/output/folder"
mkdir -p ${outputDir}_R1
mkdir -p ${outputDir}_R1/tmp
mkdir -p ${outputDir}_R1/tmp/trim
mkdir -p ${outputDir}_R1/tmp/alns
mkdir -p ${outputDir}_R1/tmp/fq

mkdir -p ${outputDir}_R2
mkdir -p ${outputDir}_R2/tmp
mkdir -p ${outputDir}_R2/tmp/trim
mkdir -p ${outputDir}_R2/tmp/alns
mkdir -p ${outputDir}_R2/tmp/fq
mkdir -p ${outputDir}
mkdir -p ${outputDir}/tmp
mkdir -p ${outputDir}/alns
mkdir -p ${outputDir}/consensus

for file in "${inputDir}"/*R1*fastq.gz ; do
	sampleR1Base=$(basename ${file})
	sampleR1="${sampleR1Base%.*.*}"
	sample="${sampleR1%_R1*}"
	sampleR2Base="${sampleR1%_R1*}"_R2"${sampleR1##*_R1}".fastq.gz
	baseEnd="${sampleR1##*_R1}"
	mkdir -p {fileoutputDir}/${sample}

		echo ANALYZING ${sample}
		echo BEGIN DATE: ${date}
		TOTAL_TIME=0

		echo ADAPTOR AND QUALITY TRIMMING R1
		${trim_galore} --trim-n \
						-q 5 \
						--phred33 \
						--length 45 \
						-o ${outputDir}_R1/tmp/trim \
						--dont_gzip ${inputDir}/${sampleR1Base}
		trimmedFQ=${outputDir}_R1/tmp/trim/"${sampleR1%_R1*}"_R1"${sampleR1##*_R1}"_trimmed.fq
		mv ${trimmedFQ} ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq
		DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		echo DONE IN ${DUR}
		TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		echo ADAPTOR AND QUALITY TRIMMING R2
		${trim_galore} --trim-n \
						-q 5 \
						--phred33 \
						--length 45 \
						-o ${outputDir}_R2/tmp/trim \
						--dont_gzip ${inputDir}/${sampleR2Base}
		trimmedFQ=${outputDir}_R2/tmp/trim/"${sampleR1%_R1*}"_R2"${sampleR1##*_R1}"_trimmed.fq 
		mv ${trimmedFQ} ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq
		DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		echo DONE IN ${DUR}
		TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		SECONDS=0
		echo BOWTIE2 ALIGNING R1 TO SIV
		bowtie2 --very-sensitive-local \
				--no-unal \
				--local \
				--phred33 \
				-p ${threads} \
				-x ${sivDB} \
				-U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq \
				-S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.siv.bowtie.sam \
				--al ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq
		DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		echo DONE IN ${DUR}
		TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		SECONDS=0
		echo BOWTIE2 ALIGNING R2 TO SIV
		bowtie2 --very-sensitive-local \
				--no-unal \
				--local \
				--phred33 \
				-p ${threads} \
				-x ${sivDB} \
				-U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq \
				-S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.siv.bowtie.sam \
				--al ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq
		DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		echo DONE IN ${DUR}
		TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
		
		# extract read names into separate file for r1
		sed -n '1~4p' ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq > ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.reads
		# extract read names into separate file for r2
		sed -n '1~4p' ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq > ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.reads

		# get list of uniq names
		cat ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.reads ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.reads | awk -F ' ' '{print $1}' | sort | uniq > ${outputDir}/${sample}${baseEnd}.reads

		# now use this list to extract data for the macaque alignments
		awk 'NR==FNR{a[$0];next} $1 in a {print; getline; print; getline; print; getline; print}' ${outputDir}/${sample}${baseEnd}.reads ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.fastq > ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq

		awk 'NR==FNR{a[$0];next} $1 in a {print; getline; print; getline; print; getline; print}' ${outputDir}/${sample}${baseEnd}.reads ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.fastq > ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq


		SECONDS=0
		echo BOWTIE2 ALIGNING R1 TO MACACA MULATTA
		bowtie2 --very-sensitive-local \
				--no-unal \
				--local \
				--phred33 \
				-p ${threads} \
				-x ${mmulDB} \
				-U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq \
				-S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.mmul.bowtie.sam
		DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		echo DONE IN ${DUR}
		TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		SECONDS=0
		echo HISAT ALIGNING R1 TO MACACA MULATTA
		hisat2 -p ${threads} \
				--very-sensitive \
				--end-to-end \
				--known-splicesite-infile ${macacaSplicing} \
				--no-unal \
				--phred33 \
				-x ${mmulDB} \
				-U ${outputDir}_R1/tmp/fq/${sample}${baseEnd}_R1.al.fastq \
				-S ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.mmul.hisat.sam
		DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		echo DONE IN ${DUR}
		TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		SECONDS=0
		echo BOWTIE2 ALIGNING R2 TO MACACA MULATTA
		bowtie2 --very-sensitive-local \
				--no-unal \
				--local \
				--phred33 \
				-p ${threads} \
				-x ${mmulDB} \
				-U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq \
				-S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.mmul.bowtie.sam
		DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		echo DONE IN ${DUR}
		TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		SECONDS=0
		echo HISAT ALIGNING R2 TO MACACA MULATTA
		hisat2 -p ${threads} \
				--very-sensitive \
				--end-to-end \
				--known-splicesite-infile ${macacaSplicing} \
				--no-unal \
				--phred33 \
				-x ${mmulDB} \
				-U ${outputDir}_R2/tmp/fq/${sample}${baseEnd}_R2.al.fastq \
				-S ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.mmul.hisat.sam
		DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		echo DONE IN ${DUR}
		TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))

		source activate my_project_env
		SECONDS=0
		echo CREATING REPORTS IN ${fileoutputDir}/${sample}${baseEnd}/${sample}
		./chimFinder.py --pathogenR1 ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.siv.bowtie.sam \
						--hostR1 ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.mmul.bowtie.sam \
						--pathogenR2 ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.siv.bowtie.sam \
						--hostR2 ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.mmul.bowtie.sam \
						--splicedR1 ${outputDir}_R1/tmp/alns/${sample}${baseEnd}.mmul.hisat.sam \
						--splicedR2 ${outputDir}_R2/tmp/alns/${sample}${baseEnd}.mmul.hisat.sam \
						-o ${fileoutputDir}/${sample}${baseEnd}/${sample} \
						-t 12 \
						--minLen 20 \
						--maxLenUnmapped 30 \
						-a ${annotation} \
						--overlap 5 \
						--gap 5 \
						--minEntropy 0.84 \
						--close 5 \
						--score 0.75 \
						-q
		DUR="$(($SECONDS / 60)) minutes and $(($SECONDS % 60)) seconds"
		echo DONE IN ${DUR}
		TOTAL_TIME=$((TOTAL_TIME + ${SECONDS}))
		SECONDS=0
		conda deactivate
done
