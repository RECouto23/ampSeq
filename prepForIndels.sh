fastQDirectory=$1
outputDirectory=$2
mkdir -p $outputDirectory'fastqcReports/'
cd $fastQDirectory
ls *.fastq* > sampleSheet.txt
ls *.fastq* | rev | cut -c 16- | rev | uniq > roots.txt
cd - 

#FastQC Check
for each in $(cat $fastQDirectory'sampleSheet.txt')
do
fastqc $fastQDirectory$each -o $outputDirectory'fastqcReports/'
done
echo '------------------------------------FASTQC DONE------------------------------------'


mkdir -p $outputDirectory'trimmed/'

#Trim sequencing adapters 
for each in $(cat $fastQDirectory'roots.txt')
do
trimmomatic PE $fastQDirectory$each'R1_001.fastq.gz' $fastQDirectory$each'R2_001.fastq.gz' -threads 4 -baseout 'PEtrimmed_'$each ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True LEADING:3
done 
echo '------------------------------------TRIM DONE------------------------------------'
mv PEtrimmed* $outputDirectory'trimmed/'

mkdir -p $outputDirectory'combined/'  
for each in $(cat $fastQDirectory'roots.txt')
do 
flash $outputDirectory'trimmed/PEtrimmed_'$each'_1P' $outputDirectory'trimmed/PEtrimmed_'$each'_2P' --max-overlap=120 --min-overlap=10 --output-prefix 'combined_'$each
done 
mv combined_* $outputDirectory'combined/'
echo '------------------------------------FLASH DONE------------------------------------'

cd $outputDirectory'combined'
ls *.extendedFrags.fastq* > combinedSamples.txt
cd - 
mkdir -p $outputDirectory'mapped/'
for each in $(cat $outputDirectory'combined/combinedSamples.txt')
do
name=${each/.*}
bwa mem -t 6 GRCr38.fa $outputDirectory'combined/'$each > $name'.sam'
mv $name'.sam' $outputDirectory'mapped/'
done 
echo '------------------------------------MAPPING D0NE (BWA)------------------------------------'

cd $outputDirectory'mapped'
ls *.sam > mappedSamples.txt
cd -
mkdir -p $outputDirectory'readyForIndels/'

for each in $(cat $outputDirectory'mapped/mappedSamples.txt')
do 
name=${each/.*}
samtools view -bS $outputDirectory'mapped/'$each | samtools sort > $outputDirectory/'readyForIndels/'$name'.bam'
samtools index $outputDirectory'readyForIndels/'$name'.bam'
done 
echo '------------------------------------BAMs GENERATED------------------------------------'


find $outputDirectory'/readyForIndels/' -type f -name '*.bam' -not -name '*.bai' | cut -d '/' -f 6 > $outputDirectory'readyForIndels/bampleSheet.txt'

mkdir -p $outputDirectory'finalOutput'
for each in $(cat $outputDirectory'readyForIndels/bampleSheet.txt')
do
name=${each/.*}
echo $name
python3 count_indels_integrations/count_indels_integrations.py --bed guide.txt --ref GRCr38.fa --bam $outputDirectory'readyForIndels/'$each --out $outputDirectory'finalOutput/INDELS_'$name'.tsv'
done 
echo '------------------------------------COMPLETE------------------------------------'

cd $outputDirectory'finalOutput/'
mkdir -p logfiles/
mv *log* logfiles/
cat *.tsv > allFiles.tsv

