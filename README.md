# G74-2017-MESTAR
Stephanie's Project

Data for this project have processed over multipe courses. The first analysis was performed on whole-exome samples from Takara-medexome with IDT UMI's. 

The tools used/tested in this project are:
SNV/Indel calling: 
- Mutect2 
- Varscan2 (though later dropped due to high noise ratio which were impractical to filter)

CNV calling:
- Sequenza
- TitanCNA
- hatched (only in the testing phase and never run so far - interesting to use together with decipher for CCF (and DCF) calculation)

Clonality analysis:
- pyclone
- sciclone (though never finished)
