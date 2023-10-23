# Interpretation


## Overview


<figure class="excalidraw">
--8<-- "docs/variantcalling/img/interpretation.excalidraw.svg"
</figure>



## Finding Causative Variants


```bash
cd /workspace/gitpod/nf-training/variantcalling/annotation/haplotypecaller/joint_variant_calling
```

Samples in VCF

```bash
zcat joint_germline_recalibrated_snpEff.ann.vcf.gz | grep "#CHROM" | cut -f 10-
```

which returns:

```bash
case_case       control_control
```

showing we case variants in field 10th and control variants in field 11th



```bash
zcat joint_germline_recalibrated_snpEff.ann.vcf.gz | grep PASS | grep HIGH | perl -nae 'if($F[10]=~/0\/0/ && $F[9]=~/1\/1/){print $_;}'
```

which results in 

```bash
chr21   32576780        rs541034925     A       AC      218.41  PASS    AC=2;AF=0.5;AN=4;DB;DP=90;ExcessHet=0;FS=0;MLEAC=2;MLEAF=0.5;MQ=60;POSITIVE_TRAIN_SITE;QD=31.2;SOR=4.174;VQSLOD=423864.23;culprit=FS;ANN=AC|frameshift_variant|HIGH|TCP10L|ENSG00000242220|transcript|ENST00000300258.8|protein_coding|5/5|c.641dupG|p.Val215fs|745/3805|641/648|214/215||,AC|frameshift_variant|HIGH|CFAP298-TCP10L|ENSG00000265590|transcript|ENST00000673807.1|protein_coding|8/8|c.1163dupG|p.Val389fs|1785/4781|1163/1170|388/389||       GT:AD:DP:GQ:PL  1/1:0,7:7:21:234,21,0   0/0:81,0:81:99:0,119,1600
```


[variant data](https://gnomad.broadinstitute.org/region/21-32576780-32576780?dataset=gnomad_r3)



```bash
zcat joint_germline_recalibrated_snpEff.ann.vcf.gz | grep PASS | grep HIGH | perl -nae 'if($F[10]=~/0\/0/ && $F[9]=~/0\/1/){print $_;}'
```

which results in:

```bash
chr21   44339194        rs769070783     T       C       45.91   PASS    AC=1;AF=0.25;AN=4;BaseQRankSum=-2.662;DB;DP=73;ExcessHet=0;FS=0;MLEAC=1;MLEAF=0.25;MQ=60;MQRankSum=0;POSITIVE_TRAIN_SITE;QD=3.53;ReadPosRankSum=1.6;SOR=1.402;VQSLOD=8.35;culprit=QD;ANN=C|start_lost|HIGH|CFAP410|ENSG00000160226|transcript|ENST00000397956.7|protein_coding|1/7|c.1A>G|p.Met1?|200/1634|1/1128|1/375||,C|upstream_gene_variant|MODIFIER|ENSG00000232969|ENSG00000232969|transcript|ENST00000426029.1|pseudogene||n.-182T>C|||||182|,C|downstream_gene_variant|MODIFIER|ENSG00000184441|ENSG00000184441|transcript|ENST00000448927.1|pseudogene||n.*3343T>C|||||3343|;LOF=(CFAP410|ENSG00000160226|1|1.00)        GT:AD:DP:GQ:PL  0/1:5,8:13:53:54,0,53   0/0:60,0:60:99:0,115,1002
chr21   44406660        rs139273180     C       T       44.91   PASS    AC=1;AF=0.25;AN=4;BaseQRankSum=-4.219;DB;DP=124;ExcessHet=0;FS=5.075;MLEAC=1;MLEAF=0.25;MQ=60;MQRankSum=0;POSITIVE_TRAIN_SITE;QD=0.67;ReadPosRankSum=0.681;SOR=0.774;VQSLOD=14.06;culprit=FS;ANN=T|stop_gained|HIGH|TRPM2|ENSG00000142185|transcript|ENST00000397932.6|protein_coding|19/33|c.2857C>T|p.Gln953*|2870/5216|2857/4662|953/1553||;LOF=(TRPM2|ENSG00000142185|1|1.00);NMD=(TRPM2|ENSG00000142185|1|1.00)   GT:AD:DP:GQ:PL  0/1:45,22:68:53:53,0,874        0/0:51,0:51:99:0,100,897
chr21   45505247        .       G       GGC     52.87   PASS    AC=1;AF=0.25;AN=4;BaseQRankSum=-0.524;DP=61;ExcessHet=0;FS=0;MLEAC=1;MLEAF=0.25;MQ=60;MQRankSum=0;POSITIVE_TRAIN_SITE;QD=13.22;ReadPosRankSum=-0.842;SOR=0.446;VQSLOD=77479.67;culprit=DP;ANN=GGC|frameshift_variant|HIGH|COL18A1|ENSG00000182871|transcript|ENST00000359759.8|protein_coding|34/41|c.4227_4228insGC|p.Pro1410fs|4228/6586|4228/5265|1410/1754||,GGC|intragenic_variant|MODIFIER|SLC19A1|ENSG00000173638|gene_variant|ENSG00000173638|||n.45505248_45505247insGC||||||;LOF=(COL18A1|ENSG00000182871|1|1.00)    GT:AD:DP:GQ:PL  0/1:3,1:5:25:61,0,25 0/0:54,0:54:99:0,102,1020
chr21   45989090        .       C       T       41.91   PASS    AC=1;AF=0.25;AN=4;BaseQRankSum=2.37;DP=86;ExcessHet=0;FS=0;MLEAC=1;MLEAF=0.25;MQ=60;MQRankSum=0;QD=2.99;ReadPosRankSum=-0.737;SOR=1.022;VQSLOD=9.09;culprit=QD;ANN=T|stop_gained|HIGH|COL6A1|ENSG00000142156|transcript|ENST00000361866.8|protein_coding|9/35|c.811C>T|p.Arg271*|892/4203|811/3087|271/1028||;LOF=(COL6A1|ENSG00000142156|1|1.00);NMD=(COL6A1|ENSG00000142156|1|1.00)     GT:AD:DP:GQ:PL       0/1:8,6:15:40:50,0,40   0/0:70,0:70:99:0,112,1494
```


Looking at [clinvar](https://www.ncbi.nlm.nih.gov/clinvar/)





## Conclusions



```bash
chr21   45989090        C       T       AC=1;AF=0.25;AN=4;BaseQRankSum=2.37;DP=86;ExcessHet=0;FS=0;MLEAC=1;MLEAF=0.25;MQ=60;MQRankSum=0;QD=2.99;ReadPosRankSum=-0.737;SOR=1.022;VQSLOD=9.09;culprit=QD;ANN=T|stop_gained|HIGH|COL6A1|ENSG00000142156|transcript|ENST00000361866.8|protein_coding|9/35|c.811C>T|p.Arg271*|892/4203|811/3087|271/1028||;LOF=(COL6A1|ENSG00000142156|1|1.00);NMD=(COL6A1|ENSG00000142156|1|1.00)     GT:AD:DP:GQ:PL       0/1:8,6:15:40:50,0,40   0/0:70,0:70:99:0,112,1494
```
