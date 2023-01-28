<div class="formalpara-title">

**Solution**

</div>

``` nextflow
reads_ch        =  Channel.fromFilePairs(params.reads) 
GATK            =  params.gatk 
```

- Create a channel using [fromFilePairs()](https://www.nextflow.io/docs/latest/channel.html#fromfilepairs).

- A variable representing the path of GATK application file.
