import groovy.json.JsonGenerator
import groovy.json.JsonGenerator.Converter

nextflow.enable.dsl=2

// comes from nf-test to store json files
params.nf_test_output  = ""

// include dependencies

include { SAMTOOLS_INDEX  } from '/workspaces/training/nf4-science/genomics/solutions/modules/gatk/haplotypecaller/tests/../../../samtools/index/main.nf'


// include test process
include { GATK_HAPLOTYPECALLER } from '/workspaces/training/nf4-science/genomics/solutions/modules/gatk/haplotypecaller/tests/../main.nf'

// define custom rules for JSON that will be generated.
def jsonOutput =
    new JsonGenerator.Options()
        .addConverter(Path) { value -> value.toAbsolutePath().toString() } // Custom converter for Path. Only filename
        .build()

def jsonWorkflowOutput = new JsonGenerator.Options().excludeNulls().build()


workflow {

    // run dependencies
    
    {
        def input = []
        
                    input[0] =  file("/workspaces/training/nf4-science/genomics/data/bam/reads_mother.bam")
                    
        SAMTOOLS_INDEX(*input)
    }
    

    // process mapping
    def input = []
    
                input[0] = SAMTOOLS_INDEX.out
                input[1] = file("/workspaces/training/nf4-science/genomics/data/ref/ref.fasta")
                input[2] = file("/workspaces/training/nf4-science/genomics/data/ref/ref.fasta.fai")
                input[3] = file("/workspaces/training/nf4-science/genomics/data/ref/ref.dict")
                input[4] = file("/workspaces/training/nf4-science/genomics/data/ref/intervals.bed")
                
    //----

    //run process
    GATK_HAPLOTYPECALLER(*input)

    if (GATK_HAPLOTYPECALLER.output){

        // consumes all named output channels and stores items in a json file
        for (def name in GATK_HAPLOTYPECALLER.out.getNames()) {
            serializeChannel(name, GATK_HAPLOTYPECALLER.out.getProperty(name), jsonOutput)
        }	  
      
        // consumes all unnamed output channels and stores items in a json file
        def array = GATK_HAPLOTYPECALLER.out as Object[]
        for (def i = 0; i < array.length ; i++) {
            serializeChannel(i, array[i], jsonOutput)
        }    	

    }
  
}

def serializeChannel(name, channel, jsonOutput) {
    def _name = name
    def list = [ ]
    channel.subscribe(
        onNext: {
            list.add(it)
        },
        onComplete: {
              def map = new HashMap()
              map[_name] = list
              def filename = "${params.nf_test_output}/output_${_name}.json"
              new File(filename).text = jsonOutput.toJson(map)		  		
        } 
    )
}


workflow.onComplete {

    def result = [
        success: workflow.success,
        exitStatus: workflow.exitStatus,
        errorMessage: workflow.errorMessage,
        errorReport: workflow.errorReport
    ]
    new File("${params.nf_test_output}/workflow.json").text = jsonWorkflowOutput.toJson(result)
    
}
