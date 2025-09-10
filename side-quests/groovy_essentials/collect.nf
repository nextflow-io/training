// Demonstrate Groovy vs Nextflow collect
def sample_ids = ['sample_001', 'sample_002', 'sample_003']

println "=== GROOVY COLLECT (transforms each item, keeps same structure) ==="
// Groovy collect: transforms each element but maintains list structure
def formatted_ids = sample_ids.collect { id ->
    id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
}
println "Original list: ${sample_ids}"
println "Groovy collect result: ${formatted_ids}"
println "Groovy collect maintains structure: ${formatted_ids.size} items (same as original)"
println ""

println "\n=== NEXTFLOW COLLECT (groups multiple items into single emission) ==="
// Nextflow collect: groups channel elements into a single emission
ch_input = Channel.of('sample_001', 'sample_002', 'sample_003')

// Show individual items before collect
ch_input.view { "Individual channel item: ${it}" }

// Collect groups all items into a single emission
ch_collected = ch_input.collect()
ch_collected.view { "Nextflow collect result: ${it} (${it.size()} items grouped together)" }
