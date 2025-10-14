def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// Channel.collect() - groups multiple channel emissions into one
ch_input = Channel.fromList(sample_ids)
ch_input.view { "Individual channel item: ${it}" }
ch_collected = ch_input.collect()
ch_collected.view { "Channel.collect() result: ${it} (${it.size()} items grouped into 1)" }

// Iterable.collect() - transforms each element, preserves structure
def formatted_ids = sample_ids.collect { id ->
    id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
}
println "Iterable.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

// Spread operator - concise property access
def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
def all_ids = sample_data*.id
println "Spread operator result: ${all_ids}"
