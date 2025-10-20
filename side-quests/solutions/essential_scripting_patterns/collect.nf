def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - groups multiple channel emissions into one
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }

// List.collect() - transforms each element, preserves structure
def formatted_ids = sample_ids.collect { id ->
    id.toUpperCase().replace('SAMPLE_', 'SPECIMEN_')
}
println "List.collect() result: ${formatted_ids} (${sample_ids.size()} items transformed into ${formatted_ids.size()})"

// Spread operator - concise property access
def sample_data = [[id: 's1', quality: 38.5], [id: 's2', quality: 42.1], [id: 's3', quality: 35.2]]
def all_ids = sample_data*.id
println "Spread operator result: ${all_ids}"
