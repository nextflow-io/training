def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// Nextflow collect() - groups multiple channel emissions into one
ch_input = Channel.fromList(sample_ids)
ch_input.view { "Individual channel item: ${it}" }
ch_collected = ch_input.collect()
ch_collected.view { "Nextflow collect() result: ${it} (${it.size()} items grouped into 1)" }
