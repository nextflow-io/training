def sample_ids = ['sample_001', 'sample_002', 'sample_003']

// channel.collect() - groups multiple channel emissions into one
ch_input = channel.fromList(sample_ids)
ch_input.view { sample -> "Individual channel item: ${sample}" }
ch_collected = ch_input.collect()
ch_collected.view { list -> "channel.collect() result: ${list} (${list.size()} items grouped into 1)" }
