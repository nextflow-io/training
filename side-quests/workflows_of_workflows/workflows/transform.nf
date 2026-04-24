include { SAY_HELLO_UPPER } from '../modules/say_hello_upper'
include { REVERSE_TEXT } from '../modules/reverse_text'

workflow {
    main:
    input_ch = channel.fromPath('results/timestamped/*.txt')

    // Apply transformations in sequence
    upper_ch = SAY_HELLO_UPPER(input_ch)
    reversed_ch = REVERSE_TEXT(upper_ch)

    publish:
    upper = upper_ch
    reversed = reversed_ch
}

output {
    upper {
    }
    reversed {
    }
}
