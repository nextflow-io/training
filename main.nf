include { GREETING_WORKFLOW } from './workflows/greeting'

workflow {
    channel.from('Alice', 'Bob', 'Charlie') | GREETING_WORKFLOW

    GREETING_WORKFLOW.out.greetings_ch.view { "Original: $it" }
    GREETING_WORKFLOW.out.timestamped_ch.view { "Timestamped: $it" }
}
