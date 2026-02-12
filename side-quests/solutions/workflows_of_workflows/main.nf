include { GREETING_WORKFLOW } from './workflows/greeting'
include { TRANSFORM_WORKFLOW } from './workflows/transform'

workflow {
    main:
    names = channel.of('Alice', 'Bob', 'Charlie')

    // Run the greeting workflow
    GREETING_WORKFLOW(names)

    // Run the transform workflow
    TRANSFORM_WORKFLOW(GREETING_WORKFLOW.out.timestamped)

    // View results
    TRANSFORM_WORKFLOW.out.upper.view { "Uppercase: ${it}" }
    TRANSFORM_WORKFLOW.out.reversed.view { "Reversed: ${it}" }

    publish:
    greetings = GREETING_WORKFLOW.out.greetings
    upper = TRANSFORM_WORKFLOW.out.upper
    reversed = TRANSFORM_WORKFLOW.out.reversed
}

output {
    greetings {
    }
    upper {
    }
    reversed {
    }
}
