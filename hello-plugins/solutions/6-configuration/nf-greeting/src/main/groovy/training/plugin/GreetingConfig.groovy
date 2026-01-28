package training.plugin

import nextflow.config.spec.ConfigOption
import nextflow.config.spec.ConfigScope
import nextflow.config.spec.ScopeName

/**
 * Configuration options for the nf-greeting plugin.
 *
 * Users configure these in nextflow.config:
 *
 *     greeting {
 *         enabled = true
 *         prefix = '>>>'
 *         suffix = '<<<'
 *         taskCounter.verbose = false
 *     }
 */
@ScopeName('greeting')
class GreetingConfig implements ConfigScope {

    GreetingConfig() {}

    GreetingConfig(Map opts) {
        this.enabled = opts.enabled as Boolean ?: true
        this.prefix = opts.prefix as String ?: '***'
        this.suffix = opts.suffix as String ?: '***'
        if (opts.taskCounter instanceof Map) {
            this.taskCounter = new TaskCounterConfig(opts.taskCounter as Map)
        }
    }

    /**
     * Enable or disable the plugin entirely.
     */
    @ConfigOption
    boolean enabled = true

    /**
     * Prefix for decorated greetings.
     */
    @ConfigOption
    String prefix = '***'

    /**
     * Suffix for decorated greetings.
     */
    @ConfigOption
    String suffix = '***'

    /**
     * Task counter configuration
     */
    TaskCounterConfig taskCounter = new TaskCounterConfig()

    static class TaskCounterConfig implements ConfigScope {
        TaskCounterConfig() {}
        TaskCounterConfig(Map opts) {
            this.verbose = opts.verbose as Boolean ?: true
        }

        @ConfigOption
        boolean verbose = true
    }
}
