package nextflow.greeting

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

/**
 * Plugin extension providing custom functions for greeting manipulation
 *
 * Functions defined here can be imported into Nextflow workflows using:
 *   include { functionName } from 'plugin/nf-greeting'
 */
@CompileStatic
class GreetingExtension extends PluginExtensionPoint {

    /**
     * Initialize the extension with the Nextflow session
     */
    @Override
    void init(Session session) {
        // Initialization logic if needed
    }

    /**
     * Reverse a greeting string
     *
     * @param greeting The original greeting
     * @return The reversed greeting
     */
    @Function
    String reverseGreeting(String greeting) {
        return greeting.reverse()
    }

    /**
     * Decorate a greeting with celebratory markers
     *
     * @param greeting The original greeting
     * @return The decorated greeting
     */
    @Function
    String decorateGreeting(String greeting) {
        return "*** ${greeting} ***"
    }

    /**
     * Convert greeting to a friendly format with exclamation
     *
     * @param greeting The original greeting
     * @param name Optional name to greet
     * @return The friendly greeting
     */
    @Function
    String friendlyGreeting(String greeting, String name = 'World') {
        return "${greeting}, ${name}!"
    }
}
