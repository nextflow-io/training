package training.plugin

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

@CompileStatic
class GreetingExtension extends PluginExtensionPoint {

    @Override
    protected void init(Session session) {
    }

    @Function
    String reverseGreeting(String greeting) {
        return greeting.reverse()
    }

    @Function
    String decorateGreeting(String greeting) {
        return "*** ${greeting} ***"
    }

    @Function
    String friendlyGreeting(String greeting, String name = 'World') {
        return "${greeting}, ${name}!"
    }
}
