package training.plugin

import groovy.transform.CompileStatic
import nextflow.Session
import nextflow.plugin.extension.Function
import nextflow.plugin.extension.PluginExtensionPoint

@CompileStatic
class GreetingExtension extends PluginExtensionPoint {

    private String prefix = '***'
    private String suffix = '***'

    @Override
    protected void init(Session session) {
        // Read configuration with defaults
        // (GreetingConfig class provides validation/documentation)
        prefix = session.config.navigate('greeting.prefix', '***') as String
        suffix = session.config.navigate('greeting.suffix', '***') as String
    }

    @Function
    String reverseGreeting(String greeting) {
        return greeting.reverse()
    }

    @Function
    String decorateGreeting(String greeting) {
        return "${prefix} ${greeting} ${suffix}"
    }

    @Function
    String friendlyGreeting(String greeting, String name = 'World') {
        return "${greeting}, ${name}!"
    }
}
