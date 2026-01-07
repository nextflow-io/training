package nextflow.greeting

import groovy.transform.CompileStatic
import groovy.util.logging.Slf4j
import nextflow.plugin.BasePlugin
import org.pf4j.PluginWrapper

/**
 * Main plugin class for nf-greeting
 *
 * This class is the entry point for the plugin and is responsible
 * for plugin lifecycle management.
 */
@Slf4j
@CompileStatic
class GreetingPlugin extends BasePlugin {

    GreetingPlugin(PluginWrapper wrapper) {
        super(wrapper)
    }

    @Override
    void start() {
        super.start()
        log.info "nf-greeting plugin started"
    }

    @Override
    void stop() {
        super.stop()
        log.info "nf-greeting plugin stopped"
    }
}
