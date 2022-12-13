import org.yaml.snakeyaml.Yaml
import groovy.json.JsonSlurper


def parseJsonFile(path) {
    def slurper = new JsonSlurper()
    return slurper.parse(file(path))
}

def parseYamlFile(path) {
    def slurper = new Yaml()
    return slurper.load(file(path))
}
