nextflow.enable.types = true

record Person {
    id: String
    name: String
    greeting: String
}

record Greeting {
    id: String
    name: String
    greeting: String
    greeting_file: Path
}

record Shouted {
    id: String
    shouted: Path
}
