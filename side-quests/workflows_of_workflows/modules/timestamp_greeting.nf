process TIMESTAMP_GREETING {
    tag "adding timestamp to greeting"

    input:
    path greeting_file

    output:
    path 'timestamped_*.txt'

    script:
    def base_name = greeting_file.baseName
    """
    echo "[\$(date '+%Y-%m-%d %H:%M:%S')] \$(cat $greeting_file)" > timestamped_${base_name}.txt
    """
}
