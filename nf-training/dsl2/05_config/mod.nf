process foo {
 label 'big_job'
 echo true
 """
 echo TASK="$task.process CPUs=$task.cpus MEM=$task.memory" QUEUE=$task.queue
 """
}

