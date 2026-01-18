# Summary

Congratulations on completing the Hello Nextflow training course! 🎉

## Your journey

You started with [...].
Over the course of six parts, you've transformed that basic workflow into [...]].

### What you built

The final form of the Hello workflow takes as input a CSV file containing greetings.

The four steps are implemented as Nextflow processes (`sayHello`, `convertToUpper`, `collectGreetings`, and `cowpy`) stored in separate module files.

The results are published to a directory called `results/`, and the final output of the pipeline (when run with default parameters) is a plain text file containing ASCII art of a turkey saying the uppercased greetings.

<figure class="excalidraw">
--8<-- "docs/hello_nextflow/img/hello_pipeline_complete.svg"
</figure>

1. **`sayHello`:** Writes each greeting to its own output file (_e.g._ "Hello-output.txt")
2. **`convertToUpper`:** Converts each greeting to uppercase (_e.g._ "HELLO")
3. **`collectGreetings`:** Collects all uppercase greetings into a single batch file
4. **`cowpy`:** Generates ASCII art using the `cowpy` tool

[add note about configuration]

### Key skills acquired

Through this hands-on course, you've learned to:

[...]

## From research script to robust pipeline

[...]

You're now equipped with the foundational knowledge to [...].

Thank you for completing this training.
Please take a minute to complete the feedback survey on the next page.

**Happy pipelining!** 🚀
