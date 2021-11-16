nums = Channel.from(1,2,3,4)         
square = nums.map { it -> it * it }  
square.view()