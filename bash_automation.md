# Automating workflows using bash

## Objectives
+ Write a shell script with multiple variables.
+ Incorporate a for loop into a shell script.
+ Understand how to make a workflow more efficient and less error-prone


We have now worked through two workflows, RNA-seq analysis and variant calling.
In both cases, we pasted all of our commands into the terminal, and used for 
loops to run each step of the workflow over our six samples. This is a great way
to understand how each tool works and to establish a workflow you want to use.

Let's imagine now that we decided to download all 672 samples from this dataset
instead of only working with our six samples. We already know how each tool in
our workflow works, and we know how to run them. We can now use bash scripting
to automate our workflow. 

We are going to automate the quality control portion of our workflow.

##  
